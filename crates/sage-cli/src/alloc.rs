//! The memory allocator of sage.exe, chosen at start-up (`--allocator`):
//! mimalloc (the default) or the system's.
//!
//! Reading Bruker timsTOF data decompresses frames on every thread at once,
//! and in Windows' allocator that slows down the longer a process runs: one
//! 1.9 GB file read in 2.0, 5.6, 9.0 and 11.3 s on 16 threads, one read after
//! another, and in 9.2 s every time on one thread (`read_timing` in
//! sage-cloudpath's `tdf.rs`).
//!
//! A Rust program has one global allocator, fixed when it is built. This one
//! hands each new block to the allocator chosen, and gives every block back to
//! the allocator it came from (mimalloc knows its own blocks), so the choice
//! can be made once the command line is read, after the first allocations.

use std::alloc::{GlobalAlloc, Layout, System};
use std::ffi::c_void;
use std::sync::atomic::{AtomicBool, Ordering};

/// The allocators to choose from
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Allocator {
    /// mimalloc (Microsoft Research, MIT licence)
    Mimalloc,
    /// The operating system's (on Windows, the process heap)
    System,
}

impl Allocator {
    /// Its name on the command line
    pub fn name(self) -> &'static str {
        match self {
            Allocator::Mimalloc => "mimalloc",
            Allocator::System => "system",
        }
    }

    pub fn parse(name: &str) -> Option<Self> {
        match name {
            "mimalloc" => Some(Allocator::Mimalloc),
            "system" => Some(Allocator::System),
            _ => None,
        }
    }
}

/// Whether new blocks come from mimalloc
static MIMALLOC: AtomicBool = AtomicBool::new(true);
/// Whether any block has come from mimalloc: until then every block is the
/// system's, and freeing one needs no look-up
static USED_MIMALLOC: AtomicBool = AtomicBool::new(false);

/// Take new blocks from `allocator` from now on
pub fn select(allocator: Allocator) {
    MIMALLOC.store(allocator == Allocator::Mimalloc, Ordering::Relaxed);
}

/// The allocator new blocks come from
pub fn selected() -> Allocator {
    if MIMALLOC.load(Ordering::Relaxed) {
        Allocator::Mimalloc
    } else {
        Allocator::System
    }
}

/// The global allocator of sage.exe (`#[global_allocator]` in `main.rs`)
pub struct Switchable;

impl Switchable {
    #[inline]
    fn take_mimalloc() -> bool {
        let on = MIMALLOC.load(Ordering::Relaxed);
        // set once: a shared flag that every thread wrote would slow them all
        if on && !USED_MIMALLOC.load(Ordering::Relaxed) {
            USED_MIMALLOC.store(true, Ordering::Relaxed);
        }
        on
    }

    /// Whether mimalloc gave out `ptr`
    #[inline]
    unsafe fn from_mimalloc(ptr: *mut u8) -> bool {
        USED_MIMALLOC.load(Ordering::Relaxed) && libmimalloc_sys::mi_is_in_heap_region(ptr as *const c_void)
    }
}

unsafe impl GlobalAlloc for Switchable {
    #[inline]
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        if Self::take_mimalloc() {
            libmimalloc_sys::mi_malloc_aligned(layout.size(), layout.align()) as *mut u8
        } else {
            System.alloc(layout)
        }
    }

    #[inline]
    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        if Self::take_mimalloc() {
            libmimalloc_sys::mi_zalloc_aligned(layout.size(), layout.align()) as *mut u8
        } else {
            System.alloc_zeroed(layout)
        }
    }

    #[inline]
    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        if Self::from_mimalloc(ptr) {
            libmimalloc_sys::mi_free(ptr as *mut c_void)
        } else {
            System.dealloc(ptr, layout)
        }
    }

    #[inline]
    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        // a block grows in the allocator it came from
        if Self::from_mimalloc(ptr) {
            libmimalloc_sys::mi_realloc_aligned(ptr as *mut c_void, new_size, layout.align()) as *mut u8
        } else {
            System.realloc(ptr, layout, new_size)
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    /// Blocks taken before and after a change of allocator are grown and
    /// freed by the allocator they came from, whatever is chosen by then
    #[test]
    fn blocks_go_back_to_their_allocator() {
        let a = Switchable;
        let layout = Layout::from_size_align(1000, 8).unwrap();
        let wide = Layout::from_size_align(4096, 64).unwrap();
        unsafe {
            select(Allocator::System);
            let system = a.alloc(layout);
            assert!(!system.is_null());
            select(Allocator::Mimalloc);
            let mi = a.alloc_zeroed(wide);
            assert!(!mi.is_null() && mi as usize % 64 == 0);
            assert!(std::slice::from_raw_parts(mi, 4096).iter().all(|&b| b == 0));
            assert!(libmimalloc_sys::mi_is_in_heap_region(mi as *const c_void));
            assert!(!libmimalloc_sys::mi_is_in_heap_region(system as *const c_void));

            // the system's block grows in the system's allocator, though mimalloc is chosen
            system.write_bytes(7, 1000);
            let system = a.realloc(system, layout, 100_000);
            assert!(!libmimalloc_sys::mi_is_in_heap_region(system as *const c_void));
            assert!(std::slice::from_raw_parts(system, 1000).iter().all(|&b| b == 7));
            // and mimalloc's in mimalloc, though the system's is chosen
            select(Allocator::System);
            let mi = a.realloc(mi, wide, 1 << 20);
            assert!(libmimalloc_sys::mi_is_in_heap_region(mi as *const c_void));
            assert_eq!(mi as usize % 64, 0);

            a.dealloc(system, Layout::from_size_align(100_000, 8).unwrap());
            a.dealloc(mi, Layout::from_size_align(1 << 20, 64).unwrap());
        }
        select(Allocator::Mimalloc);
        assert_eq!(selected(), Allocator::Mimalloc);
        assert_eq!(Allocator::parse("system"), Some(Allocator::System));
        assert_eq!(Allocator::parse("jemalloc"), None);
    }
}
