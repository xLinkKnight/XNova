//! The processors of this computer, and work spread over them in order.
//!
//! Windows puts at most 64 logical processors in a processor group, and the
//! usual counts (`GetSystemInfo`, which `num_cpus` uses, and Rust's
//! `available_parallelism`) see only the calling thread's group: 64 of a
//! 128-core, 256-thread EPYC. [`logical_cpus`] counts every group, and
//! [`spread_thread`] places a worker thread in a group, so that a pool of
//! more than 64 threads runs on all of them. Windows 11 and Windows Server
//! 2022 and later spread a process's threads over the groups themselves;
//! earlier versions keep them in one group unless told otherwise.

use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};

/// Logical processors (hardware threads) of this computer, every processor
/// group counted
pub fn logical_cpus() -> usize {
    let std = std::thread::available_parallelism().map_or(1, |n| n.get());
    #[cfg(windows)]
    {
        let all: usize = windows::active_groups().iter().map(|g| g.count).sum();
        all.max(std)
    }
    #[cfg(not(windows))]
    std
}

/// Processor groups of this computer (always 1 outside Windows)
pub fn processor_groups() -> usize {
    #[cfg(windows)]
    {
        windows::active_groups().len().max(1)
    }
    #[cfg(not(windows))]
    1
}

/// Put the calling thread, worker `index` of a pool of `threads`, in the
/// processor group its share of the processors falls in: the workers spread
/// evenly over every group. Returns the group, or `None` with a single group
/// (or when Windows refuses, e.g. in a job limited to some processors; the
/// thread then stays where it was).
pub fn spread_thread(index: usize, threads: usize) -> Option<u16> {
    #[cfg(windows)]
    {
        let groups = windows::active_groups();
        let counts: Vec<usize> = groups.iter().map(|g| g.count).collect();
        let group = &groups[group_of(index, threads, &counts)?];
        windows::join_group(group).then_some(group.number)
    }
    #[cfg(not(windows))]
    {
        let _ = (index, threads);
        None
    }
}

/// The group worker `index` of `threads` belongs in, the groups holding
/// `counts` processors: worker `i` takes the processor at `i * total /
/// threads` (every processor in turn when there are more workers than
/// processors). `None` with fewer than two groups.
fn group_of(index: usize, threads: usize, counts: &[usize]) -> Option<usize> {
    let total: usize = counts.iter().sum();
    if counts.len() < 2 || total == 0 || threads == 0 {
        return None;
    }
    let place = if threads >= total { index % total } else { index * total / threads };
    let mut start = 0;
    for (group, &count) in counts.iter().enumerate() {
        if place < start + count {
            return Some(group);
        }
        start += count;
    }
    Some(counts.len() - 1)
}

/// `f` of every item, returned in the items' order, run on the rayon pool
/// with the items taken in order, a few at a time: the first items are done
/// first, so the spectra of a batch's first file are all searched long before
/// those of its last one (a parallel iterator splits the items into ranges and
/// finishes every range near the end). Each thread takes the next few items
/// as soon as it is free, so the work stays as balanced as with a parallel
/// iterator.
pub fn par_map_in_order<T: Sync, R: Send>(items: &[T], f: impl Fn(&T) -> R + Sync) -> Vec<R> {
    let n = items.len();
    if n == 0 {
        return Vec::new();
    }
    let threads = rayon::current_num_threads().max(1);
    // a few items at a time: in order, without the threads queueing for fast items
    let step = (n / (threads * 64)).clamp(1, 16);
    let next = AtomicUsize::new(0);
    let parts: Vec<Vec<(usize, R)>> = (0..threads)
        .into_par_iter()
        .map(|_| {
            let mut out = Vec::new();
            loop {
                let start = next.fetch_add(step, Ordering::Relaxed);
                if start >= n {
                    break out;
                }
                for i in start..(start + step).min(n) {
                    out.push((i, f(&items[i])));
                }
            }
        })
        .collect();
    let mut slots: Vec<Option<R>> = std::iter::repeat_with(|| None).take(n).collect();
    for (i, r) in parts.into_iter().flatten() {
        slots[i] = Some(r);
    }
    slots.into_iter().map(|r| r.expect("every item is mapped once")).collect()
}

#[cfg(windows)]
mod windows {
    use std::ffi::c_void;

    /// A processor group: its number, its active processors and their mask
    pub struct Group {
        pub number: u16,
        pub count: usize,
        pub mask: usize,
    }

    /// GROUP_AFFINITY
    #[repr(C)]
    struct GroupAffinity {
        mask: usize,
        group: u16,
        reserved: [u16; 3],
    }

    /// LOGICAL_PROCESSOR_RELATIONSHIP RelationGroup
    const RELATION_GROUP: u32 = 4;

    #[link(name = "kernel32")]
    extern "system" {
        fn GetLogicalProcessorInformationEx(relationship: u32, buffer: *mut u8, length: *mut u32) -> i32;
        fn GetCurrentThread() -> *mut c_void;
        fn SetThreadGroupAffinity(
            thread: *mut c_void,
            affinity: *const GroupAffinity,
            previous: *mut GroupAffinity,
        ) -> i32;
    }

    /// The active processor groups, as Windows reports them
    pub fn active_groups() -> Vec<Group> {
        let mut len = 0u32;
        // Safety: asking for the size (null buffer); then a buffer of that size, 8-byte aligned
        unsafe { GetLogicalProcessorInformationEx(RELATION_GROUP, std::ptr::null_mut(), &mut len) };
        if len == 0 {
            return Vec::new();
        }
        let mut buffer = vec![0u64; (len as usize).div_ceil(8)];
        if unsafe { GetLogicalProcessorInformationEx(RELATION_GROUP, buffer.as_mut_ptr() as *mut u8, &mut len) } == 0
        {
            return Vec::new();
        }
        let bytes = unsafe { std::slice::from_raw_parts(buffer.as_ptr() as *const u8, len as usize) };
        parse_groups(bytes)
    }

    /// The groups of a `GetLogicalProcessorInformationEx(RelationGroup)`
    /// answer: SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX records (relationship
    /// u32, size u32, then GROUP_RELATIONSHIP: maximum and active group counts
    /// (u16 each), 20 reserved bytes, and per active group a 48-byte
    /// PROCESSOR_GROUP_INFO: maximum and active processor counts (u8 each), 38
    /// reserved bytes, the active processor mask (u64))
    pub fn parse_groups(bytes: &[u8]) -> Vec<Group> {
        let u16_at = |at: usize| u16::from_le_bytes([bytes[at], bytes[at + 1]]);
        let u32_at = |at: usize| u32::from_le_bytes(bytes[at..at + 4].try_into().expect("4 bytes"));
        let mut out = Vec::new();
        let mut at = 0;
        while at + 8 <= bytes.len() {
            let (relation, size) = (u32_at(at), u32_at(at + 4) as usize);
            if size == 0 || at + size > bytes.len() {
                break;
            }
            if relation == RELATION_GROUP && size >= 32 {
                let active = u16_at(at + 10) as usize;
                for group in 0..active {
                    let info = at + 32 + 48 * group;
                    if info + 48 > at + size {
                        break;
                    }
                    out.push(Group {
                        number: group as u16,
                        count: bytes[info + 1] as usize,
                        mask: u64::from_le_bytes(bytes[info + 40..info + 48].try_into().expect("8 bytes")) as usize,
                    });
                }
            }
            at += size;
        }
        out
    }

    /// Run the calling thread on `group`'s processors
    pub fn join_group(group: &Group) -> bool {
        let affinity = GroupAffinity { mask: group.mask, group: group.number, reserved: [0; 3] };
        // Safety: the pseudo-handle of the calling thread, and a valid GROUP_AFFINITY
        unsafe { SetThreadGroupAffinity(GetCurrentThread(), &affinity, std::ptr::null_mut()) != 0 }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    /// Workers spread evenly over the groups; one group: no placing at all
    #[test]
    fn threads_spread_over_processor_groups() {
        let four = [64, 64, 64, 64];
        // as many workers as processors: 64 to a group
        let groups: Vec<usize> = (0..256).map(|i| group_of(i, 256, &four).unwrap()).collect();
        assert!(groups[..64].iter().all(|&g| g == 0) && groups[192..].iter().all(|&g| g == 3));
        // fewer: 16 to a group, not 64 in the first
        let mut per_group = [0; 4];
        for i in 0..64 {
            per_group[group_of(i, 64, &four).unwrap()] += 1;
        }
        assert_eq!(per_group, [16; 4]);
        // more workers than processors: around again
        assert_eq!(group_of(256, 300, &four), Some(0));
        assert_eq!(group_of(299, 300, &four), Some(0));
        // groups of different sizes (e.g. 36 + 36 + 20 processors)
        let uneven = [36, 36, 20];
        let mut per_group = [0; 3];
        for i in 0..92 {
            per_group[group_of(i, 92, &uneven).unwrap()] += 1;
        }
        assert_eq!(per_group, [36, 36, 20]);
        assert_eq!(group_of(0, 64, &[64]), None);
        assert_eq!(group_of(0, 64, &[]), None);
    }

    #[cfg(windows)]
    #[test]
    fn processor_groups_are_parsed() {
        // two groups of 64 and 32 processors, as GetLogicalProcessorInformationEx answers
        let mut bytes = vec![0u8; 32 + 48 * 2];
        let size = bytes.len() as u32;
        bytes[0..4].copy_from_slice(&4u32.to_le_bytes());
        bytes[4..8].copy_from_slice(&size.to_le_bytes());
        bytes[8..10].copy_from_slice(&2u16.to_le_bytes());
        bytes[10..12].copy_from_slice(&2u16.to_le_bytes());
        for (g, count) in [(0usize, 64u8), (1, 32)] {
            let info = 32 + 48 * g;
            bytes[info] = count;
            bytes[info + 1] = count;
            let mask: u64 = if count == 64 { u64::MAX } else { (1u64 << count) - 1 };
            bytes[info + 40..info + 48].copy_from_slice(&mask.to_le_bytes());
        }
        let groups = windows::parse_groups(&bytes);
        assert_eq!(groups.len(), 2);
        assert_eq!((groups[0].number, groups[0].count, groups[0].mask), (0, 64, usize::MAX));
        assert_eq!((groups[1].number, groups[1].count, groups[1].mask), (1, 32, (1usize << 32) - 1));
        // this computer: at least what Rust counts, in at least one group
        assert!(logical_cpus() >= std::thread::available_parallelism().unwrap().get());
        assert!(processor_groups() >= 1);
        eprintln!("this computer: {} logical CPUs in {} group(s)", logical_cpus(), processor_groups());
    }

    /// Results come back in the items' order, each item mapped once, and the
    /// items are started roughly in order
    #[test]
    fn par_map_in_order_keeps_the_order() {
        let items: Vec<usize> = (0..10_000).collect();
        let started = AtomicUsize::new(0);
        let mapped: Vec<(usize, usize)> = par_map_in_order(&items, |&i| (i * 2, started.fetch_add(1, Ordering::Relaxed)));
        assert_eq!(mapped.iter().map(|m| m.0).collect::<Vec<_>>(), items.iter().map(|i| i * 2).collect::<Vec<_>>());
        // nearly every item starts within a few threads' worth of its place
        // (a thread the system pauses mid-batch starts its next items late)
        let bound = 16 * (rayon::current_num_threads() + 1);
        let late = mapped.iter().enumerate().filter(|(i, m)| m.1.abs_diff(*i) > bound).count();
        assert!(late < items.len() / 20, "{late} of {} items started out of order", items.len());
        assert!(par_map_in_order(&Vec::<u8>::new(), |&x| x).is_empty());
    }
}
