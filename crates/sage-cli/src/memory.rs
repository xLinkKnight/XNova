//! sage.exe's memory: a limit (`max_memory_gb` / `--max-memory`), the peak
//! use it reports, and the database size estimate (`--estimate`).
//!
//! A background thread samples the process's memory (its resident set, the
//! working set on Windows) every 250 ms; allocations themselves are not
//! instrumented. With a limit, reaching it stops the search with a message
//! and exit code 137, instead of the computer running out of memory.

use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::Arc;
use std::thread::JoinHandle;
use std::time::Duration;
use sysinfo::{ProcessExt, System, SystemExt};

/// Bytes in a GB, as sage-gui counts them (GiB)
pub const GIB: f64 = 1024.0 * 1024.0 * 1024.0;
const POLL: Duration = Duration::from_millis(250);
/// The exit code of a search stopped at its memory limit
pub const EXIT_CODE: i32 = 137;

pub fn gib(bytes: u64) -> f64 {
    bytes as f64 / GIB
}

/// Samples the process's memory until dropped
pub struct MemoryGuard {
    stop: Arc<AtomicBool>,
    peak: Arc<AtomicU64>,
    handle: Option<JoinHandle<()>>,
}

impl MemoryGuard {
    /// Start sampling; above `limit_gb` (GiB, if positive) the process stops
    pub fn start(limit_gb: Option<f64>) -> Self {
        let stop = Arc::new(AtomicBool::new(false));
        let peak = Arc::new(AtomicU64::new(0));
        let limit = limit_gb.filter(|gb| *gb > 0.0).map(|gb| (gb * GIB) as u64);
        let handle = sysinfo::get_current_pid().ok().and_then(|pid| {
            let (stop, peak) = (stop.clone(), peak.clone());
            std::thread::Builder::new()
                .name("memory-guard".into())
                .spawn(move || {
                    let mut system = System::new();
                    while !stop.load(Ordering::Relaxed) {
                        if system.refresh_process(pid) {
                            if let Some(process) = system.process(pid) {
                                let used = process.memory();
                                peak.fetch_max(used, Ordering::Relaxed);
                                if let Some(limit) = limit.filter(|&limit| used >= limit) {
                                    exceeded(used, limit);
                                }
                            }
                        }
                        std::thread::sleep(POLL);
                    }
                })
                .ok()
        });
        Self { stop, peak, handle }
    }

    /// The most memory seen so far (bytes)
    pub fn peak_bytes(&self) -> u64 {
        self.peak.load(Ordering::Relaxed)
    }
}

impl Drop for MemoryGuard {
    fn drop(&mut self) {
        self.stop.store(true, Ordering::Relaxed);
        if let Some(handle) = self.handle.take() {
            let _ = handle.join();
        }
    }
}

fn exceeded(used: u64, limit: u64) -> ! {
    log::error!(
        "sage.exe uses {:.2} GB of memory, its limit is {:.2} GB (max_memory_gb / --max-memory): the search is \
         stopped. Search fewer files at once (--batch-size), fewer variable modifications or missed cleavages, \
         turn on the database prefilter for a large FASTA, or raise the limit.",
        gib(used),
        gib(limit)
    );
    std::process::exit(EXIT_CODE)
}

/// Total size of the spectrum files (a Bruker .d folder: its files), bytes;
/// files that are not local count as 0
pub fn spectra_bytes(paths: &[url::Url]) -> u64 {
    fn size(path: &std::path::Path) -> u64 {
        match std::fs::metadata(path) {
            Ok(m) if m.is_dir() => std::fs::read_dir(path)
                .map(|entries| entries.filter_map(Result::ok).map(|e| size(&e.path())).sum())
                .unwrap_or(0),
            Ok(m) => m.len(),
            Err(_) => 0,
        }
    }
    paths
        .iter()
        .filter(|url| url.scheme() == "file")
        .filter_map(|url| url.to_file_path().ok())
        .map(|path| size(&path))
        .sum()
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn the_peak_is_sampled() {
        let guard = MemoryGuard::start(None);
        // something to see: 64 MB touched
        let block = vec![1u8; 64 << 20];
        std::thread::sleep(Duration::from_millis(600));
        assert!(block.iter().map(|&b| b as u64).sum::<u64>() > 0);
        assert!(guard.peak_bytes() > 64 << 20, "{}", guard.peak_bytes());
        drop(guard);
        assert!((gib(3 * (1 << 30)) - 3.0).abs() < 1e-12);
    }
}
