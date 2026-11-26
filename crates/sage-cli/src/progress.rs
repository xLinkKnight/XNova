//! What sage.exe is doing, in the log while a step runs long.
//!
//! A step (reading a batch of files, searching it, the FDR and reports)
//! registers itself with [`step`]; a background thread ([`Heartbeat`]) writes
//! a `progress:` line about it every [`INTERVAL`] once it has run that long,
//! so the log never goes quiet for more than that. Steps that finish sooner
//! leave no progress line, only their own result line.

use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread::JoinHandle;
use std::time::{Duration, Instant};

/// How often a long step is reported
pub const INTERVAL: Duration = Duration::from_secs(15);

struct Step {
    label: String,
    /// What the step is busy with now (e.g. the file being read)
    detail: Arc<Mutex<String>>,
    unit: &'static str,
    total: usize,
    done: Arc<AtomicUsize>,
    start: Instant,
    /// Distinguishes this step from the next one with the same label
    id: usize,
}

static STEP: Mutex<Option<Step>> = Mutex::new(None);
static NEXT_ID: AtomicUsize = AtomicUsize::new(0);

/// A step under way: count what it has done with [`Tracker::add`]; it ends
/// when dropped
pub struct Tracker {
    done: Arc<AtomicUsize>,
    detail: Arc<Mutex<String>>,
    id: usize,
}

impl Tracker {
    /// `n` more items done
    pub fn add(&self, n: usize) {
        self.done.fetch_add(n, Ordering::Relaxed);
    }

    /// What the step is busy with now, shown after its counts
    pub fn detail(&self, text: impl Into<String>) {
        if let Ok(mut d) = self.detail.lock() {
            *d = text.into();
        }
    }

    /// Counts this step's items from code that does not know the step (a
    /// report writer)
    pub fn counter(&self) -> Counter {
        Counter(Some(self.done.clone()))
    }
}

/// Counts a step's items from inside the code doing them; the default one
/// counts nothing (tests, callers without a step)
#[derive(Clone, Default)]
pub struct Counter(Option<Arc<AtomicUsize>>);

impl Counter {
    pub fn add(&self, n: usize) {
        if let Some(done) = &self.0 {
            done.fetch_add(n, Ordering::Relaxed);
        }
    }
}

impl Drop for Tracker {
    fn drop(&mut self) {
        if let Ok(mut step) = STEP.lock() {
            if step.as_ref().is_some_and(|s| s.id == self.id) {
                *step = None;
            }
        }
    }
}

/// Start a step of `total` `unit`s (0: no count, only the time it has taken)
pub fn step(label: impl Into<String>, total: usize, unit: &'static str) -> Tracker {
    let done = Arc::new(AtomicUsize::new(0));
    let detail = Arc::new(Mutex::new(String::new()));
    let id = NEXT_ID.fetch_add(1, Ordering::Relaxed);
    if let Ok(mut s) = STEP.lock() {
        *s = Some(Step {
            label: label.into(),
            detail: detail.clone(),
            unit,
            total,
            done: done.clone(),
            start: Instant::now(),
            id,
        });
    }
    Tracker { done, detail, id }
}

/// `1 h 02 min`, `4 min 05 s`, `12.3 s`, `85 ms`
pub fn duration(d: Duration) -> String {
    let s = d.as_secs_f64();
    if s < 1.0 {
        format!("{} ms", d.as_millis())
    } else if s < 60.0 {
        format!("{s:.1} s")
    } else if s < 3600.0 {
        format!("{} min {:02} s", d.as_secs() / 60, d.as_secs() % 60)
    } else {
        format!("{} h {:02} min", d.as_secs() / 3600, (d.as_secs() % 3600) / 60)
    }
}

/// `1,234,567`
pub fn thousands(n: usize) -> String {
    let digits = n.to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(c);
    }
    out
}

/// The `progress:` line of a step that has run for `elapsed`
fn report(label: &str, detail: &str, unit: &str, total: usize, done: usize, elapsed: Duration) -> String {
    let mut line = format!("progress: {label}");
    if total > 0 {
        let pct = 100.0 * done as f64 / total as f64;
        line += &format!(": {} of {} {unit} ({pct:.0}%)", thousands(done.min(total)), thousands(total));
        let secs = elapsed.as_secs_f64();
        if done > 0 && secs > 0.0 {
            let rate = done as f64 / secs;
            // a rate of files a second says little; spectra a second do
            if rate >= 1.0 {
                line += &format!(", {}/s", thousands(rate.round() as usize));
            }
            if done < total {
                let left = Duration::from_secs_f64((total - done) as f64 / rate);
                line += &format!(", about {} left", duration(left));
            }
        }
    }
    if !detail.is_empty() {
        line += &format!(", {detail}");
    }
    line + &format!(" ({} so far)", duration(elapsed))
}

/// Writes a `progress:` line about the step under way every [`INTERVAL`],
/// until dropped
pub struct Heartbeat {
    stop: Arc<AtomicBool>,
    handle: Option<JoinHandle<()>>,
}

impl Heartbeat {
    pub fn start() -> Self {
        let stop = Arc::new(AtomicBool::new(false));
        let flag = stop.clone();
        let handle = std::thread::Builder::new()
            .name("progress".into())
            .spawn(move || {
                // the step last reported, and when
                let mut last: Option<(usize, Instant)> = None;
                while !flag.load(Ordering::Relaxed) {
                    std::thread::sleep(Duration::from_millis(250));
                    let Ok(guard) = STEP.lock() else { continue };
                    let Some(step) = guard.as_ref() else { continue };
                    let since = match last {
                        Some((id, at)) if id == step.id => at,
                        _ => step.start,
                    };
                    if since.elapsed() < INTERVAL {
                        continue;
                    }
                    let detail = step.detail.lock().map(|d| d.clone()).unwrap_or_default();
                    let done = step.done.load(Ordering::Relaxed);
                    log::info!("{}", report(&step.label, &detail, step.unit, step.total, done, step.start.elapsed()));
                    last = Some((step.id, Instant::now()));
                }
            })
            .ok();
        Self { stop, handle }
    }
}

impl Drop for Heartbeat {
    fn drop(&mut self) {
        self.stop.store(true, Ordering::Relaxed);
        if let Some(handle) = self.handle.take() {
            let _ = handle.join();
        }
    }
}

/// Each spectrum file of a batch, logged as its last spectrum is searched,
/// with how many files the search has done in all:
/// `- cas9_b.raw: searched for crosslinks (2 of 3 files done)` (sage-gui
/// counts the files from these lines)
pub struct FilesSearched<'a> {
    what: &'static str,
    /// file id -> its place in `left` and `names`
    slots: std::collections::HashMap<usize, usize>,
    /// spectra of each file still to search
    left: Vec<AtomicUsize>,
    names: Vec<String>,
    /// files the search has done, over every batch
    done: &'a AtomicUsize,
    total: usize,
}

impl<'a> FilesSearched<'a> {
    /// The batch's files (id, name) and the file id of each spectrum to
    /// search; `done` counts the search's files over every batch, of `total`.
    /// A file with no spectrum to search is done at once.
    pub fn new(
        what: &'static str,
        files: Vec<(usize, String)>,
        spectra: impl Iterator<Item = usize>,
        done: &'a AtomicUsize,
        total: usize,
    ) -> Self {
        let slots: std::collections::HashMap<usize, usize> =
            files.iter().enumerate().map(|(slot, (id, _))| (*id, slot)).collect();
        let mut counts = vec![0usize; files.len()];
        for id in spectra {
            if let Some(&slot) = slots.get(&id) {
                counts[slot] += 1;
            }
        }
        let searched = FilesSearched {
            what,
            slots,
            left: counts.iter().map(|&n| AtomicUsize::new(n)).collect(),
            names: files.into_iter().map(|(_, name)| name).collect(),
            done,
            total,
        };
        for (slot, &n) in counts.iter().enumerate() {
            if n == 0 {
                searched.finished(slot);
            }
        }
        searched
    }

    /// A spectrum of file `file_id` has been searched
    pub fn spectrum_done(&self, file_id: usize) {
        if let Some(&slot) = self.slots.get(&file_id) {
            if self.left[slot].fetch_sub(1, Ordering::AcqRel) == 1 {
                self.finished(slot);
            }
        }
    }

    fn finished(&self, slot: usize) {
        let n = self.done.fetch_add(1, Ordering::AcqRel) + 1;
        log::info!("- {}: searched for {} ({} of {} files done)", self.names[slot], self.what, n, self.total);
    }
}

#[cfg(test)]
mod test {
    use super::*;

    /// A file is done with its last spectrum; one with none to search at once
    #[test]
    fn files_are_counted_as_their_last_spectrum_is_searched() {
        let done = AtomicUsize::new(1);
        let files = vec![(1, "b.raw".to_string()), (2, "c.raw".to_string()), (3, "d.raw".to_string())];
        let spectra = [1, 1, 2];
        let searched = FilesSearched::new("crosslinks", files, spectra.iter().copied(), &done, 4);
        // d.raw had nothing to search
        assert_eq!(done.load(Ordering::Relaxed), 2);
        searched.spectrum_done(1);
        assert_eq!(done.load(Ordering::Relaxed), 2);
        searched.spectrum_done(2);
        searched.spectrum_done(1);
        assert_eq!(done.load(Ordering::Relaxed), 4);
        // a spectrum of another batch's file changes nothing
        searched.spectrum_done(9);
        assert_eq!(done.load(Ordering::Relaxed), 4);
    }

    #[test]
    fn progress_lines_read_well() {
        assert_eq!(duration(Duration::from_millis(85)), "85 ms");
        assert_eq!(duration(Duration::from_secs_f64(12.34)), "12.3 s");
        assert_eq!(duration(Duration::from_secs(245)), "4 min 05 s");
        assert_eq!(duration(Duration::from_secs(3720)), "1 h 02 min");
        assert_eq!(thousands(1_234_567), "1,234,567");
        assert_eq!(thousands(999), "999");
        assert_eq!(
            report("crosslink search, files 1-8 of 21", "", "spectra", 1_000_000, 250_000, Duration::from_secs(30)),
            "progress: crosslink search, files 1-8 of 21: 250,000 of 1,000,000 spectra (25%), 8,333/s, about 1 min 30 s left (30.0 s so far)"
        );
        assert_eq!(
            report("reading files 9-16 of 21", "now file 11: F11.d", "files", 8, 2, Duration::from_secs(75)),
            "progress: reading files 9-16 of 21: 2 of 8 files (25%), about 3 min 45 s left, now file 11: F11.d (1 min 15 s so far)"
        );
        assert_eq!(report("FDR and reports", "", "", 0, 0, Duration::from_secs(20)), "progress: FDR and reports (20.0 s so far)");
    }

    /// Only the step under way is reported; a finished one leaves nothing
    #[test]
    fn a_finished_step_is_forgotten() {
        let a = step("one", 10, "files");
        a.add(3);
        assert!(STEP.lock().unwrap().as_ref().is_some_and(|s| s.label == "one"));
        let b = step("two", 0, "");
        drop(a); // an older step ending does not end the newer one
        assert!(STEP.lock().unwrap().as_ref().is_some_and(|s| s.label == "two"));
        drop(b);
        assert!(STEP.lock().unwrap().is_none());
    }
}
