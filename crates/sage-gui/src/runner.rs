//! Runs external programs (RAW conversion, then `sage.exe`) on a background
//! thread, collecting their console output for the log panel.

use std::ffi::OsString;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Duration;

const MAX_LOG_LINES: usize = 20_000;

/// One program invocation in a job
pub struct Step {
    pub title: String,
    pub program: PathBuf,
    pub args: Vec<OsString>,
    /// Skip the step when this file already exists (converted RAW files)
    pub skip_if_exists: Option<PathBuf>,
    /// Cap on the program's memory in bytes (see [`crate::system::MemoryCap`])
    pub memory_limit: Option<u64>,
}

#[derive(Clone, PartialEq, Eq)]
pub enum Status {
    Running(String),
    Finished,
    Failed(String),
    Cancelled,
}

pub struct Job {
    /// The latest output lines (older ones are dropped beyond `MAX_LOG_LINES`)
    pub log: Arc<Mutex<Vec<String>>>,
    /// Lines ever written to `log`, so readers can tell which ones are new
    pub total: Arc<AtomicUsize>,
    pub status: Arc<Mutex<Status>>,
    cancel: Arc<AtomicBool>,
}

impl Job {
    pub fn start(steps: Vec<Step>, ctx: eframe::egui::Context) -> Job {
        let job = Job {
            log: Arc::new(Mutex::new(Vec::new())),
            total: Arc::new(AtomicUsize::new(0)),
            status: Arc::new(Mutex::new(Status::Running("starting".into()))),
            cancel: Arc::new(AtomicBool::new(false)),
        };
        let log = job.log.clone();
        let total = job.total.clone();
        let status = job.status.clone();
        let cancel = job.cancel.clone();

        std::thread::spawn(move || {
            let push = |line: String| push_line(&log, &total, &ctx, line);
            let result = run_steps(&steps, &push, &status, &cancel);
            let end = match result {
                Ok(()) => Status::Finished,
                Err(_) if cancel.load(Ordering::Relaxed) => Status::Cancelled,
                Err(e) => Status::Failed(e),
            };
            match &end {
                Status::Finished => push("=== Search finished ===".into()),
                Status::Cancelled => push("=== Cancelled ===".into()),
                Status::Failed(e) => push(format!("=== Failed: {e} ===")),
                Status::Running(_) => {}
            }
            *status.lock().unwrap() = end;
            ctx.request_repaint();
        });
        job
    }

    pub fn cancel(&self) {
        self.cancel.store(true, Ordering::Relaxed);
    }

    pub fn status(&self) -> Status {
        self.status.lock().unwrap().clone()
    }

    pub fn is_running(&self) -> bool {
        matches!(self.status(), Status::Running(_))
    }
}

fn push_line(log: &Mutex<Vec<String>>, total: &AtomicUsize, ctx: &eframe::egui::Context, line: String) {
    let mut log = log.lock().unwrap();
    log.push(line);
    total.fetch_add(1, Ordering::Relaxed);
    if log.len() > MAX_LOG_LINES {
        let excess = log.len() - MAX_LOG_LINES;
        log.drain(..excess);
    }
    ctx.request_repaint();
}

fn run_steps(
    steps: &[Step],
    push: &(dyn Fn(String) + Sync),
    status: &Mutex<Status>,
    cancel: &AtomicBool,
) -> Result<(), String> {
    for step in steps {
        if cancel.load(Ordering::Relaxed) {
            return Err("cancelled".into());
        }
        if let Some(out) = &step.skip_if_exists {
            if out.exists() {
                push(format!("--- {}: {} already exists, skipped", step.title, out.display()));
                continue;
            }
        }
        *status.lock().unwrap() = Status::Running(step.title.clone());
        push(format!("--- {} ---", step.title));
        push(format!("> {}", command_line(&step.program, &step.args)));
        run_one(step, push, cancel)?;
    }
    Ok(())
}

fn command_line(program: &Path, args: &[OsString]) -> String {
    let quote = |s: String| if s.contains(' ') { format!("\"{s}\"") } else { s };
    std::iter::once(program.as_os_str())
        .chain(args.iter().map(OsString::as_os_str))
        .map(|s| quote(s.to_string_lossy().into_owned()))
        .collect::<Vec<_>>()
        .join(" ")
}

fn run_one(step: &Step, push: &(dyn Fn(String) + Sync), cancel: &AtomicBool) -> Result<(), String> {
    let mut cmd = Command::new(&step.program);
    cmd.args(&step.args)
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    #[cfg(windows)]
    {
        use std::os::windows::process::CommandExt;
        const CREATE_NO_WINDOW: u32 = 0x0800_0000;
        cmd.creation_flags(CREATE_NO_WINDOW);
    }
    let mut child: Child = cmd
        .spawn()
        .map_err(|e| format!("could not start {}: {e}", step.program.display()))?;
    let cap = match step.memory_limit {
        Some(limit) => match crate::system::MemoryCap::apply(&child, limit) {
            Ok(cap) => {
                push(format!("RAM limit: {}", crate::system::format_memory(limit)));
                Some(cap)
            }
            Err(e) => {
                let _ = child.kill();
                return Err(e);
            }
        },
        None => None,
    };
    let over_limit = |code| match &cap {
        Some(cap) if cap.reached() => format!(
            "{} ran out of memory: it needed more than the {} RAM limit (Files > Resources); \
             raise the limit, or search fewer files at once",
            step.title,
            crate::system::format_memory(cap.limit())
        ),
        _ => format!("{} exited with {code}", step.title),
    };

    let stdout = child.stdout.take().unwrap();
    let stderr = child.stderr.take().unwrap();

    std::thread::scope(|scope| {
        scope.spawn(|| read_lines(stdout, push));
        scope.spawn(|| read_lines(stderr, push));
        loop {
            if cancel.load(Ordering::Relaxed) {
                let _ = child.kill();
                let _ = child.wait();
                return Err("cancelled".to_string());
            }
            match child.try_wait() {
                Ok(Some(code)) if code.success() => return Ok(()),
                Ok(Some(code)) => return Err(over_limit(code)),
                Ok(None) => std::thread::sleep(Duration::from_millis(100)),
                Err(e) => return Err(e.to_string()),
            }
        }
    })
}

/// Forward output line by line. Progress lines that redraw themselves with
/// `\r` (msconvert) are reduced to their final state.
fn read_lines(source: impl Read, push: &(dyn Fn(String) + Sync)) {
    for line in BufReader::new(source).split(b'\n').map_while(Result::ok) {
        let text = String::from_utf8_lossy(&line);
        let last = text
            .trim_end_matches('\r')
            .rsplit('\r')
            .next()
            .unwrap_or_default()
            .to_string();
        if !last.trim().is_empty() {
            push(last);
        }
    }
}

/// `sage.exe` next to this program (both are built into `target\release`,
/// and released in one folder)
pub fn find_sage() -> Option<PathBuf> {
    let candidate = sage_beside()?;
    candidate.is_file().then_some(candidate)
}

/// Where `sage.exe` belongs: next to this program
pub fn sage_beside() -> Option<PathBuf> {
    let exe = std::env::current_exe().ok()?;
    let name = if cfg!(windows) { "sage.exe" } else { "sage" };
    Some(exe.parent()?.join(name))
}

/// The version a `sage.exe` reports for `--version` (`sage 1.8.6` → `1.8.6`);
/// `None` if it cannot be run or does not say
pub fn sage_version(sage: &Path) -> Option<String> {
    if !sage.is_file() {
        return None;
    }
    let mut cmd = Command::new(sage);
    cmd.arg("--version");
    #[cfg(windows)]
    {
        use std::os::windows::process::CommandExt;
        const CREATE_NO_WINDOW: u32 = 0x0800_0000;
        cmd.creation_flags(CREATE_NO_WINDOW);
    }
    let out = cmd.output().ok()?;
    parse_sage_version(&String::from_utf8_lossy(&out.stdout))
}

/// The version in `sage --version`'s answer (`sage 1.8.6`)
pub fn parse_sage_version(text: &str) -> Option<String> {
    text.lines().find_map(|line| {
        let mut words = line.split_whitespace();
        match (words.next(), words.next()) {
            (Some("sage"), Some(v)) if v.chars().next().is_some_and(|c| c.is_ascii_digit()) => Some(v.to_string()),
            _ => None,
        }
    })
}

/// What is wrong with the `sage.exe` sage-gui is set to run, as `(error?,
/// message)`: errors stop a search, warnings do not. `chosen` is the one
/// set (Files > Programs), `beside` the one next to sage-gui.exe if there is
/// one, `beside_place` where that would be, `version` what `chosen` reports
/// and `own` sage-gui's version. Empty when all is well.
pub fn sage_problems(chosen: &Path, beside: Option<&Path>, beside_place: Option<&Path>, version: Option<&str>, own: &str) -> Vec<(bool, String)> {
    let folder = beside_place.and_then(Path::parent).map_or_else(|| "sage-gui's folder".to_string(), |p| p.display().to_string());
    let release = format!("the sage.exe of this release (XNova {own})");
    if chosen.as_os_str().is_empty() || !chosen.is_file() {
        return vec![(
            true,
            format!("sage.exe not found: put {release} next to sage-gui.exe ({folder}), or set its location under Files > Programs."),
        )];
    }
    let mut out = Vec::new();
    match beside {
        None => out.push((
            false,
            format!(
                "There is no sage.exe next to sage-gui.exe ({folder}), so {} is used. sage.exe and \
                 sage-gui.exe of one release belong in the same folder.",
                chosen.display()
            ),
        )),
        Some(beside) if !same_file(beside, chosen) => out.push((
            false,
            format!("Using {} (set under Files > Programs), not the sage.exe next to sage-gui.exe ({}).", chosen.display(), beside.display()),
        )),
        Some(_) => {}
    }
    match version {
        Some(v) if v == own => {}
        Some(v) => out.push((
            true,
            format!(
                "{} is version {v}, but this sage-gui is version {own}: sage.exe may not read the \
                 settings sage-gui writes (an older one stops at settings it does not know). Use {release}: \
                 put it next to sage-gui.exe, or set it under Files > Programs.",
                chosen.display()
            ),
        )),
        None => out.push((
            false,
            format!(
                "{} did not report its version, so it may be from another release (this sage-gui is version \
                 {own}); if the search stops at settings it does not know, use {release}.",
                chosen.display()
            ),
        )),
    }
    out
}

/// Are two paths the same file (whatever their spelling)?
fn same_file(a: &Path, b: &Path) -> bool {
    match (a.canonicalize(), b.canonicalize()) {
        (Ok(x), Ok(y)) => x == y,
        _ => a == b,
    }
}

/// The newest ProteoWizard msconvert in the usual install locations
pub fn find_msconvert() -> Option<PathBuf> {
    let mut roots: Vec<PathBuf> = Vec::new();
    for var in ["ProgramFiles", "ProgramFiles(x86)"] {
        if let Ok(dir) = std::env::var(var) {
            roots.push(PathBuf::from(dir).join("ProteoWizard"));
        }
    }
    if let Ok(dir) = std::env::var("LOCALAPPDATA") {
        roots.push(PathBuf::from(dir).join("Apps"));
    }
    let mut found: Vec<PathBuf> = Vec::new();
    for root in roots {
        let direct = root.join("msconvert.exe");
        if direct.is_file() {
            found.push(direct);
        }
        let Ok(entries) = std::fs::read_dir(&root) else {
            continue;
        };
        for entry in entries.flatten() {
            let path = entry.path().join("msconvert.exe");
            let is_pwiz = entry.file_name().to_string_lossy().contains("ProteoWizard");
            if (is_pwiz || root.ends_with("ProteoWizard")) && path.is_file() {
                found.push(path);
            }
        }
    }
    // Install folders carry the version, so the last name is the newest
    found.sort();
    found.pop()
}

/// The command that converts one Thermo .raw file to mzML in `out_dir`.
/// Supports ProteoWizard msconvert and ThermoRawFileParser; `noise` asks
/// ThermoRawFileParser for its noise data (TMT signal-to-noise).
pub fn conversion_step(converter: &Path, raw: &Path, out_dir: &Path, noise: bool) -> Step {
    let output = converted_path(raw, out_dir);
    let args: Vec<OsString> = if is_thermorawfileparser(converter) {
        let mut args: Vec<OsString> = vec![
            "-i".into(),
            raw.into(),
            "-o".into(),
            out_dir.into(),
            "-f".into(),
            "2".into(),
        ];
        if noise {
            args.push("-N".into());
        }
        args
    } else {
        vec![
            raw.into(),
            "--mzML".into(),
            "--zlib".into(),
            "--filter".into(),
            "peakPicking vendor msLevel=1-".into(),
            "-o".into(),
            out_dir.into(),
        ]
    };
    Step {
        title: format!("Convert {}", raw.file_name().unwrap_or_default().to_string_lossy()),
        program: converter.to_path_buf(),
        args,
        skip_if_exists: Some(output),
        memory_limit: None,
    }
}

/// Is this converter ThermoRawFileParser (else msconvert)?
pub fn is_thermorawfileparser(converter: &Path) -> bool {
    converter
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .to_ascii_lowercase()
        .contains("thermorawfileparser")
}

pub fn converted_path(raw: &Path, out_dir: &Path) -> PathBuf {
    let stem = raw.file_stem().unwrap_or_default().to_string_lossy();
    out_dir.join(format!("{stem}.mzML"))
}

#[cfg(test)]
mod raw_conversion_test {
    use super::*;

    /// ThermoRawFileParser writes noise data (-N) only when asked, for TMT
    /// signal-to-noise; msconvert is run as before
    #[test]
    fn converter_arguments() {
        let raw = Path::new("C:/data/run1.raw");
        let out = Path::new("C:/out/converted_mzML");
        let trfp = Path::new("C:/Tools/ThermoRawFileParser/ThermoRawFileParser.exe");
        let args = |step: Step| step.args.iter().map(|a| a.to_string_lossy().into_owned()).collect::<Vec<_>>();
        assert!(is_thermorawfileparser(trfp));
        assert!(!is_thermorawfileparser(Path::new("C:/ProteoWizard/msconvert.exe")));
        let plain = args(conversion_step(trfp, raw, out, false));
        assert_eq!(plain.iter().filter(|a| *a == "-f").count(), 1);
        assert!(!plain.contains(&"-N".to_string()));
        assert!(args(conversion_step(trfp, raw, out, true)).ends_with(&["-N".to_string()]));
        let msconvert = args(conversion_step(Path::new("C:/ProteoWizard/msconvert.exe"), raw, out, true));
        assert!(msconvert.contains(&"--mzML".to_string()) && !msconvert.contains(&"-N".to_string()));
        assert_eq!(conversion_step(trfp, raw, out, false).skip_if_exists, Some(out.join("run1.mzML")));
    }
}

#[cfg(test)]
mod sage_exe_tests {
    use super::*;

    #[test]
    fn sage_versions_are_read() {
        assert_eq!(parse_sage_version("sage 1.8.6\n").as_deref(), Some("1.8.6"));
        assert_eq!(parse_sage_version("[log] something\nsage 1.0.0\r\n").as_deref(), Some("1.0.0"));
        assert_eq!(parse_sage_version("error: unexpected argument"), None);
        assert_eq!(parse_sage_version("sage"), None);
    }

    /// The sage.exe next to sage-gui is the one to use; a missing one is an
    /// error, one elsewhere a warning, one of another version an error
    #[test]
    fn sage_exe_problems() {
        let dir = std::env::temp_dir().join(format!("xnova_sage_check_{}", std::process::id()));
        let (gui_dir, other_dir) = (dir.join("release"), dir.join("old"));
        std::fs::create_dir_all(&gui_dir).unwrap();
        std::fs::create_dir_all(&other_dir).unwrap();
        let (beside, other) = (gui_dir.join("sage.exe"), other_dir.join("sage.exe"));
        std::fs::write(&beside, b"").unwrap();
        std::fs::write(&other, b"").unwrap();
        let errors = |p: &[(bool, String)]| p.iter().filter(|(e, _)| *e).count();
        let warnings = |p: &[(bool, String)]| p.iter().filter(|(e, _)| !*e).count();

        // next to sage-gui, same version: nothing to say
        assert!(sage_problems(&beside, Some(&beside), Some(&beside), Some("1.8.6"), "1.8.6").is_empty());
        // of another version: an error that names both versions
        let p = sage_problems(&beside, Some(&beside), Some(&beside), Some("1.0.0"), "1.8.6");
        assert_eq!((errors(&p), warnings(&p)), (1, 0));
        assert!(p[0].1.contains("version 1.0.0") && p[0].1.contains("version 1.8.6"), "{}", p[0].1);
        // none next to sage-gui, an old one elsewhere: a warning and an error
        let p = sage_problems(&other, None, Some(&beside.with_file_name("missing.exe")), Some("1.0.0"), "1.8.6");
        assert_eq!((errors(&p), warnings(&p)), (1, 1));
        assert!(p[0].1.contains("no sage.exe next to sage-gui.exe"), "{}", p[0].1);
        // another one chosen while one is next to sage-gui: a warning
        let p = sage_problems(&other, Some(&beside), Some(&beside), Some("1.8.6"), "1.8.6");
        assert_eq!((errors(&p), warnings(&p)), (0, 1));
        // no version reported: a warning
        let p = sage_problems(&beside, Some(&beside), Some(&beside), None, "1.8.6");
        assert_eq!((errors(&p), warnings(&p)), (0, 1));
        // none at all: an error, and nothing else
        let p = sage_problems(Path::new(""), None, Some(&beside), None, "1.8.6");
        assert_eq!((errors(&p), warnings(&p)), (1, 0));
        assert!(p[0].1.starts_with("sage.exe not found") && p[0].1.contains("XNova 1.8.6"), "{}", p[0].1);
        let _ = std::fs::remove_dir_all(&dir);
    }
}
