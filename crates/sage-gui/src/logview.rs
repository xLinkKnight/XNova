//! The log as shown on the Progress tab, copied and saved: each line with the
//! local time it arrived (in the user's format, see [`crate::clock`]), a line
//! with the date where a day begins, and sage.exe's messages without its own
//! UTC timestamp and module name.

use crate::clock::{self, Stamp};
use crate::runner::{Log, NOTE};

/// What stands before each line
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum TimeShown {
    /// The time of day (`1:32:22 PM`)
    #[default]
    Time,
    /// The date and the time (`11/19/2026 1:32:22 PM`)
    DateTime,
    /// The time since the search began (`0:12:05`)
    SinceStart,
    /// Nothing
    None,
}

impl TimeShown {
    pub const ALL: [TimeShown; 4] = [TimeShown::Time, TimeShown::DateTime, TimeShown::SinceStart, TimeShown::None];

    pub fn label(self) -> &'static str {
        match self {
            TimeShown::Time => "time of day",
            TimeShown::DateTime => "date and time",
            TimeShown::SinceStart => "time since start",
            TimeShown::None => "no time",
        }
    }

    /// Its name in sage-gui's saved state
    pub fn key(self) -> &'static str {
        match self {
            TimeShown::Time => "time",
            TimeShown::DateTime => "datetime",
            TimeShown::SinceStart => "since_start",
            TimeShown::None => "none",
        }
    }

    pub fn from_key(key: &str) -> TimeShown {
        TimeShown::ALL.into_iter().find(|t| t.key() == key).unwrap_or_default()
    }

    /// For a saved log: a clock time in any case (the time of day when the
    /// screen shows none or the time since the start)
    pub fn with_clock(self) -> TimeShown {
        match self {
            TimeShown::DateTime => TimeShown::DateTime,
            _ => TimeShown::Time,
        }
    }

    /// The text before a line that arrived at `stamp`
    pub fn text(self, stamp: &Stamp) -> String {
        match self {
            TimeShown::Time => stamp.time(),
            TimeShown::DateTime => format!("{} {}", stamp.short_date(), stamp.time()),
            TimeShown::SinceStart => clock::since(stamp.since_start),
            TimeShown::None => String::new(),
        }
    }
}

/// What a line is
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Info,
    /// sage.exe's detail (Verbose)
    Debug,
    Warn,
    Error,
    /// A stage beginning (`=== ... ===`, `--- Sage search ---`)
    Banner,
    /// sage.exe's progress on a long step
    Progress,
    /// sage-gui's own note
    Note,
    /// Other output (the command line, sage.exe's settings at its end)
    Plain,
}

/// The kind of `line` and the text to show: sage.exe's log lines
/// (`[2026-11-19T20:32:22Z INFO  sage_cli::runner] message`) lose their UTC
/// timestamp and module name
pub fn parse(line: &str) -> (Kind, &str) {
    if let Some(note) = line.strip_prefix(NOTE) {
        return (Kind::Note, note);
    }
    if line.starts_with("---") || line.starts_with("===") {
        let kind = if line.starts_with("=== Failed") { Kind::Error } else { Kind::Banner };
        return (kind, line);
    }
    if line.starts_with("Error") || line.contains("panicked at") {
        return (Kind::Error, line);
    }
    if let Some(rest) = line.strip_prefix('[') {
        if let Some((head, message)) = rest.split_once("] ") {
            let level = head.split_whitespace().nth(1).unwrap_or("");
            let kind = match level {
                "INFO" if message.starts_with("===") => Kind::Banner,
                "INFO" if message.starts_with("progress:") => Kind::Progress,
                "INFO" => Kind::Info,
                "WARN" => Kind::Warn,
                "ERROR" => Kind::Error,
                "DEBUG" | "TRACE" => Kind::Debug,
                _ => return (Kind::Plain, line),
            };
            return (kind, message);
        }
    }
    (Kind::Plain, line)
}

/// A row of the log as shown: a date line where a day begins, or a log line
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Row {
    /// The day of the log line at this index begins
    Day(usize),
    Line(usize),
}

/// The rows of `log`: each line, after a date line where its day begins
pub fn rows(log: &Log) -> Vec<Row> {
    let mut rows = Vec::with_capacity(log.lines.len() + 2);
    let mut day = None;
    for (i, stamp) in log.stamps.iter().enumerate().take(log.lines.len()) {
        if day != Some(stamp.date()) {
            rows.push(Row::Day(i));
            day = Some(stamp.date());
        }
        rows.push(Row::Line(i));
    }
    rows
}

/// The date line of a day (`──── Thursday, November 19, 2026 ────`)
pub fn day_text(stamp: &Stamp, ascii: bool) -> String {
    let rule = if ascii { "----" } else { "────" };
    format!("{rule} {} {rule}", stamp.long_date())
}

/// The marker before a line's text
pub fn marker(kind: Kind) -> &'static str {
    match kind {
        Kind::Warn => "⚠ ",
        Kind::Error => "✖ ",
        _ => "",
    }
}

/// The whole log as text, as on screen (with ASCII date lines): for "Copy
/// log", "Save log…" and the log saved next to the results
pub fn text(log: &Log, time: TimeShown) -> String {
    let width = time_width(log, time);
    let mut out = String::new();
    for row in rows(log) {
        match row {
            Row::Day(i) => {
                if time != TimeShown::SinceStart && time != TimeShown::None {
                    out += &day_text(&log.stamps[i], true);
                    out.push('\n');
                }
            }
            Row::Line(i) => {
                let (kind, message) = parse(&log.lines[i]);
                let stamp = time.text(&log.stamps[i]);
                if time != TimeShown::None {
                    out += &format!("{stamp:>width$}  ");
                }
                let mark = match kind {
                    Kind::Warn => "WARNING: ",
                    Kind::Error if !message.starts_with("Error") && !message.starts_with("===") => "ERROR: ",
                    _ => "",
                };
                out += mark;
                out += message;
                out.push('\n');
            }
        }
    }
    out
}

/// Width of the time column: the longest time text in the log
pub fn time_width(log: &Log, time: TimeShown) -> usize {
    if time == TimeShown::None {
        return 0;
    }
    // the first and the last line are enough: times grow in length only
    // from 9:59 to 10:00, and a log rarely spans more than a day
    let mut width = 0;
    for stamp in log.stamps.first().into_iter().chain(log.stamps.last()) {
        width = width.max(time.text(stamp).chars().count());
    }
    width.max(match time {
        TimeShown::Time => 8,
        TimeShown::DateTime => 18,
        TimeShown::SinceStart => 7,
        TimeShown::None => 0,
    })
}

/// sage.exe's warnings and errors in the log
pub fn problems(log: &Log) -> (usize, usize) {
    let mut warnings = 0;
    let mut errors = 0;
    for line in &log.lines {
        match parse(line).0 {
            Kind::Warn => warnings += 1,
            Kind::Error => errors += 1,
            _ => {}
        }
    }
    (warnings, errors)
}

#[cfg(test)]
mod test {
    use super::*;
    use std::time::Duration;

    fn stamp(day: u8, hour: u8, second: u8) -> Stamp {
        Stamp { year: 2026, month: 11, day, weekday: 4, hour, minute: 32, second, since_start: Duration::from_secs(second as u64) }
    }

    #[test]
    fn sage_lines_lose_their_utc_time_and_module() {
        assert_eq!(
            parse("[2026-11-19T20:32:22Z INFO  sage_cli::runner] generated 105116 fragments"),
            (Kind::Info, "generated 105116 fragments")
        );
        assert_eq!(parse("[2026-11-19T20:32:22Z WARN  sage_core::ml] Mobility failed"), (Kind::Warn, "Mobility failed"));
        assert_eq!(parse("[2026-11-19T20:32:22Z DEBUG sage_cli::runner]   file 1 of 2"), (Kind::Debug, "  file 1 of 2"));
        assert_eq!(parse("[2026-11-19T20:32:22Z INFO  sage] === crosslink results ==="), (Kind::Banner, "=== crosslink results ==="));
        assert_eq!(parse("[2026-11-19T20:32:22Z INFO  sage_cli::progress] progress: reading"), (Kind::Progress, "progress: reading"));
        assert_eq!(parse("# CPU 95%"), (Kind::Note, "CPU 95%"));
        assert_eq!(parse("--- Sage search ---"), (Kind::Banner, "--- Sage search ---"));
        assert_eq!(parse("=== Failed: sage.exe exited with 1 ==="), (Kind::Error, "=== Failed: sage.exe exited with 1 ==="));
        assert_eq!(parse("Error: no spectra"), (Kind::Error, "Error: no spectra"));
        assert_eq!(parse("{"), (Kind::Plain, "{"));
        assert_eq!(parse("[not a log line"), (Kind::Plain, "[not a log line"));
    }

    /// A date line before the first line and where the next day begins; the
    /// saved text has the times, the date lines and the warnings marked
    #[test]
    fn days_begin_with_a_date_line() {
        let log = Log {
            lines: vec![
                "--- Sage search ---".into(),
                "[2026-11-19T23:59:58Z WARN  sage] slow".into(),
                "[2026-11-20T00:00:01Z INFO  sage] next day".into(),
            ],
            stamps: vec![stamp(19, 23, 58), stamp(19, 23, 59), stamp(20, 0, 1)],
        };
        assert_eq!(rows(&log), [Row::Day(0), Row::Line(0), Row::Line(1), Row::Day(2), Row::Line(2)]);
        assert_eq!(problems(&log), (1, 0));
        let text = text(&log, TimeShown::SinceStart);
        assert_eq!(text, "0:00:58  --- Sage search ---\n0:00:59  WARNING: slow\n0:00:01  next day\n");
        let with_dates = super::text(&log, TimeShown::Time);
        assert_eq!(with_dates.lines().filter(|l| l.starts_with("---- ")).count(), 2, "{with_dates}");
    }
}
