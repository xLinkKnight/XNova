//! The local date and time of log lines, written the way the user's Windows
//! is set to write them (Settings > Time & language > Language & region >
//! Regional format): `11/19/2026  1:32:22 PM` in the US, `19/11/2026
//! 13:32:22` in the UK, `19.11.2026 13:32:22` in Germany, a 24-hour clock if
//! chosen there. Elsewhere than Windows: ISO dates and a 24-hour clock, in UTC.

use std::time::Duration;

/// When a log line arrived: the local calendar date and wall-clock time, and
/// the time since the job began
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Stamp {
    pub year: u16,
    pub month: u8,
    pub day: u8,
    /// 0 = Sunday
    pub weekday: u8,
    pub hour: u8,
    pub minute: u8,
    pub second: u8,
    pub since_start: Duration,
}

impl Stamp {
    /// Now, `since_start` after the job began
    pub fn now(since_start: Duration) -> Stamp {
        let (year, month, day, weekday, hour, minute, second) = local_now();
        Stamp { year, month, day, weekday, hour, minute, second, since_start }
    }

    /// The calendar day, to tell when a new one begins
    pub fn date(&self) -> (u16, u8, u8) {
        (self.year, self.month, self.day)
    }

    /// The time of day with seconds, in the user's format (`1:32:22 PM`)
    pub fn time(&self) -> String {
        format_time(self)
    }

    /// The date in the user's short format (`11/19/2026`)
    pub fn short_date(&self) -> String {
        format_date(self, false)
    }

    /// The date in the user's long format (`Thursday, November 19, 2026`)
    pub fn long_date(&self) -> String {
        format_date(self, true)
    }

    /// `2026-11-19_13-32-22`: sorts by time in any locale, for file names
    pub fn file_stamp(&self) -> String {
        format!(
            "{:04}-{:02}-{:02}_{:02}-{:02}-{:02}",
            self.year, self.month, self.day, self.hour, self.minute, self.second
        )
    }
}

/// `1:02:03` (hours, minutes, seconds) since the job began
pub fn since(d: Duration) -> String {
    let s = d.as_secs();
    format!("{}:{:02}:{:02}", s / 3600, (s / 60) % 60, s % 60)
}

#[cfg(windows)]
fn local_now() -> (u16, u8, u8, u8, u8, u8, u8) {
    let t = unsafe { windows::Win32::System::SystemInformation::GetLocalTime() };
    (t.wYear, t.wMonth as u8, t.wDay as u8, t.wDayOfWeek as u8, t.wHour as u8, t.wMinute as u8, t.wSecond as u8)
}

#[cfg(windows)]
fn system_time(s: &Stamp) -> windows::Win32::Foundation::SYSTEMTIME {
    windows::Win32::Foundation::SYSTEMTIME {
        wYear: s.year,
        wMonth: s.month as u16,
        wDayOfWeek: s.weekday as u16,
        wDay: s.day as u16,
        wHour: s.hour as u16,
        wMinute: s.minute as u16,
        wSecond: s.second as u16,
        wMilliseconds: 0,
    }
}

/// The user's locale (a null name) and format (a null format: the one set in
/// Windows), written into a buffer by `call`
#[cfg(windows)]
fn windows_text(call: impl Fn(Option<&mut [u16]>) -> i32) -> Option<String> {
    let mut buffer = [0u16; 128];
    let n = call(Some(&mut buffer));
    (n > 1).then(|| String::from_utf16_lossy(&buffer[..n as usize - 1]))
}

#[cfg(windows)]
fn format_time(s: &Stamp) -> String {
    use windows::core::PCWSTR;
    use windows::Win32::Globalization::{GetTimeFormatEx, TIME_FORMAT_FLAGS};
    let t = system_time(s);
    windows_text(|buffer| unsafe { GetTimeFormatEx(PCWSTR::null(), TIME_FORMAT_FLAGS(0), Some(&t), PCWSTR::null(), buffer) })
        .unwrap_or_else(|| format!("{:02}:{:02}:{:02}", s.hour, s.minute, s.second))
}

#[cfg(windows)]
fn format_date(s: &Stamp, long: bool) -> String {
    use windows::core::PCWSTR;
    use windows::Win32::Globalization::{GetDateFormatEx, DATE_LONGDATE, DATE_SHORTDATE};
    let t = system_time(s);
    let flags = if long { DATE_LONGDATE } else { DATE_SHORTDATE };
    windows_text(|buffer| unsafe { GetDateFormatEx(PCWSTR::null(), flags, Some(&t), PCWSTR::null(), buffer, PCWSTR::null()) })
        .unwrap_or_else(|| format!("{:04}-{:02}-{:02}", s.year, s.month, s.day))
}

/// UTC, from the system clock: days since 1970 to a calendar date
/// (H. Hinnant's `civil_from_days`)
#[cfg(not(windows))]
fn local_now() -> (u16, u8, u8, u8, u8, u8, u8) {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0) as i64;
    let (days, rem) = (secs.div_euclid(86_400), secs.rem_euclid(86_400));
    let weekday = ((days + 4).rem_euclid(7)) as u8; // 1970-01-01 was a Thursday
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = (doy - (153 * mp + 2) / 5 + 1) as u8;
    let month = if mp < 10 { mp + 3 } else { mp - 9 } as u8;
    let year = (yoe + era * 400 + if month <= 2 { 1 } else { 0 }) as u16;
    (year, month, day, weekday, (rem / 3600) as u8, ((rem / 60) % 60) as u8, (rem % 60) as u8)
}

#[cfg(not(windows))]
fn format_time(s: &Stamp) -> String {
    format!("{:02}:{:02}:{:02}", s.hour, s.minute, s.second)
}

#[cfg(not(windows))]
fn format_date(s: &Stamp, long: bool) -> String {
    let date = format!("{:04}-{:02}-{:02}", s.year, s.month, s.day);
    if long {
        let names = ["Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"];
        format!("{}, {date}", names[s.weekday as usize % 7])
    } else {
        date
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn stamps_read_as_dates_and_times() {
        let s = Stamp {
            year: 2026,
            month: 11,
            day: 19,
            weekday: 4,
            hour: 13,
            minute: 2,
            second: 3,
            since_start: Duration::from_secs(3725),
        };
        assert_eq!(s.file_stamp(), "2026-11-19_13-02-03");
        assert_eq!(since(s.since_start), "1:02:05");
        // whatever the locale, the year and the minutes are in it
        assert!(s.short_date().contains("2026") || s.short_date().contains("26"), "{}", s.short_date());
        assert!(s.long_date().contains("2026"), "{}", s.long_date());
        assert!(s.time().contains("02"), "{}", s.time());
        let now = Stamp::now(Duration::ZERO);
        assert!(now.year >= 2024 && (1..=12).contains(&now.month) && (1..=31).contains(&now.day));
    }
}
