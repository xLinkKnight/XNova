//! Live CPU, RAM, GPU and vRAM utilisation, sampled once a second on a
//! background thread while someone is looking: XNova's own (sage-gui and every
//! program it started, such as sage.exe and the RAW converters), and the whole
//! computer's for comparison.
//!
//! Sources match Task Manager. Whole computer: CPU from `GetSystemTimes`, RAM
//! from `GlobalMemoryStatusEx`, GPU from the "GPU Engine" and "GPU Adapter
//! Memory" performance counters, vRAM totals from DXGI. XNova: its processes
//! from a process snapshot (this process and its descendants), their CPU time
//! (`GetProcessTimes`) and working sets (`GetProcessMemoryInfo`), and their
//! instances of the "GPU Engine" and "GPU Process Memory" counters (Task
//! Manager's GPU and dedicated GPU memory columns).

use crate::system::Gpu;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Duration;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct GpuUsage {
    pub name: String,
    /// Busiest engine type (3D, compute, copy, video…), 0–100, as Task Manager reports
    pub util: f32,
    pub vram_used: u64,
    pub vram_total: u64,
    /// XNova's processes alone: their busiest engine type, 0–100
    pub own_util: f32,
    /// ...and the dedicated video memory they use
    pub own_vram: u64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Usage {
    /// Whole computer, 0–100
    pub cpu: f32,
    pub ram_used: u64,
    pub ram_total: u64,
    pub gpus: Vec<GpuUsage>,
    /// XNova's processes (sage-gui and what it started), how many
    pub own_processes: usize,
    /// ...their processor use, 0–100 of the whole computer (as Task Manager's CPU column)
    pub own_cpu: f32,
    /// ...and their working sets: the memory they hold in RAM
    pub own_ram: u64,
}

impl Usage {
    pub fn ram_pct(&self) -> f32 {
        if self.ram_total == 0 {
            0.0
        } else {
            100.0 * self.ram_used as f32 / self.ram_total as f32
        }
    }

    /// XNova's share of the RAM, 0–100
    pub fn own_ram_pct(&self) -> f32 {
        if self.ram_total == 0 {
            0.0
        } else {
            100.0 * self.own_ram as f32 / self.ram_total as f32
        }
    }

    /// The GPU worth showing: the one with the most dedicated memory (the discrete card)
    pub fn main_gpu(&self) -> Option<&GpuUsage> {
        self.gpus.iter().max_by_key(|g| g.vram_total)
    }
}

impl GpuUsage {
    pub fn vram_pct(&self) -> Option<f32> {
        (self.vram_total > 0).then(|| 100.0 * self.vram_used as f32 / self.vram_total as f32)
    }

    /// XNova's share of the video memory, 0–100
    pub fn own_vram_pct(&self) -> Option<f32> {
        (self.vram_total > 0).then(|| 100.0 * self.own_vram as f32 / self.vram_total as f32)
    }
}

/// Samples in the background; `latest` returns the newest sample
pub struct Monitor {
    latest: Arc<Mutex<Option<Usage>>>,
    /// Set by the GUI each frame it shows the numbers; the sampler idles otherwise
    wanted: Arc<AtomicBool>,
}

impl Monitor {
    pub fn start(gpus: Vec<Gpu>) -> Self {
        let monitor = Monitor { latest: Arc::default(), wanted: Arc::default() };
        let (latest, wanted) = (monitor.latest.clone(), monitor.wanted.clone());
        std::thread::spawn(move || {
            let mut sampler = Sampler::new(gpus);
            loop {
                if wanted.swap(false, Ordering::Relaxed) {
                    *latest.lock().unwrap() = sampler.sample();
                    std::thread::sleep(Duration::from_millis(1000));
                } else {
                    std::thread::sleep(Duration::from_millis(250));
                }
            }
        });
        monitor
    }

    /// The newest sample; also keeps sampling going and the screen updating
    pub fn latest(&self, ctx: &eframe::egui::Context) -> Option<Usage> {
        self.wanted.store(true, Ordering::Relaxed);
        ctx.request_repaint_after(Duration::from_millis(1000));
        self.latest.lock().unwrap().clone()
    }
}

/// The process id of a per-process counter instance (`pid_123_luid_…`)
fn instance_pid(name: &str) -> Option<u32> {
    let rest = name.strip_prefix("pid_")?;
    rest[..rest.find('_').unwrap_or(rest.len())].parse().ok()
}

/// Whether a counter instance belongs to one of `pids` (all instances count when `None`)
fn counts(name: &str, pids: Option<&[u32]>) -> bool {
    match pids {
        None => true,
        Some(pids) => instance_pid(name).is_some_and(|pid| pids.contains(&pid)),
    }
}

/// Task Manager's GPU figure from per-engine counter instances
/// (`pid_123_luid_0x…_0x…_phys_0_eng_0_engtype_3D` → value): per adapter, sum
/// each engine type over the processes (all, or `pids`), then take the busiest type.
fn gpu_utilisation(items: &[(String, f64)], luid: &str, pids: Option<&[u32]>) -> f32 {
    let mut by_type: std::collections::HashMap<&str, f64> = Default::default();
    for (name, value) in items {
        if !name.to_ascii_lowercase().contains(&luid.to_ascii_lowercase()) || !counts(name, pids) {
            continue;
        }
        let engtype = name.rsplit_once("engtype_").map_or("", |(_, t)| t);
        *by_type.entry(engtype).or_default() += value;
    }
    by_type.values().copied().fold(0.0, f64::max).min(100.0) as f32
}

/// Dedicated memory in use on one adapter, from "GPU Adapter Memory"
/// instances (whole computer) or "GPU Process Memory" instances of `pids`
fn vram_used(items: &[(String, f64)], luid: &str, pids: Option<&[u32]>) -> u64 {
    items
        .iter()
        .filter(|(name, _)| name.to_ascii_lowercase().contains(&luid.to_ascii_lowercase()) && counts(name, pids))
        .map(|(_, v)| *v)
        .sum::<f64>() as u64
}

/// `root` and every process descended from it, from `(pid, parent pid,
/// creation time)` of every process. A process counts as a child only when it
/// started after its parent: Windows reuses the ids of ended processes, so a
/// recorded parent id may now name an unrelated, newer process.
fn descendants(root: u32, processes: &[(u32, u32, u64)]) -> Vec<u32> {
    let created = |pid: u32| processes.iter().find(|p| p.0 == pid).map(|p| p.2);
    let mut out = vec![root];
    let mut i = 0;
    while i < out.len() {
        let parent = out[i];
        let since = created(parent).unwrap_or(0);
        for &(pid, ppid, at) in processes {
            if ppid == parent && pid != parent && at >= since && !out.contains(&pid) {
                out.push(pid);
            }
        }
        i += 1;
    }
    out
}

/// One XNova process at one sample
#[cfg(windows)]
#[derive(Clone, Copy)]
struct ProcessTimes {
    /// creation time and CPU time (kernel + user), in 100 ns units
    created: u64,
    busy: u64,
    working_set: u64,
}

#[cfg(windows)]
struct Sampler {
    gpus: Vec<Gpu>,
    /// Previous (idle, busy) CPU time of the whole computer
    cpu_prev: Option<(u64, u64)>,
    /// Previous sample of XNova's processes: when (system time, 100 ns) and each process's times
    own_prev: Option<(u64, std::collections::HashMap<u32, ProcessTimes>)>,
    query: isize,
    engine: isize,
    memory: isize,
    process_memory: isize,
}

#[cfg(windows)]
fn ticks(t: windows::Win32::Foundation::FILETIME) -> u64 {
    ((t.dwHighDateTime as u64) << 32) | t.dwLowDateTime as u64
}

#[cfg(windows)]
impl Sampler {
    fn new(gpus: Vec<Gpu>) -> Self {
        use windows::core::w;
        use windows::Win32::System::Performance::*;
        let (mut query, mut engine, mut memory, mut process_memory) = (0isize, 0isize, 0isize, 0isize);
        unsafe {
            if PdhOpenQueryW(None, 0, &mut query) == 0 {
                PdhAddEnglishCounterW(query, w!("\\GPU Engine(*)\\Utilization Percentage"), 0, &mut engine);
                PdhAddEnglishCounterW(query, w!("\\GPU Adapter Memory(*)\\Dedicated Usage"), 0, &mut memory);
                PdhAddEnglishCounterW(query, w!("\\GPU Process Memory(*)\\Dedicated Usage"), 0, &mut process_memory);
                // rates need two collections; the first one primes them
                PdhCollectQueryData(query);
            } else {
                query = 0;
            }
        }
        Sampler { gpus, cpu_prev: None, own_prev: None, query, engine, memory, process_memory }
    }

    fn cpu(&mut self) -> f32 {
        use windows::Win32::Foundation::FILETIME;
        use windows::Win32::System::Threading::GetSystemTimes;
        let (mut idle, mut kernel, mut user) = (FILETIME::default(), FILETIME::default(), FILETIME::default());
        if unsafe { GetSystemTimes(Some(&mut idle), Some(&mut kernel), Some(&mut user)) }.is_err() {
            return 0.0;
        }
        // kernel time includes idle time
        let (idle, busy) = (ticks(idle), ticks(kernel) + ticks(user) - ticks(idle));
        let pct = match self.cpu_prev {
            Some((pi, pb)) if (idle - pi) + (busy - pb) > 0 => {
                100.0 * (busy - pb) as f32 / ((idle - pi) + (busy - pb)) as f32
            }
            _ => 0.0,
        };
        self.cpu_prev = Some((idle, busy));
        pct
    }

    /// XNova's processes: this one and its descendants
    fn own_pids() -> Vec<u32> {
        use windows::Win32::Foundation::CloseHandle;
        use windows::Win32::System::Diagnostics::ToolHelp::*;
        let me = std::process::id();
        let mut processes = Vec::new();
        unsafe {
            let Ok(snapshot) = CreateToolhelp32Snapshot(TH32CS_SNAPPROCESS, 0) else {
                return vec![me];
            };
            let mut entry = PROCESSENTRY32W { dwSize: std::mem::size_of::<PROCESSENTRY32W>() as u32, ..Default::default() };
            let mut more = Process32FirstW(snapshot, &mut entry).is_ok();
            while more {
                processes.push((entry.th32ProcessID, entry.th32ParentProcessID));
                more = Process32NextW(snapshot, &mut entry).is_ok();
            }
            let _ = CloseHandle(snapshot);
        }
        // creation times only of the candidates: this process and those whose parent is in the tree
        let mut with_times: Vec<(u32, u32, u64)> = Vec::new();
        let mut tree = vec![me];
        let mut i = 0;
        while i < tree.len() {
            let parent = tree[i];
            if let Some(t) = Self::times(parent) {
                with_times.push((parent, processes.iter().find(|p| p.0 == parent).map_or(0, |p| p.1), t.created));
            }
            for &(pid, ppid) in &processes {
                if ppid == parent && pid != parent && !tree.contains(&pid) {
                    tree.push(pid);
                }
            }
            i += 1;
        }
        descendants(me, &with_times)
    }

    /// A process's creation time, CPU time and working set
    fn times(pid: u32) -> Option<ProcessTimes> {
        use windows::Win32::Foundation::{CloseHandle, FILETIME};
        use windows::Win32::System::ProcessStatus::{GetProcessMemoryInfo, PROCESS_MEMORY_COUNTERS};
        use windows::Win32::System::Threading::*;
        unsafe {
            let handle = OpenProcess(PROCESS_QUERY_LIMITED_INFORMATION | PROCESS_VM_READ, false, pid).ok()?;
            let (mut created, mut exited, mut kernel, mut user) =
                (FILETIME::default(), FILETIME::default(), FILETIME::default(), FILETIME::default());
            let got = GetProcessTimes(handle, &mut created, &mut exited, &mut kernel, &mut user);
            let mut memory = PROCESS_MEMORY_COUNTERS { cb: std::mem::size_of::<PROCESS_MEMORY_COUNTERS>() as u32, ..Default::default() };
            let working_set = match GetProcessMemoryInfo(handle, &mut memory, memory.cb) {
                Ok(()) => memory.WorkingSetSize as u64,
                Err(_) => 0,
            };
            let _ = CloseHandle(handle);
            got.ok()?;
            Some(ProcessTimes { created: ticks(created), busy: ticks(kernel) + ticks(user), working_set })
        }
    }

    /// XNova's processes, their CPU use since the last sample (0–100 of the
    /// whole computer) and their working sets
    fn own(&mut self) -> (Vec<u32>, f32, u64) {
        use windows::Win32::System::SystemInformation::GetSystemTimeAsFileTime;
        let pids = Self::own_pids();
        let now = ticks(unsafe { GetSystemTimeAsFileTime() });
        let current: std::collections::HashMap<u32, ProcessTimes> =
            pids.iter().filter_map(|&pid| Self::times(pid).map(|t| (pid, t))).collect();
        let ram = current.values().map(|t| t.working_set).sum();
        let cpu = match &self.own_prev {
            Some((then, previous)) if now > *then => {
                let mut busy = 0u64;
                for (pid, t) in &current {
                    busy += match previous.get(pid) {
                        // the same process as before
                        Some(p) if p.created == t.created => t.busy.saturating_sub(p.busy),
                        // started since the last sample: all its time is new
                        _ if t.created >= *then => t.busy,
                        // seen for the first time but older: no baseline yet
                        _ => 0,
                    };
                }
                let threads = std::thread::available_parallelism().map_or(1, |n| n.get()) as f64;
                (100.0 * busy as f64 / ((now - then) as f64 * threads)).clamp(0.0, 100.0) as f32
            }
            _ => 0.0,
        };
        self.own_prev = Some((now, current));
        (pids, cpu, ram)
    }

    /// Every instance of a wildcard counter as (instance name, value)
    fn counter_items(&self, counter: isize) -> Vec<(String, f64)> {
        use windows::Win32::System::Performance::*;
        if counter == 0 {
            return Vec::new();
        }
        unsafe {
            let (mut size, mut count) = (0u32, 0u32);
            let status = PdhGetFormattedCounterArrayW(counter, PDH_FMT_DOUBLE, &mut size, &mut count, None);
            if status != PDH_MORE_DATA || size == 0 {
                return Vec::new();
            }
            // the buffer holds the items followed by their names
            let mut buf = vec![0u64; (size as usize).div_ceil(8)];
            let items = buf.as_mut_ptr() as *mut PDH_FMT_COUNTERVALUE_ITEM_W;
            if PdhGetFormattedCounterArrayW(counter, PDH_FMT_DOUBLE, &mut size, &mut count, Some(items)) != 0 {
                return Vec::new();
            }
            std::slice::from_raw_parts(items, count as usize)
                .iter()
                .filter(|item| item.FmtValue.CStatus == 0)
                .map(|item| (item.szName.to_string().unwrap_or_default(), item.FmtValue.Anonymous.doubleValue))
                .collect()
        }
    }

    fn sample(&mut self) -> Option<Usage> {
        use windows::Win32::System::Performance::PdhCollectQueryData;
        use windows::Win32::System::SystemInformation::{GlobalMemoryStatusEx, MEMORYSTATUSEX};
        let cpu = self.cpu();
        let (pids, own_cpu, own_ram) = self.own();
        let mut mem = MEMORYSTATUSEX { dwLength: std::mem::size_of::<MEMORYSTATUSEX>() as u32, ..Default::default() };
        unsafe { GlobalMemoryStatusEx(&mut mem) }.ok()?;
        let (engines, memory, process_memory) = if self.query != 0 && unsafe { PdhCollectQueryData(self.query) } == 0 {
            (
                self.counter_items(self.engine),
                self.counter_items(self.memory),
                self.counter_items(self.process_memory),
            )
        } else {
            (Vec::new(), Vec::new(), Vec::new())
        };
        let gpus = self
            .gpus
            .iter()
            .map(|g| GpuUsage {
                name: g.name.clone(),
                util: gpu_utilisation(&engines, &g.luid, None),
                vram_used: vram_used(&memory, &g.luid, None),
                vram_total: g.vram,
                own_util: gpu_utilisation(&engines, &g.luid, Some(&pids)),
                own_vram: vram_used(&process_memory, &g.luid, Some(&pids)),
            })
            .collect();
        Some(Usage {
            cpu,
            ram_used: mem.ullTotalPhys - mem.ullAvailPhys,
            ram_total: mem.ullTotalPhys,
            gpus,
            own_processes: pids.len(),
            own_cpu,
            own_ram,
        })
    }
}

#[cfg(windows)]
impl Drop for Sampler {
    fn drop(&mut self) {
        if self.query != 0 {
            unsafe {
                windows::Win32::System::Performance::PdhCloseQuery(self.query);
            }
        }
    }
}

/// No usage figures outside Windows
#[cfg(not(windows))]
struct Sampler;

#[cfg(not(windows))]
impl Sampler {
    fn new(_gpus: Vec<Gpu>) -> Self {
        Sampler
    }

    fn sample(&mut self) -> Option<Usage> {
        None
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn gpu_figures_follow_task_manager() {
        let luid = "luid_0x00000000_0x0000D1C2";
        let items: Vec<(String, f64)> = [
            ("pid_10_luid_0x00000000_0x0000D1C2_phys_0_eng_0_engtype_3D", 30.0),
            ("pid_11_luid_0x00000000_0x0000D1C2_phys_0_eng_0_engtype_3D", 25.0),
            ("pid_11_luid_0x00000000_0x0000D1C2_phys_0_eng_5_engtype_Compute_0", 40.0),
            // another adapter
            ("pid_12_luid_0x00000000_0x0000AAAA_phys_0_eng_0_engtype_3D", 90.0),
        ]
        .iter()
        .map(|(n, v)| (n.to_string(), *v))
        .collect();
        // 3D: 30 + 25 = 55 beats Compute 40; the other adapter is ignored
        assert_eq!(gpu_utilisation(&items, luid, None), 55.0);
        // case of the hex digits does not matter
        assert_eq!(gpu_utilisation(&items, &luid.to_ascii_lowercase(), None), 55.0);
        // one process alone: its compute engine is its busiest
        assert_eq!(gpu_utilisation(&items, luid, Some(&[11])), 40.0);
        assert_eq!(gpu_utilisation(&items, luid, Some(&[10, 99])), 30.0);
        assert_eq!(gpu_utilisation(&items, luid, Some(&[12])), 0.0);

        let mem = vec![
            ("luid_0x00000000_0x0000D1C2_phys_0".to_string(), 1.5e9),
            ("luid_0x00000000_0x0000AAAA_phys_0".to_string(), 9.0e9),
        ];
        assert_eq!(vram_used(&mem, luid, None), 1_500_000_000);
        let per_process = vec![
            ("pid_10_luid_0x00000000_0x0000D1C2_phys_0".to_string(), 0.5e9),
            ("pid_11_luid_0x00000000_0x0000D1C2_phys_0".to_string(), 2.0e9),
            ("pid_11_luid_0x00000000_0x0000AAAA_phys_0".to_string(), 7.0e9),
        ];
        assert_eq!(vram_used(&per_process, luid, Some(&[11])), 2_000_000_000);
        assert_eq!(vram_used(&per_process, luid, Some(&[10, 11])), 2_500_000_000);
        assert_eq!(instance_pid("pid_4242_luid_0x0_0x1_phys_0"), Some(4242));
        assert_eq!(instance_pid("luid_0x0_0x1_phys_0"), None);
    }

    /// A process and its descendants; a process whose recorded parent id was
    /// reused by a newer process is not one of them
    #[test]
    fn descendants_follow_the_process_tree() {
        // (pid, parent, created)
        let processes = [(100, 1, 50), (200, 100, 60), (300, 200, 70), (400, 1, 55), (500, 400, 56), (600, 100, 40)];
        let mut tree = descendants(100, &processes);
        tree.sort();
        // 600 names 100 as its parent but started before it: an older process's child
        assert_eq!(tree, vec![100, 200, 300]);
        assert_eq!(descendants(400, &processes), vec![400, 500]);
    }

    /// A real sample on this machine: plausible numbers, matching the GPU
    /// list, with this test process counted as XNova's own
    #[cfg(windows)]
    #[test]
    fn samples_this_computer() {
        let gpus = crate::system::gpus();
        let mut s = Sampler::new(gpus.clone());
        s.sample();
        // keep a core busy so that the process's own CPU time shows
        let until = std::time::Instant::now() + Duration::from_millis(600);
        let mut x = 0u64;
        while std::time::Instant::now() < until {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1);
        }
        std::hint::black_box(x);
        let u = s.sample().unwrap();
        assert!((0.0..=100.0).contains(&u.cpu), "{u:?}");
        assert!(u.ram_total > 0 && u.ram_used <= u.ram_total);
        assert_eq!(u.gpus.len(), gpus.len());
        for g in &u.gpus {
            assert!((0.0..=100.0).contains(&g.util), "{g:?}");
            assert!(g.own_util <= 100.0 && g.own_vram <= g.vram_total.max(g.own_vram), "{g:?}");
        }
        assert!(u.own_processes >= 1, "{u:?}");
        assert!(u.own_ram > 0 && u.own_ram < u.ram_used, "{u:?}");
        assert!(u.own_cpu > 0.0 && u.own_cpu <= u.cpu + 5.0, "{u:?}");
        eprintln!("{u:#?}");
    }
}
