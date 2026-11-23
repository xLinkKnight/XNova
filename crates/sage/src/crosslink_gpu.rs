//! GPU implementation of the chain-scoring kernel (see `crosslink_kernel`).
//!
//! Built with `--features gpu`; uses `wgpu`, so the same code runs on DirectX 12
//! (Windows), Vulkan (Linux/Android) and Metal (macOS).
//!
//! One GPU thread scores one peptide against one query (a spectrum and one of
//! its precursor hypotheses), exactly as
//! [`crate::crosslink_kernel::score_peptide`] does on the CPU. The ion table is
//! uploaded once and stays resident. Queries come in batches of up to
//! [`MAX_BATCH`]: one round trip to the device scores every peptide against
//! every query of the batch, so the round trip's fixed cost (a few tenths of a
//! millisecond, more than the work itself) is paid once per batch instead of
//! once per query.
//!
//! The device also narrows each query's peptides before anything is copied
//! back. A histogram of the scores (64 bins per doubling) finds the bin that
//! holds the query's `max_survivors`-th best score; only peptides in that bin
//! or above are copied back. That is every survivor and a few more (the rest
//! of the bin), and the CPU then keeps exactly the best `max_survivors` with
//! [`top_survivors`], as when every score was copied back.
//!
//! Each fragment window is looked up in a spectrum's peaks starting from a
//! table: for every half-dalton bin, the first peak in that bin or above. The
//! window's exact edges then decide which peaks count, as a binary search
//! over all the peaks did (see `peak_bin` in the shader), so the scores are
//! the same; the table only saves the search's steps.
//!
//! Tables larger than the device's buffer limits are split into chunks of
//! consecutive peptides and dispatched one after another.

use crate::crosslink_kernel::{
    lnfact_table, top_survivors, ChainKernel, IonTable, KernelParams, Survivor,
};
use crate::spectrum::ProcessedSpectrum;
use rayon::prelude::*;
use std::borrow::Cow;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Mutex;
use std::time::Instant;
use wgpu::util::DeviceExt;

const WORKGROUP: u32 = 64;
const MAX_WORKGROUPS: u32 = 65_535;

/// Most queries in one batch
pub const MAX_BATCH: usize = 2048;
/// Device memory for the scores of a batch (every peptide against every query)
const SCORE_BUDGET: usize = 256 << 20;
/// Survivors the copy-back buffer holds at first (it grows when needed)
const FIRST_CAPACITY: usize = 1 << 20;
/// Score histogram bins per query, and the bit shift that makes a score's
/// bin: 64 bins per doubling from 1.0 up (lower scores share the first bin,
/// scores of 2^16 and above the last). Non-negative floats order like their
/// bits, so a bin never holds a better score than the bins above it.
const BINS: u32 = 1024;
const BIN_SHIFT: u32 = 17;
const BIN_BASE: u32 = 0x3F80_0000 >> BIN_SHIFT;

const SHADER: &str = r#"
// One chunk of peptides and one batch of queries
struct Batch {
    n_queries: u32,
    n_pep: u32,
    ion_base: u32,
    site_base: u32,
    // database index of the chunk's first peptide
    pep_base: u32,
    // survivors `entries` holds
    capacity: u32,
    // spectra in the batch (`spectra`)
    n_spectra: u32,
    pad1: u32,
}

// One spectrum and precursor hypothesis
struct Query {
    n_peaks: u32,
    peak_off: u32,
    max_charge: u32,
    min_matched: u32,
    precursor_mass: f32,
    // 1 + rel of the fragment tolerance (LinearTol::factors)
    lo_factor: f32,
    hi_factor: f32,
    abs_lo: f32,
    abs_hi: f32,
    max_survivors: u32,
    // the spectrum's peak table: where in `table`, and its bins
    table_off: u32,
    n_bins: u32,
    // per fragment charge z = 1, 2, ...: x = 1/z, y = slack/z (KernelParams::charge_terms)
    charges: array<vec2<f32>, 32>,
}

const BINS: u32 = __BINS__u;
const BIN_SHIFT: u32 = __BIN_SHIFT__u;
const BIN_BASE: u32 = __BIN_BASE__u;

@group(0) @binding(0) var<uniform> batch: Batch;
@group(0) @binding(1) var<storage, read> queries: array<Query>;
// (mass, intensity) of every query's spectrum, each sorted by mass
@group(0) @binding(2) var<storage, read> peaks: array<vec2<f32>>;
// (ion_off, site_off) of each peptide of the chunk, and one past the last
@group(0) @binding(3) var<storage, read> offs: array<vec2<u32>>;
// (mass bits, meta) of each ion
@group(0) @binding(4) var<storage, read> ions: array<vec2<u32>>;
@group(0) @binding(5) var<storage, read> sites: array<u32>;
@group(0) @binding(6) var<storage, read> pep_mass: array<f32>;
@group(0) @binding(7) var<storage, read> lnfact: array<f32>;
// (score bits, site << 16 | matched) of each query and peptide, query-major
@group(0) @binding(8) var<storage, read_write> scores: array<vec2<u32>>;
@group(0) @binding(9) var<storage, read_write> hist: array<atomic<u32>>;
@group(0) @binding(10) var<storage, read_write> thresh: array<u32>;
@group(0) @binding(11) var<storage, read_write> counter: array<atomic<u32>>;
// (query, peptide, score bits, info) of the peptides copied back
@group(0) @binding(12) var<storage, read_write> entries: array<vec4<u32>>;
// per spectrum: for each bin, the first of its peaks whose bin is that or above
@group(0) @binding(13) var<storage, read_write> table: array<u32>;
// per spectrum: (peak_off, n_peaks, table_off, n_bins)
@group(0) @binding(14) var<storage, read> spectra: array<vec4<u32>>;

// The histogram bin of a score (its bits)
fn bin_of(bits: u32) -> u32 {
    if (bits >= 0x80000000u) {
        return 0u;
    }
    return clamp(bits >> BIN_SHIFT, BIN_BASE, BIN_BASE + BINS - 1u) - BIN_BASE;
}

// The peak-table bin of a mass: floor(2 x mass), 0 below 0, at most
// n_bins - 1. Doubling is exact and floor is exact, so every device computes
// the same bin, and the bin never decreases as the mass grows: that is all
// the lookup below relies on.
fn peak_bin(mass: f32, n_bins: u32) -> u32 {
    let v = floor(mass * 2.0);
    if (!(v >= 0.0)) {
        return 0u;
    }
    if (v >= f32(n_bins - 1u)) {
        return n_bins - 1u;
    }
    return u32(v);
}

// Each spectrum's peak table: entry b is the first peak whose bin is b or
// more. x = bin, y = spectrum of the batch
@compute @workgroup_size(64)
fn build_table(@builtin(global_invocation_id) gid: vec3<u32>) {
    if (gid.y >= batch.n_spectra) {
        return;
    }
    let s = spectra[gid.y];
    let b = gid.x;
    if (b >= s.w) {
        return;
    }
    var l = 0u;
    var r = s.y;
    while (l < r) {
        let mid = (l + r) / 2u;
        if (peak_bin(peaks[s.x + mid].x, s.w) < b) {
            l = mid + 1u;
        } else {
            r = mid;
        }
    }
    table[s.z + b] = l;
}

// Intensity of the most intense peak in [lo, hi] (the later peak wins ties),
// or -1.0 when there is none. Peaks are sorted by mass.
//
// The first peak at lo or above, as a binary search finds it: every peak
// before the table's entry for lo's bin is in a lower bin, so lighter than lo
// (a bin never decreases with the mass); the few peaks of lo's own bin below
// lo are then stepped over.
fn best_peak(off: u32, n: u32, table_off: u32, n_bins: u32, lo: f32, hi: f32) -> f32 {
    var i = table[table_off + peak_bin(lo, n_bins)];
    while (i < n && peaks[off + i].x < lo) {
        i = i + 1u;
    }
    var best = -1.0;
    while (i < n && peaks[off + i].x <= hi) {
        if (peaks[off + i].y >= best) {
            best = peaks[off + i].y;
        }
        i = i + 1u;
    }
    return best;
}

// Every peptide of the chunk against every query: x = peptide, y = query
@compute @workgroup_size(64)
fn score_peptides(@builtin(global_invocation_id) gid: vec3<u32>) {
    let p = gid.x;
    let q = gid.y;
    if (p >= batch.n_pep) {
        return;
    }
    let cell = q * batch.n_pep + p;
    scores[cell] = vec2<u32>(bitcast<u32>(-1.0f), 0u);

    let o0 = offs[p];
    let o1 = offs[p + 1u];
    let i0 = o0.x - batch.ion_base;
    let i1 = o1.x - batch.ion_base;
    let s0 = o0.y - batch.site_base;
    let s1 = o1.y - batch.site_base;
    if (i0 == i1 || s0 == s1) {
        return;
    }
    let shift = queries[q].precursor_mass - pep_mass[p];
    if (shift <= 0.0) {
        return;
    }
    let n_sites = min(s1 - s0, 16u);
    let peak_off = queries[q].peak_off;
    let n_peaks = queries[q].n_peaks;
    let table_off = queries[q].table_off;
    let n_bins = queries[q].n_bins;
    let max_charge = queries[q].max_charge;
    let lo_factor = queries[q].lo_factor;
    let hi_factor = queries[q].hi_factor;
    let abs_lo = queries[q].abs_lo;
    let abs_hi = queries[q].abs_hi;

    var cnt_b: array<u32, 16>;
    var cnt_y: array<u32, 16>;
    var int_b: array<f32, 16>;
    var int_y: array<f32, 16>;
    for (var k = 0u; k < 16u; k = k + 1u) {
        cnt_b[k] = 0u;
        cnt_y[k] = 0u;
        int_b[k] = 0.0;
        int_y[k] = 0.0;
    }

    for (var ion = i0; ion < i1; ion = ion + 1u) {
        let ion_entry = ions[ion];
        let mass = bitcast<f32>(ion_entry.x);
        var un = 0u;
        var ui = 0.0;
        var sn = 0u;
        var si = 0.0;
        // Window edges exactly as LinearTol::bounds computes them, so that a
        // peak on an edge is in or out on every device (see its docs)
        for (var c = 1u; c < max_charge; c = c + 1u) {
            let k = queries[q].charges[c - 1u];
            let t = mass * k.x;
            let lo = t * lo_factor + abs_lo;
            let hi = t * hi_factor + abs_hi;
            let a = best_peak(peak_off, n_peaks, table_off, n_bins, lo, hi);
            if (a >= 0.0) {
                un = un + 1u;
                ui = ui + a;
            }
            let ts = (mass + shift) * k.x;
            let slo = ts * lo_factor + abs_lo - k.y;
            let shi = ts * hi_factor + abs_hi + k.y;
            let b = best_peak(peak_off, n_peaks, table_off, n_bins, slo, shi);
            if (b >= 0.0) {
                sn = sn + 1u;
                si = si + b;
            }
        }
        let ion_bits = ion_entry.y;
        let n_term = (ion_bits & 1u) == 1u;
        let idx = ion_bits >> 1u;
        for (var k = 0u; k < n_sites; k = k + 1u) {
            let site = sites[s0 + k];
            var contains = idx < site;
            if (n_term) {
                contains = idx >= site;
            }
            var n = un;
            var inten = ui;
            if (contains) {
                n = sn;
                inten = si;
            }
            if (n_term) {
                cnt_b[k] = cnt_b[k] + n;
                int_b[k] = int_b[k] + inten;
            } else {
                cnt_y[k] = cnt_y[k] + n;
                int_y[k] = int_y[k] + inten;
            }
        }
    }

    var best_score = -1.0;
    var best_info = 0u;
    for (var k = 0u; k < n_sites; k = k + 1u) {
        let matched = cnt_b[k] + cnt_y[k];
        let score = log((int_b[k] + 1.0) * (int_y[k] + 1.0))
            + lnfact[min(cnt_b[k], 255u)] + lnfact[min(cnt_y[k], 255u)];
        if (matched >= queries[q].min_matched && (best_score < 0.0 || score > best_score)) {
            best_score = score;
            best_info = (sites[s0 + k] << 16u) | min(matched, 65535u);
        }
    }
    if (best_score >= 0.0) {
        let bits = bitcast<u32>(best_score);
        scores[cell] = vec2<u32>(bits, best_info);
        atomicAdd(&hist[q * BINS + bin_of(bits)], 1u);
    }
}

// Per query: the bin holding its max_survivors-th best score (the first bin
// when it has fewer survivors)
@compute @workgroup_size(64)
fn threshold(@builtin(global_invocation_id) gid: vec3<u32>) {
    let q = gid.x;
    if (q >= batch.n_queries) {
        return;
    }
    let want = queries[q].max_survivors;
    var seen = 0u;
    var t = 0u;
    for (var b = BINS; b > 0u; b = b - 1u) {
        seen = seen + atomicLoad(&hist[q * BINS + b - 1u]);
        if (seen >= want) {
            t = b - 1u;
            break;
        }
    }
    thresh[q] = t;
}

// The survivors of every query at or above its bin, in any order
@compute @workgroup_size(64)
fn compact(@builtin(global_invocation_id) gid: vec3<u32>) {
    let p = gid.x;
    let q = gid.y;
    if (p >= batch.n_pep) {
        return;
    }
    let s = scores[q * batch.n_pep + p];
    if (!(bitcast<f32>(s.x) >= 0.0) || bin_of(s.x) < thresh[q]) {
        return;
    }
    let slot = atomicAdd(&counter[0], 1u);
    if (slot < batch.capacity) {
        entries[slot] = vec4<u32>(q, batch.pep_base + p, s.x, s.y);
    }
}
"#;

/// Fragment charges the shader's `charges` table holds
const MAX_CHARGE_TERMS: usize = 32;

/// Uniform block, must match `Batch` in the shader
#[repr(C)]
#[derive(Copy, Clone)]
struct GpuBatch {
    n_queries: u32,
    n_pep: u32,
    ion_base: u32,
    site_base: u32,
    pep_base: u32,
    capacity: u32,
    n_spectra: u32,
    pad1: u32,
}

/// Must match `Query` in the shader (12 four-byte words, then the charges)
#[repr(C)]
#[derive(Copy, Clone)]
struct GpuQuery {
    n_peaks: u32,
    peak_off: u32,
    max_charge: u32,
    min_matched: u32,
    precursor_mass: f32,
    lo_factor: f32,
    hi_factor: f32,
    abs_lo: f32,
    abs_hi: f32,
    max_survivors: u32,
    table_off: u32,
    n_bins: u32,
    charges: [[f32; 2]; MAX_CHARGE_TERMS],
}

impl GpuQuery {
    fn new(spectrum: &GpuSpectrum, params: &KernelParams) -> Self {
        // Fragment charges beyond the table (precursors above 33+) are not tried
        let mut charges = [[0f32; 2]; MAX_CHARGE_TERMS];
        let mut n_charges = 0;
        for (slot, (inv_z, slack_z)) in charges.iter_mut().zip(params.charge_terms()) {
            *slot = [inv_z, slack_z];
            n_charges += 1;
        }
        let (lo_factor, hi_factor) = params.fragment_tol.factors();
        Self {
            n_peaks: spectrum.n_peaks,
            peak_off: spectrum.peak_off,
            max_charge: n_charges as u32 + 1,
            min_matched: params.min_matched as u32,
            precursor_mass: params.precursor_mass,
            lo_factor,
            hi_factor,
            abs_lo: params.fragment_tol.abs_lo,
            abs_hi: params.fragment_tol.abs_hi,
            max_survivors: params.max_survivors.min(u32::MAX as usize) as u32,
            table_off: spectrum.table_off,
            n_bins: spectrum.n_bins,
            charges,
        }
    }
}

/// A spectrum of a batch, as `spectra` in the shader holds it
#[repr(C)]
#[derive(Copy, Clone, Debug, PartialEq)]
struct GpuSpectrum {
    peak_off: u32,
    n_peaks: u32,
    table_off: u32,
    n_bins: u32,
}

impl GpuSpectrum {
    /// Its peaks from `peak_off` and its table from `table_off`: a bin for
    /// each half-dalton up to its heaviest peak, and one above
    fn new(masses: &[f32], peak_off: usize, table_off: usize) -> Self {
        let heaviest = masses.iter().copied().fold(0.0f32, f32::max);
        // (at most 2^20 bins: masses beyond 524,288 Da share the last; the lookup stays exact)
        let n_bins = if heaviest.is_finite() { ((heaviest * 2.0).floor() as u32).min(1 << 20) + 2 } else { 2 };
        GpuSpectrum { peak_off: peak_off as u32, n_peaks: masses.len() as u32, table_off: table_off as u32, n_bins }
    }
}

fn bytes_of<T: Copy>(v: &[T]) -> &[u8] {
    // Safety: plain-old-data slices of u32/f32/#[repr(C)] structs of the same
    unsafe { std::slice::from_raw_parts(v.as_ptr() as *const u8, std::mem::size_of_val(v)) }
}

/// A run of consecutive peptides small enough for the device's buffer limits
struct Chunk {
    /// first peptide (database index) in this chunk, and how many
    first: usize,
    count: usize,
    ion_base: u32,
    site_base: u32,
    offs: wgpu::Buffer,
    ions: wgpu::Buffer,
    sites: wgpu::Buffer,
    pep_mass: wgpu::Buffer,
}

/// A batch's spectra and their peak tables (see `build_table` in the shader)
struct PeakTables {
    table: wgpu::Buffer,
    spectra: wgpu::Buffer,
    /// spectra in the batch, and the most bins of any
    count: u32,
    most_bins: u32,
}

/// Device buffers of a batch, kept from one batch to the next
struct Work {
    /// per chunk: every query's score of every peptide
    scores: Vec<wgpu::Buffer>,
    hist: wgpu::Buffer,
    thresh: wgpu::Buffer,
    counter: wgpu::Buffer,
    read_counter: wgpu::Buffer,
    entries: wgpu::Buffer,
    read_entries: wgpu::Buffer,
    /// survivors `entries` and `read_entries` hold
    capacity: usize,
}

/// Adapters on the backends the kernel uses, as `name | backend | device type`.
/// Software adapters are listed too; the kernel skips them unless asked.
pub fn adapters() -> Vec<String> {
    all_adapters()
        .into_iter()
        .map(|a| {
            let info = a.get_info();
            format!("{} | {:?} | {:?}", info.name, info.backend, info.device_type)
        })
        .collect()
}

/// Graphics APIs used: DirectX 12 on Windows, where it is the native API (some
/// Vulkan drivers there, e.g. Intel UHD 31.0.101.2145, crash opening a device);
/// Vulkan and Metal elsewhere. Never GL: slow and noisy.
fn backends() -> wgpu::Backends {
    if cfg!(windows) {
        wgpu::Backends::DX12
    } else {
        wgpu::Backends::PRIMARY
    }
}

/// One instance for the whole program: adapters and devices must not outlive it
fn instance() -> &'static wgpu::Instance {
    static INSTANCE: std::sync::OnceLock<wgpu::Instance> = std::sync::OnceLock::new();
    INSTANCE.get_or_init(|| {
        let mut descriptor = wgpu::InstanceDescriptor::new_without_display_handle();
        descriptor.backends = backends();
        wgpu::Instance::new(descriptor)
    })
}

fn all_adapters() -> Vec<wgpu::Adapter> {
    pollster::block_on(instance().enumerate_adapters(backends()))
}

/// Which of `names` (the adapters as `--gpu-info` lists them, in its order)
/// `want` picks: `#n`, the n-th; else a name equal to it, then the one name
/// containing it (case is ignored)
pub fn pick_adapter(names: &[String], want: &str) -> Result<usize, String> {
    let want = want.trim();
    let listed = || names.iter().enumerate().map(|(i, n)| format!("#{} {n}", i + 1)).collect::<Vec<_>>().join("; ");
    if let Some(n) = want.strip_prefix('#').and_then(|n| n.trim().parse::<usize>().ok()) {
        return match n.checked_sub(1).filter(|i| *i < names.len()) {
            Some(i) => Ok(i),
            None => Err(format!("there is no GPU #{n} on this computer ({})", listed())),
        };
    }
    let lower = want.to_lowercase();
    if let Some(i) = names.iter().position(|n| n.to_lowercase() == lower) {
        return Ok(i);
    }
    let matching: Vec<usize> = (0..names.len()).filter(|&i| names[i].to_lowercase().contains(&lower)).collect();
    match matching.as_slice() {
        [i] => Ok(*i),
        [] => Err(format!("no GPU named `{want}` on this computer ({})", listed())),
        _ => Err(format!("`{want}` names several GPUs ({}); give #n or the whole name", listed())),
    }
}

pub struct GpuKernel {
    device: wgpu::Device,
    queue: wgpu::Queue,
    /// the shader's steps: the spectra's peak tables, score, find each
    /// query's threshold, copy back
    build_table: wgpu::ComputePipeline,
    score: wgpu::ComputePipeline,
    threshold: wgpu::ComputePipeline,
    compact: wgpu::ComputePipeline,
    chunks: Vec<Chunk>,
    lnfact: wgpu::Buffer,
    adapter_name: String,
    label: String,
    /// most queries in one batch (the scores of a batch take `SCORE_BUDGET`)
    max_batch: usize,
    /// most survivors one batch may copy back; a batch with more is halved
    max_entries: usize,
    /// one batch at a time uses the device and these buffers
    work: Mutex<Work>,
    /// where the time of the chain filter goes, logged when the kernel ends
    stats: Stats,
}

/// Counts, and nanoseconds spent in each part of [`GpuKernel::survivors_many`]
#[derive(Default)]
struct Stats {
    queries: AtomicU64,
    batches: AtomicU64,
    /// waiting for another caller's turn on the device
    queued: AtomicU64,
    /// gathering the peaks and parameters, uploading them, binding the buffers
    setup: AtomicU64,
    /// from submitting the work to its results being readable: the device's
    /// work and the copies back
    device: AtomicU64,
    /// sorting the survivors copied back to their queries and keeping the best
    ranking: AtomicU64,
}

impl Stats {
    fn add(counter: &AtomicU64, since: Instant) -> Instant {
        let now = Instant::now();
        counter.fetch_add((now - since).as_nanos() as u64, Ordering::Relaxed);
        now
    }
}

impl Drop for GpuKernel {
    fn drop(&mut self) {
        let n = self.stats.queries.load(Ordering::Relaxed);
        if n == 0 {
            return;
        }
        let ms = |c: &AtomicU64| c.load(Ordering::Relaxed) as f64 / 1e6;
        let peptides: usize = self.chunks.iter().map(|c| c.count).sum();
        let (device, ranking, setup) = (ms(&self.stats.device), ms(&self.stats.ranking), ms(&self.stats.setup));
        log::info!(
            "gpu chain filter on {}: {n} queries (spectra x precursor hypotheses) x {peptides} peptides in {} \
             batches of up to {}; {:.1} s on the device (work and copies back), {:.1} s ranking the survivors \
             on the CPU, {:.1} s setting up: {:.4} ms per query",
            self.adapter_name,
            self.stats.batches.load(Ordering::Relaxed),
            self.max_batch,
            device / 1e3,
            ranking / 1e3,
            setup / 1e3,
            (device + ranking + setup) / n as f64,
        );
    }
}

impl GpuKernel {
    /// Upload `table` to an adapter: `device` (see [`pick_adapter`]), or the
    /// one the driver ranks fastest. Fails (rather than silently using a CPU
    /// adapter) when there is none.
    ///
    /// Software adapters (llvmpipe, WARP, ...) are CPUs in disguise and much
    /// slower than the `cpu` backend, so they are used only when `allow_software`
    /// is set (explicit `"backend": "gpu"`, or testing) or chosen by name.
    pub fn new(table: IonTable, allow_software: bool, device: Option<&str>) -> Result<Self, String> {
        let adapter = match device {
            Some(want) => {
                let mut adapters = all_adapters();
                let names: Vec<String> = adapters.iter().map(|a| a.get_info().name).collect();
                let i = pick_adapter(&names, want)?;
                let adapter = adapters.swap_remove(i);
                log::info!("GPU: {} (chosen with --gpu {want})", names[i]);
                adapter
            }
            None => {
                let adapter = pollster::block_on(instance().request_adapter(&wgpu::RequestAdapterOptions {
                    power_preference: wgpu::PowerPreference::HighPerformance,
                    ..Default::default()
                }))
                .map_err(|e| format!("no GPU adapter found: {e}"))?;
                let info = adapter.get_info();
                if info.device_type == wgpu::DeviceType::Cpu && !allow_software {
                    return Err(format!("only a software adapter ({}) is available", info.name));
                }
                // several GPUs: say which one, and how to choose another
                let others: Vec<String> = all_adapters()
                    .iter()
                    .map(|a| a.get_info())
                    .filter(|i| i.device_type != wgpu::DeviceType::Cpu && i.name != info.name)
                    .map(|i| i.name)
                    .collect();
                if !others.is_empty() {
                    log::info!(
                        "GPU: {} (the driver's fastest); also on this computer: {} (choose one with --gpu, or in \
                         sage-gui: Files > Resources > GPU to search on)",
                        info.name,
                        others.join(", ")
                    );
                }
                adapter
            }
        };
        Self::on_adapter(table, adapter)
    }

    /// One kernel per adapter on this machine, software ones included, so
    /// that tests can check every device and driver
    #[cfg(test)]
    pub(crate) fn on_every_adapter(table: impl Fn() -> IonTable) -> Vec<Result<Self, String>> {
        all_adapters().into_iter().map(|a| Self::on_adapter(table(), a)).collect()
    }

    /// Small batches and buffers, so that tests reach the paths large
    /// searches take: several batches, a copy-back buffer that grows, a batch
    /// halved, one query's scores read whole
    #[cfg(test)]
    pub(crate) fn squeeze(&mut self, max_batch: usize, capacity: usize, max_entries: usize) {
        self.max_batch = max_batch.clamp(1, self.max_batch);
        self.max_entries = max_entries.max(1);
        let work = self.work.get_mut().unwrap_or_else(|e| e.into_inner());
        Self::resize(&self.device, work, capacity.max(1));
    }

    fn on_adapter(table: IonTable, adapter: wgpu::Adapter) -> Result<Self, String> {
        let info = adapter.get_info();
        let limits = adapter.limits();
        let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
            label: Some("sage-crosslink"),
            required_limits: wgpu::Limits {
                max_storage_buffer_binding_size: limits.max_storage_buffer_binding_size,
                max_buffer_size: limits.max_buffer_size,
                max_storage_buffers_per_shader_stage: 12.max(
                    wgpu::Limits::default().max_storage_buffers_per_shader_stage,
                ).min(limits.max_storage_buffers_per_shader_stage),
                ..wgpu::Limits::default()
            },
            ..Default::default()
        }))
        .map_err(|e| format!("could not open GPU device: {e}"))?;

        let source = SHADER
            .replace("__BINS__", &BINS.to_string())
            .replace("__BIN_SHIFT__", &BIN_SHIFT.to_string())
            .replace("__BIN_BASE__", &BIN_BASE.to_string());
        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("crosslink-chain"),
            source: wgpu::ShaderSource::Wgsl(Cow::Owned(source)),
        });
        let pipeline = |entry: &str| {
            device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some(entry),
                layout: None,
                module: &module,
                entry_point: Some(entry),
                compilation_options: Default::default(),
                cache: None,
            })
        };
        let (build_table, score, threshold, compact) =
            (pipeline("build_table"), pipeline("score_peptides"), pipeline("threshold"), pipeline("compact"));

        let max_binding = limits.max_storage_buffer_binding_size.min(limits.max_buffer_size) as usize;
        let chunks = build_chunks(&device, &table, max_binding)?;
        let lnfact = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("lnfact"),
            contents: bytes_of(&lnfact_table()),
            usage: wgpu::BufferUsages::STORAGE,
        });

        // a batch's scores: 8 bytes per query and peptide, in one buffer per chunk
        let total: usize = chunks.iter().map(|c| c.count).sum::<usize>().max(1);
        let largest = chunks.iter().map(|c| c.count).max().unwrap_or(1).max(1);
        let max_batch = (SCORE_BUDGET / (total * 8)).min(max_binding / (largest * 8)).clamp(1, MAX_BATCH);
        let max_entries = (max_binding / 16).min(1 << 24);
        let work = Work::new(&device, &chunks, max_batch, FIRST_CAPACITY.min(max_entries));

        let label = format!(
            "gpu ({}, {:?}, {} chunk(s), batches of up to {max_batch} queries)",
            info.name,
            info.backend,
            chunks.len()
        );
        Ok(Self {
            device,
            queue,
            build_table,
            score,
            threshold,
            compact,
            chunks,
            lnfact,
            adapter_name: info.name,
            label,
            max_batch,
            max_entries,
            work: Mutex::new(work),
            stats: Stats::default(),
        })
    }

    pub fn adapter_name(&self) -> &str {
        &self.adapter_name
    }

    /// New copy-back buffers holding `capacity` survivors
    fn resize(device: &wgpu::Device, work: &mut Work, capacity: usize) {
        let size = (capacity * 16) as u64;
        work.entries = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("entries"),
            size,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        work.read_entries = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("read_entries"),
            size,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        work.capacity = capacity;
    }

    fn bind(&self, pipeline: &wgpu::ComputePipeline, buffers: &[(u32, &wgpu::Buffer)]) -> wgpu::BindGroup {
        let entries: Vec<wgpu::BindGroupEntry> = buffers
            .iter()
            .map(|&(binding, buffer)| wgpu::BindGroupEntry {
                binding,
                resource: buffer.as_entire_binding(),
            })
            .collect();
        self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("crosslink-chain"),
            layout: &pipeline.get_bind_group_layout(0),
            entries: &entries,
        })
    }

    /// Each chunk's `Batch` uniform
    fn uniforms(&self, n_queries: usize, n_spectra: usize, capacity: usize) -> Vec<wgpu::Buffer> {
        self.chunks
            .iter()
            .map(|c| {
                let batch = GpuBatch {
                    n_queries: n_queries as u32,
                    n_pep: c.count as u32,
                    ion_base: c.ion_base,
                    site_base: c.site_base,
                    pep_base: c.first as u32,
                    capacity: capacity.min(u32::MAX as usize) as u32,
                    n_spectra: n_spectra as u32,
                    pad1: 0,
                };
                self.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some("batch"),
                    contents: bytes_of(std::slice::from_ref(&batch)),
                    usage: wgpu::BufferUsages::UNIFORM,
                })
            })
            .collect()
    }

    /// Score every peptide against every query, and find each query's threshold
    fn encode_scoring(
        &self,
        encoder: &mut wgpu::CommandEncoder,
        work: &Work,
        uniforms: &[wgpu::Buffer],
        queries: &wgpu::Buffer,
        peaks: &wgpu::Buffer,
        tables: &PeakTables,
        n_queries: usize,
    ) {
        encoder.clear_buffer(&work.hist, 0, Some((n_queries * BINS as usize * 4) as u64));
        let build = self.bind(&self.build_table, &[(0, &uniforms[0]), (2, peaks), (13, &tables.table), (14, &tables.spectra)]);
        let groups: Vec<wgpu::BindGroup> = self
            .chunks
            .iter()
            .zip(uniforms)
            .zip(&work.scores)
            .map(|((chunk, uniform), scores)| {
                self.bind(
                    &self.score,
                    &[
                        (0, uniform),
                        (1, queries),
                        (2, peaks),
                        (3, &chunk.offs),
                        (4, &chunk.ions),
                        (5, &chunk.sites),
                        (6, &chunk.pep_mass),
                        (7, &self.lnfact),
                        (8, scores),
                        (9, &work.hist),
                        (13, &tables.table),
                    ],
                )
            })
            .collect();
        let thresholds = self.bind(&self.threshold, &[(0, &uniforms[0]), (1, queries), (9, &work.hist), (10, &work.thresh)]);
        let mut pass = encoder.begin_compute_pass(&Default::default());
        pass.set_pipeline(&self.build_table);
        pass.set_bind_group(0, &build, &[]);
        pass.dispatch_workgroups(tables.most_bins.div_ceil(WORKGROUP), tables.count, 1);
        pass.set_pipeline(&self.score);
        for (chunk, group) in self.chunks.iter().zip(&groups) {
            pass.set_bind_group(0, group, &[]);
            pass.dispatch_workgroups((chunk.count as u32).div_ceil(WORKGROUP), n_queries as u32, 1);
        }
        pass.set_pipeline(&self.threshold);
        pass.set_bind_group(0, &thresholds, &[]);
        pass.dispatch_workgroups((n_queries as u32).div_ceil(WORKGROUP), 1, 1);
    }

    /// Gather the survivors at or above each query's threshold, and copy
    /// their number back
    fn encode_compact(&self, encoder: &mut wgpu::CommandEncoder, work: &Work, uniforms: &[wgpu::Buffer], n_queries: usize) {
        encoder.clear_buffer(&work.counter, 0, None);
        let groups: Vec<wgpu::BindGroup> = uniforms
            .iter()
            .zip(&work.scores)
            .map(|(uniform, scores)| {
                self.bind(
                    &self.compact,
                    &[(0, uniform), (8, scores), (10, &work.thresh), (11, &work.counter), (12, &work.entries)],
                )
            })
            .collect();
        {
            let mut pass = encoder.begin_compute_pass(&Default::default());
            pass.set_pipeline(&self.compact);
            for (chunk, group) in self.chunks.iter().zip(&groups) {
                pass.set_bind_group(0, group, &[]);
                pass.dispatch_workgroups((chunk.count as u32).div_ceil(WORKGROUP), n_queries as u32, 1);
            }
        }
        encoder.copy_buffer_to_buffer(&work.counter, 0, &work.read_counter, 0, 4);
    }

    /// Wait for the device, then map the first `bytes` of `buffer`
    fn map_read(&self, buffer: &wgpu::Buffer, bytes: u64) {
        let (tx, rx) = std::sync::mpsc::channel();
        buffer.slice(0..bytes).map_async(wgpu::MapMode::Read, move |r| {
            let _ = tx.send(r);
        });
        let _ = self.device.poll(wgpu::PollType::wait_indefinitely());
        rx.recv().expect("gpu map channel").expect("gpu buffer map failed");
    }

    fn read_counter(&self, work: &Work) -> usize {
        self.map_read(&work.read_counter, 4);
        let n = {
            let view = work.read_counter.slice(0..4).get_mapped_range().expect("mapped range");
            cast::<u32>(&view)[0] as usize
        };
        work.read_counter.unmap();
        n
    }

    /// Every survivor of one batch of queries, in their order, at or above
    /// each query's threshold: a superset of the best `max_survivors`, which
    /// [`Self::rank`] then keeps
    fn run_batch(&self, work: &mut Work, batch: &[(&ProcessedSpectrum, KernelParams)]) -> Vec<Vec<Survivor>> {
        let stats = &self.stats;
        stats.batches.fetch_add(1, Ordering::Relaxed);
        let n_queries = batch.len();
        let mut t = Instant::now();

        // each spectrum's peaks and table once: its hypotheses come one after another
        let mut peaks: Vec<[f32; 2]> = Vec::new();
        let mut spectra: Vec<GpuSpectrum> = Vec::new();
        let mut queries = Vec::with_capacity(n_queries);
        let mut last: Option<*const ProcessedSpectrum> = None;
        let mut table_len = 0usize;
        for (spectrum, params) in batch {
            let here: *const ProcessedSpectrum = *spectrum;
            if last != Some(here) {
                let s = GpuSpectrum::new(&spectrum.masses, peaks.len(), table_len);
                table_len += s.n_bins as usize;
                peaks.extend(spectrum.masses.iter().zip(&spectrum.intensities).map(|(&m, &i)| [m, i]));
                spectra.push(s);
                last = Some(here);
            }
            queries.push(GpuQuery::new(spectra.last().expect("a spectrum"), params));
        }
        let st = wgpu::BufferUsages::STORAGE;
        let queries = make_buffer(&self.device, "queries", bytes_of(&queries), st);
        let peaks = make_buffer(&self.device, "peaks", bytes_of(&peaks), st);
        let tables = PeakTables {
            table: self.device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("peak tables"),
                size: (table_len.max(1) * 4) as u64,
                usage: st,
                mapped_at_creation: false,
            }),
            spectra: make_buffer(&self.device, "spectra", bytes_of(&spectra), st),
            count: spectra.len() as u32,
            most_bins: spectra.iter().map(|s| s.n_bins).max().unwrap_or(1),
        };
        let uniforms = self.uniforms(n_queries, spectra.len(), work.capacity);
        t = Stats::add(&stats.setup, t);

        let mut encoder = self.device.create_command_encoder(&Default::default());
        self.encode_scoring(&mut encoder, work, &uniforms, &queries, &peaks, &tables, n_queries);
        self.encode_compact(&mut encoder, work, &uniforms, n_queries);
        self.queue.submit(Some(encoder.finish()));
        let mut found = self.read_counter(work);

        // more survivors than the copy-back buffer holds: a larger one, and
        // the survivors gathered again (the scores are still on the device)
        while found > work.capacity {
            if found > self.max_entries {
                Stats::add(&stats.device, t);
                if n_queries > 1 {
                    let (a, b) = batch.split_at(n_queries / 2);
                    let mut out = self.run_batch(work, a);
                    out.extend(self.run_batch(work, b));
                    return out;
                }
                return vec![self.read_scores(work)];
            }
            Self::resize(&self.device, work, found.next_power_of_two().min(self.max_entries));
            let uniforms = self.uniforms(n_queries, spectra.len(), work.capacity);
            let mut encoder = self.device.create_command_encoder(&Default::default());
            self.encode_compact(&mut encoder, work, &uniforms, n_queries);
            self.queue.submit(Some(encoder.finish()));
            found = self.read_counter(work);
        }

        let bytes = (found * 16) as u64;
        if found > 0 {
            let mut encoder = self.device.create_command_encoder(&Default::default());
            encoder.copy_buffer_to_buffer(&work.entries, 0, &work.read_entries, 0, bytes);
            self.queue.submit(Some(encoder.finish()));
            self.map_read(&work.read_entries, bytes);
        }
        t = Stats::add(&stats.device, t);

        let mut per_query: Vec<Vec<Survivor>> = vec![Vec::new(); n_queries];
        if found > 0 {
            {
                let view = work.read_entries.slice(0..bytes).get_mapped_range().expect("mapped range");
                let entries: &[[u32; 4]] = cast(&view);
                let mut sizes = vec![0usize; n_queries];
                for e in entries {
                    sizes[e[0] as usize] += 1;
                }
                for (survivors, &n) in per_query.iter_mut().zip(&sizes) {
                    survivors.reserve_exact(n);
                }
                for e in entries {
                    per_query[e[0] as usize].push(Survivor {
                        peptide: e[1],
                        score: f32::from_bits(e[2]),
                        site: e[3] >> 16,
                        matched: (e[3] & 0xFFFF) as u16,
                    });
                }
            }
            work.read_entries.unmap();
        }
        Stats::add(&stats.ranking, t);
        per_query
    }

    /// The best `max_survivors` of each query's survivors, in peptide order
    fn rank(&self, survivors: Vec<Vec<Survivor>>, batch: &[(&ProcessedSpectrum, KernelParams)]) -> Vec<Vec<Survivor>> {
        let t = Instant::now();
        let out = survivors
            .into_par_iter()
            .zip(batch.par_iter())
            .map(|(survivors, (_, params))| top_survivors(survivors, params.max_survivors))
            .collect();
        Stats::add(&self.stats.ranking, t);
        out
    }

    /// The survivors of a lone query with more than `max_entries`: every
    /// peptide's score copied back, chunk by chunk
    fn read_scores(&self, work: &Work) -> Vec<Survivor> {
        let t = Instant::now();
        let mut all = Vec::new();
        for (chunk, scores) in self.chunks.iter().zip(&work.scores) {
            let bytes = (chunk.count * 8) as u64;
            let read = self.device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("read_scores"),
                size: bytes,
                usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            });
            let mut encoder = self.device.create_command_encoder(&Default::default());
            encoder.copy_buffer_to_buffer(scores, 0, &read, 0, bytes);
            self.queue.submit(Some(encoder.finish()));
            self.map_read(&read, bytes);
            {
                let view = read.slice(..).get_mapped_range().expect("mapped range");
                let cells: &[[u32; 2]] = cast(&view);
                for (p, cell) in cells.iter().enumerate() {
                    let score = f32::from_bits(cell[0]);
                    if score >= 0.0 {
                        all.push(Survivor {
                            peptide: (chunk.first + p) as u32,
                            score,
                            site: cell[1] >> 16,
                            matched: (cell[1] & 0xFFFF) as u16,
                        });
                    }
                }
            }
            read.unmap();
        }
        Stats::add(&self.stats.device, t);
        all
    }
}

impl Work {
    fn new(device: &wgpu::Device, chunks: &[Chunk], max_batch: usize, capacity: usize) -> Self {
        let buffer = |label: &str, size: usize, usage: wgpu::BufferUsages| {
            device.create_buffer(&wgpu::BufferDescriptor {
                label: Some(label),
                size: size.max(16) as u64,
                usage,
                mapped_at_creation: false,
            })
        };
        let st = wgpu::BufferUsages::STORAGE;
        let mut work = Work {
            // (copied back whole only for a lone query with too many survivors)
            scores: chunks
                .iter()
                .map(|c| buffer("scores", max_batch * c.count * 8, st | wgpu::BufferUsages::COPY_SRC))
                .collect(),
            hist: buffer("hist", max_batch * BINS as usize * 4, st | wgpu::BufferUsages::COPY_DST),
            thresh: buffer("thresh", max_batch * 4, st),
            counter: buffer("counter", 4, st | wgpu::BufferUsages::COPY_SRC | wgpu::BufferUsages::COPY_DST),
            read_counter: buffer("read_counter", 4, wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST),
            entries: buffer("entries", 16, st),
            read_entries: buffer("read_entries", 16, wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST),
            capacity: 0,
        };
        GpuKernel::resize(device, &mut work, capacity);
        work
    }
}

fn make_buffer(device: &wgpu::Device, label: &str, data: &[u8], usage: wgpu::BufferUsages) -> wgpu::Buffer {
    // Storage bindings may not be empty; pad to one word
    let padded;
    let data = if data.is_empty() {
        padded = [0u8; 4];
        &padded[..]
    } else {
        data
    };
    device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some(label),
        contents: data,
        usage,
    })
}

fn build_chunks(device: &wgpu::Device, table: &IonTable, max_binding: usize) -> Result<Vec<Chunk>, String> {
    let n = table.len();
    // an ion takes 8 bytes (mass, meta), a peptide's offsets 8
    let max_ions = max_binding / 8;
    let max_peps = (MAX_WORKGROUPS * WORKGROUP) as usize;
    let mut chunks = Vec::new();
    let mut first = 0usize;
    while first < n {
        let mut end = first + 1;
        let ion_start = table.ion_off[first] as usize;
        if (table.ion_off[end] as usize) - ion_start > max_ions {
            return Err(format!("peptide {first} alone exceeds the GPU buffer limit"));
        }
        while end < n
            && end - first < max_peps
            && (table.ion_off[end + 1] as usize) - ion_start <= max_ions
            && (end + 1 - first + 1) * 8 <= max_binding
        {
            end += 1;
        }
        let (i0, i1) = (table.ion_off[first] as usize, table.ion_off[end] as usize);
        let (s0, s1) = (table.site_off[first] as usize, table.site_off[end] as usize);
        let st = wgpu::BufferUsages::STORAGE;
        let offs: Vec<[u32; 2]> = (first..=end).map(|p| [table.ion_off[p], table.site_off[p]]).collect();
        let ions: Vec<[u32; 2]> = (i0..i1).map(|i| [table.ion_mass[i].to_bits(), table.ion_meta[i]]).collect();
        chunks.push(Chunk {
            first,
            count: end - first,
            ion_base: i0 as u32,
            site_base: s0 as u32,
            offs: make_buffer(device, "offs", bytes_of(&offs), st),
            ions: make_buffer(device, "ions", bytes_of(&ions), st),
            sites: make_buffer(device, "sites", bytes_of(&table.sites[s0..s1]), st),
            pep_mass: make_buffer(device, "pep_mass", bytes_of(&table.peptide_mass[first..end]), st),
        });
        first = end;
    }
    Ok(chunks)
}

impl ChainKernel for GpuKernel {
    fn name(&self) -> &str {
        &self.label
    }

    fn survivors(&self, query: &ProcessedSpectrum, params: &KernelParams) -> Vec<Survivor> {
        self.survivors_many(&[(query, *params)]).pop().unwrap_or_default()
    }

    fn survivors_many(&self, queries: &[(&ProcessedSpectrum, KernelParams)]) -> Vec<Vec<Survivor>> {
        let asked = Instant::now();
        let mut work = self.work.lock().unwrap_or_else(|e| e.into_inner());
        Stats::add(&self.stats.queued, asked);
        self.stats.queries.fetch_add(queries.len() as u64, Ordering::Relaxed);
        let work = &mut *work;
        // the CPU ranks a batch's survivors while the device works on the next
        let mut out = Vec::with_capacity(queries.len());
        let mut pending: Option<(Vec<Vec<Survivor>>, &[(&ProcessedSpectrum, KernelParams)])> = None;
        for batch in queries.chunks(self.max_batch) {
            let previous = pending.take();
            let (survivors, ranked) = rayon::join(
                || self.run_batch(work, batch),
                || previous.map(|(survivors, batch)| self.rank(survivors, batch)),
            );
            out.extend(ranked.into_iter().flatten());
            pending = Some((survivors, batch));
        }
        out.extend(pending.map(|(survivors, batch)| self.rank(survivors, batch)).into_iter().flatten());
        out
    }

    fn batches(&self) -> bool {
        true
    }
}

fn cast<T: Copy>(bytes: &[u8]) -> &[T] {
    // Safety: mapped ranges are 4-byte aligned and sized in whole words
    unsafe { std::slice::from_raw_parts(bytes.as_ptr() as *const T, bytes.len() / std::mem::size_of::<T>()) }
}

#[cfg(test)]
mod pick_test {
    use super::pick_adapter;

    /// `--gpu`: a number, a whole name, or part of one
    #[test]
    fn a_gpu_is_picked_by_number_or_name() {
        let names: Vec<String> = ["NVIDIA T400", "NVIDIA RTX 5000 Ada Generation", "Microsoft Basic Render Driver"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        assert_eq!(pick_adapter(&names, "#2"), Ok(1));
        assert_eq!(pick_adapter(&names, "nvidia rtx 5000 ada generation"), Ok(1));
        assert_eq!(pick_adapter(&names, "RTX 5000"), Ok(1));
        assert_eq!(pick_adapter(&names, "T400"), Ok(0));
        let several = pick_adapter(&names, "NVIDIA").unwrap_err();
        assert!(several.contains("several") && several.contains("#2 NVIDIA RTX 5000"), "{several}");
        assert!(pick_adapter(&names, "#4").unwrap_err().contains("no GPU #4"));
        assert!(pick_adapter(&names, "A100").unwrap_err().contains("no GPU named `A100`"));
        // two identical cards: the name picks the first, #n either
        let twins: Vec<String> = vec!["NVIDIA RTX 5000 Ada Generation".into(), "NVIDIA RTX 5000 Ada Generation".into()];
        assert_eq!(pick_adapter(&twins, "RTX 5000 Ada Generation"), Err(
            "`RTX 5000 Ada Generation` names several GPUs (#1 NVIDIA RTX 5000 Ada Generation; #2 NVIDIA RTX 5000 Ada Generation); give #n or the whole name".into()
        ));
        assert_eq!(pick_adapter(&twins, "NVIDIA RTX 5000 Ada Generation"), Ok(0));
        assert_eq!(pick_adapter(&twins, "#2"), Ok(1));
    }
}
