//! Thermo Fisher `.raw` files, read directly, without converting them to mzML.
//!
//! The file format is decoded by `sage-plus-raw`, a pure-Rust reader (no
//! Thermo libraries; Apache-2.0), Sage Plus's fork of Sigilweaver OpenTFRaw.
//! This module is adapted from Sage Plus's `thermoraw.rs` (MIT,
//! <https://github.com/pgarrett-scripps/sage-plus>).
//!
//! Spectra get the Thermo native ids msconvert writes
//! (`controllerType=0 controllerNumber=1 scan=N`), so results read the same
//! as those of a converted mzML. Each scan's centroid peak list is used (the
//! instrument's own centroids when the scan was recorded in profile mode),
//! and the precursor m/z is Thermo's monoisotopic m/z when it has one, else
//! the isolation target, as msconvert reports them.

use rayon::prelude::*;
use sage_core::{
    mass::Tolerance,
    spectrum::{Precursor, RawSpectrum, Representation},
};
use sage_plus_raw::{scan_metadata, try_extract_spectrum, PrecursorInfo, RawFileReader, SpectrumRecord};
use std::collections::HashSet;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

const ID_PREFIX: &str = "controllerType=0 controllerNumber=1 scan=";

/// Scans decoded by one thread at a time, each with its own file handle
const CHUNK: usize = 256;

/// The id of a scan: its Thermo native id, as msconvert writes it
pub fn scan_id(scan: u32) -> String {
    format!("{ID_PREFIX}{scan}")
}

/// The scan number in a Thermo native id
pub fn scan_number(id: &str) -> Option<u32> {
    id.strip_prefix(ID_PREFIX)?.parse().ok()
}

/// A Thermo RAW file opened for reading: its header, scan index and scan
/// events are read once; spectra are decoded on demand
pub struct ThermoRaw {
    path: PathBuf,
    raw: RawFileReader,
    file_id: usize,
}

impl ThermoRaw {
    pub fn open(path: &Path, file_id: usize) -> sage_plus_raw::Result<Self> {
        Ok(Self { path: path.to_path_buf(), raw: RawFileReader::open_path(path)?, file_id })
    }

    pub fn num_scans(&self) -> u32 {
        self.raw.num_scans
    }

    /// Every scan, in scan order, decoded in parallel. MS1 scans only with
    /// `ms1` (LFQ and SILAC need them, a search does not). A scan that cannot
    /// be decoded is skipped with a warning; it does not fail the file.
    pub fn read_all(&self, ms1: bool) -> sage_plus_raw::Result<Vec<RawSpectrum>> {
        let wanted: Vec<u32> = (0..self.raw.num_scans)
            .filter(|&idx| ms1 || scan_metadata(&self.raw, idx).map_or(true, |m| m.ms_level != 1))
            .collect();
        let decoded = wanted
            .par_chunks(CHUNK)
            .map(|chunk| {
                let mut source = BufReader::new(File::open(&self.path)?);
                Ok(chunk
                    .iter()
                    .map(|&idx| try_extract_spectrum(&self.raw, &mut source, idx, false))
                    .collect::<Vec<_>>())
            })
            .collect::<std::io::Result<Vec<_>>>()?;
        Ok(self.finish(decoded.into_iter().flatten(), wanted.len()))
    }

    /// The scans with these ids (Thermo native ids), without decoding the rest
    /// of the file. Ids that are not scans of this file are left out.
    pub fn read_ids(&self, ids: &HashSet<String>) -> sage_plus_raw::Result<Vec<RawSpectrum>> {
        let first = self.raw.run_header.sample_info.first_scan_number;
        let mut indices: Vec<u32> = ids
            .iter()
            .filter_map(|id| scan_number(id))
            .filter_map(|scan| scan.checked_sub(first))
            .filter(|&idx| idx < self.raw.num_scans)
            .collect();
        indices.sort_unstable();
        indices.dedup();
        let mut source = BufReader::new(File::open(&self.path)?);
        let decoded: Vec<_> = indices
            .iter()
            .map(|&idx| try_extract_spectrum(&self.raw, &mut source, idx, false))
            .collect();
        Ok(self.finish(decoded.into_iter(), indices.len()))
    }

    /// Convert decoded scans, warning once about scans that failed and peaks
    /// that were dropped
    fn finish(
        &self,
        decoded: impl Iterator<Item = sage_plus_raw::Result<SpectrumRecord>>,
        expected: usize,
    ) -> Vec<RawSpectrum> {
        let mut spectra = Vec::with_capacity(expected);
        let (mut failed, mut first_error, mut dropped, mut unsearchable) = (0usize, None, 0usize, 0usize);
        for result in decoded {
            match result {
                Ok(record) => {
                    dropped += record.dropped_peaks;
                    if record.ms_level == 2 && record.precursor.as_ref().and_then(plausible_mz).is_none() {
                        unsearchable += 1;
                    }
                    spectra.push(convert(record, self.file_id));
                }
                Err(error) => {
                    failed += 1;
                    first_error.get_or_insert(error);
                }
            }
        }
        let path = self.path.display();
        if let Some(error) = first_error {
            log::warn!("skipped {failed} of {expected} scans in {path} that could not be decoded; first error: {error}");
        }
        if dropped > 0 {
            log::warn!("dropped {dropped} peaks with an invalid m/z or intensity from {path}");
        }
        if unsearchable > 0 {
            log::warn!(
                "{unsearchable} MS2 scans in {path} have no plausible precursor m/z and are not searched; \
                 convert the file to mzML with msconvert to search them"
            );
        }
        spectra
    }
}

/// Read a whole Thermo RAW file
pub fn read(path: &Path, file_id: usize, ms1: bool) -> sage_plus_raw::Result<Vec<RawSpectrum>> {
    ThermoRaw::open(path, file_id)?.read_all(ms1)
}

fn convert(record: SpectrumRecord, file_id: usize) -> RawSpectrum {
    let ms_level = record.ms_level;
    let precursor = record.precursor.and_then(|value| {
        // SPS-MS3 reporter-ion quantification needs only the link to the
        // MS2 scan, which is kept even when the precursor m/z is unknown
        let mz = match plausible_mz(&value) {
            Some(mz) => mz as f32,
            None if ms_level > 2 && value.master_scan_number.is_some() => 0.0,
            None => return None,
        };
        Some(Precursor {
            mz,
            charge: value.charge.and_then(|charge| u8::try_from(charge).ok()).filter(|&z| z > 0),
            spectrum_ref: value.master_scan_number.map(scan_id),
            isolation_window: value.isolation_width.map(|width| {
                let half = width as f32 / 2.0;
                Tolerance::Da(-half, half)
            }),
            ..Default::default()
        })
    });
    RawSpectrum {
        file_id,
        ms_level: ms_level.min(u8::MAX as u32) as u8,
        id: scan_id(record.scan_number),
        precursors: precursor.into_iter().collect(),
        // decoded with `include_profile: false`: every scan is its centroid list
        representation: Representation::Centroid,
        scan_start_time: record.retention_time_min as f32,
        ion_injection_time: record.ion_injection_time_ms.unwrap_or_default() as f32,
        total_ion_current: record.total_ion_current as f32,
        mz: record.mz.into_iter().map(|mz| mz as f32).collect(),
        intensity: record.intensity,
        mobility: None,
    }
}

fn plausible_precursor_mz(mz: f64) -> bool {
    (50.0..20_000.0).contains(&mz)
}

/// The monoisotopic (selected) m/z, or else the isolation target, when plausible
fn plausible_mz(precursor: &PrecursorInfo) -> Option<f64> {
    [precursor.selected_mz, precursor.target_mz]
        .into_iter()
        .flatten()
        .find(|&mz| plausible_precursor_mz(mz))
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn native_ids() {
        assert_eq!(scan_id(1234), "controllerType=0 controllerNumber=1 scan=1234");
        assert_eq!(scan_number("controllerType=0 controllerNumber=1 scan=1234"), Some(1234));
        assert_eq!(scan_number("scan=1234"), None);
        assert_eq!(scan_number("controllerType=0 controllerNumber=1 scan=x"), None);
    }

    #[test]
    fn precursor_mz_falls_back_to_the_isolation_target() {
        let p = |selected: Option<f64>, target: Option<f64>| PrecursorInfo {
            target_mz: target,
            selected_mz: selected,
            ..Default::default()
        };
        assert_eq!(plausible_mz(&p(Some(722.3242), Some(722.5))), Some(722.3242));
        assert_eq!(plausible_mz(&p(Some(0.0), Some(722.5))), Some(722.5));
        assert_eq!(plausible_mz(&p(None, Some(1.0))), None);
    }

    /// ProteoWizard's Thermo test files read directly and as ProteoWizard
    /// wrote them to mzML give the same spectrum. The files are not part of
    /// the repository: set `XNOVA_THERMO_FIXTURES` to a folder with
    /// `BSA-FT-HCD.raw`, `BSA-FT-HCD-centroid.mzML`, `BSA-FT-ETD.raw` and
    /// `BSA-FT-ETD-centroid.mzML` from github.com/ProteoWizard/pwiz
    /// (pwiz/data/vendor_readers/Thermo/Reader_Thermo_Test.data).
    #[test]
    #[ignore]
    fn same_spectra_as_proteowizard() {
        let dir = PathBuf::from(std::env::var("XNOVA_THERMO_FIXTURES").expect("set XNOVA_THERMO_FIXTURES"));
        for name in ["BSA-FT-HCD", "BSA-FT-ETD"] {
            let raw_path = dir.join(format!("{name}.raw"));
            let url = url::Url::from_file_path(dir.join(format!("{name}-centroid.mzML"))).unwrap();
            let expected = crate::util::read_mzml(&url, 0, None).unwrap();
            // through the same entry point as the search
            let raw_url = url::Url::from_file_path(&raw_path).unwrap();
            let found = crate::util::read_spectra(&raw_url, 0, None, Default::default(), true).unwrap();
            assert_eq!(found.len(), expected.len(), "{name}: spectra");
            for (a, b) in found.iter().zip(&expected) {
                assert_eq!((&a.id, a.ms_level), (&b.id, b.ms_level), "{name}");
                assert!((a.scan_start_time - b.scan_start_time).abs() < 1e-4, "{name}: rt {} {}", a.scan_start_time, b.scan_start_time);
                assert_eq!(a.precursors.len(), b.precursors.len(), "{name}: precursors");
                for (p, q) in a.precursors.iter().zip(&b.precursors) {
                    assert!((p.mz - q.mz).abs() < 1e-4, "{name}: precursor m/z {} {}", p.mz, q.mz);
                    assert_eq!(p.charge, q.charge, "{name}: charge");
                    assert_eq!(p.isolation_window, q.isolation_window, "{name}: isolation window");
                }
                assert_eq!(a.mz.len(), b.mz.len(), "{name}: peaks");
                for (x, y) in a.mz.iter().zip(&b.mz) {
                    assert!((x - y).abs() / y < 1e-6, "{name}: m/z {x} {y}");
                }
                for (x, y) in a.intensity.iter().zip(&b.intensity) {
                    assert!((x - y).abs() <= 1e-3 * y.abs().max(1.0), "{name}: intensity {x} {y}");
                }
                println!(
                    "{name}: {} MS{} precursor {:.4} {:?}+, {} peaks: same as ProteoWizard's mzML",
                    a.id, a.ms_level, a.precursors[0].mz, a.precursors[0].charge, a.mz.len()
                );
            }
            // a scan read on its own is the same as in the whole file
            let run = ThermoRaw::open(&raw_path, 0).unwrap();
            let ids: HashSet<String> = found.iter().map(|s| s.id.clone()).collect();
            let one = run.read_ids(&ids).unwrap();
            assert_eq!(one.len(), found.len());
            assert_eq!((&one[0].id, &one[0].mz), (&found[0].id, &found[0].mz));
            assert!(run.read_ids(&HashSet::from(["controllerType=0 controllerNumber=1 scan=999999".to_string()])).unwrap().is_empty());
        }
    }
}
