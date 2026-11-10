//! timsTOF ion mobility (1/K0) from the instrument's own calibration.
//!
//! timsrust converts scan numbers to 1/K0 linearly between the acquisition
//! limits (`OneOverK0AcqRangeLower` / `Upper`). Bruker's software instead
//! applies the calibration stored for each frame in the `TimsCalibration`
//! table of `analysis.tdf`; the two scales differ by up to about 1-2% of 1/K0.
//! The calibrated scale is the default (`bruker_config.ion_mobility_scale`).
//!
//! ModelType 2 is the model timsTOF instruments write. For a scan number `s`:
//!
//! ```text
//! V    = C2 + (C2 - C3) / C1 * (C4 + C0 - s)
//! 1/K0 = V / (C6 * V + C7)
//! ```
//!
//! Outside `C8 <= V <= C9` the model continues as a straight line with the
//! slope at the limit. This reproduces the Bruker SDK's
//! `tims_scannum_to_oneoverk0`, including fractional scans.
//!
//! Adapted from Sage Plus (`crates/sage-cloudpath/src/tims_mobility.rs`, MIT,
//! <https://github.com/pgarrett-scripps/sage-plus>).

use rusqlite::{Connection, OpenFlags};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// Which scan-to-1/K0 scale timsTOF mobilities are reported on
#[derive(Deserialize, Serialize, Debug, Clone, Copy, Default, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum MobilityScale {
    /// Each frame's `TimsCalibration` model, as Bruker's software reports
    #[default]
    Calibrated,
    /// timsrust's interpolation between the acquisition limits (earlier versions)
    Linear,
}

/// One ModelType 2 `TimsCalibration` row, coefficients `C0` to `C9`
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MobilityModel {
    c: [f64; 10],
}

impl MobilityModel {
    pub fn new(coefficients: [f64; 10]) -> Self {
        Self { c: coefficients }
    }

    /// 1/K0 of a (possibly fractional) scan number
    pub fn one_over_k0(&self, scan: f64) -> f64 {
        let [c0, c1, c2, c3, c4, _, c6, c7, c8, c9] = self.c;
        let v = c2 + (c2 - c3) / c1 * (c4 + c0 - scan);
        let model = |v: f64| v / (c6 * v + c7);
        let slope = |v: f64| c7 / (c6 * v + c7).powi(2);
        if v < c8 {
            model(c8) + slope(c8) * (v - c8)
        } else if v > c9 {
            model(c9) + slope(c9) * (v - c9)
        } else {
            model(v)
        }
    }
}

/// The calibration models of one acquisition, the model of each frame, and
/// each DDA precursor's parent frame and fractional scan
#[derive(Debug, Clone)]
pub struct MobilityCalibration {
    models: Vec<MobilityModel>,
    frame_models: HashMap<usize, usize>,
    /// The model used by the most frames, for values without a frame
    dominant: usize,
    precursors: HashMap<usize, (usize, f64)>,
}

/// The `analysis.tdf` of a `.d` folder (or of a file inside it)
pub fn analysis_tdf(path: &Path) -> Option<PathBuf> {
    path.ancestors()
        .filter(|dir| !dir.as_os_str().is_empty())
        .map(|dir| dir.join("analysis.tdf"))
        .find(|tdf| tdf.is_file())
}

impl MobilityCalibration {
    /// Read the calibration of a `.d` folder
    pub fn read(path: &Path) -> Result<Self, String> {
        let tdf = analysis_tdf(path).ok_or("no analysis.tdf")?;
        let con = Connection::open_with_flags(&tdf, OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX)
            .map_err(|e| e.to_string())?;
        let sql = |e: rusqlite::Error| e.to_string();

        let mut rows: Vec<(i64, i64, [f64; 10])> = Vec::new();
        let mut stmt = con
            .prepare("SELECT Id, ModelType, C0, C1, C2, C3, C4, C5, C6, C7, C8, C9 FROM TimsCalibration ORDER BY Id")
            .map_err(sql)?;
        let mut query = stmt.query([]).map_err(sql)?;
        while let Some(row) = query.next().map_err(sql)? {
            let mut c = [0.0; 10];
            for (k, value) in c.iter_mut().enumerate() {
                *value = row.get::<_, f64>(k + 2).map_err(sql)?;
            }
            rows.push((row.get(0).map_err(sql)?, row.get(1).map_err(sql)?, c));
        }
        let mut ids = HashMap::new();
        let mut models = Vec::new();
        for (id, model, c) in rows {
            if model != 2 {
                return Err(format!("TimsCalibration {id} uses ModelType {model}, which XNova does not read"));
            }
            ids.insert(id, models.len());
            models.push(MobilityModel::new(c));
        }
        if models.is_empty() {
            return Err("no TimsCalibration rows".into());
        }

        let mut frame_models = HashMap::new();
        let mut counts = vec![0usize; models.len()];
        let mut stmt = con.prepare("SELECT Id, TimsCalibration FROM Frames").map_err(sql)?;
        let mut query = stmt.query([]).map_err(sql)?;
        while let Some(row) = query.next().map_err(sql)? {
            let (frame, id): (i64, i64) = (row.get(0).map_err(sql)?, row.get(1).map_err(sql)?);
            let model = *ids.get(&id).ok_or(format!("frame {frame} references missing TimsCalibration {id}"))?;
            frame_models.insert(frame as usize, model);
            counts[model] += 1;
        }
        let dominant = (0..models.len()).max_by_key(|&m| (counts[m], std::cmp::Reverse(m))).unwrap_or(0);

        // DDA precursors: parent frame and fractional scan number
        let mut precursors = HashMap::new();
        let has_precursors: bool = con
            .query_row("SELECT COUNT(*) FROM sqlite_master WHERE type = 'table' AND name = 'Precursors'", [], |r| {
                r.get::<_, i64>(0)
            })
            .map_err(sql)?
            > 0;
        if has_precursors {
            let mut stmt = con
                .prepare("SELECT Id, Parent, ScanNumber FROM Precursors WHERE Parent IS NOT NULL AND ScanNumber IS NOT NULL")
                .map_err(sql)?;
            let mut query = stmt.query([]).map_err(sql)?;
            while let Some(row) = query.next().map_err(sql)? {
                let (id, parent, scan): (i64, i64, f64) =
                    (row.get(0).map_err(sql)?, row.get(1).map_err(sql)?, row.get(2).map_err(sql)?);
                precursors.insert(id as usize, (parent as usize, scan));
            }
        }
        Ok(Self { models, frame_models, dominant, precursors })
    }

    /// The model of `frame` (its `Frames.Id`); a frame without one uses the
    /// model of most frames
    pub fn frame(&self, frame: usize) -> &MobilityModel {
        &self.models[self.frame_models.get(&frame).copied().unwrap_or(self.dominant)]
    }

    /// The model used by the most frames
    pub fn dominant(&self) -> &MobilityModel {
        &self.models[self.dominant]
    }

    /// 1/K0 of a DDA precursor (its `Precursors.Id`, recorded in `parent`
    /// frame): the parent frame's model at the precursor's fractional scan
    pub fn precursor(&self, id: usize, parent: usize) -> Option<f64> {
        let &(frame, scan) = self.precursors.get(&id)?;
        (frame == parent).then(|| self.frame(frame).one_over_k0(scan))
    }
}

/// The calibration for a `.d` folder under `scale`: none for the linear
/// scale, or (with a warning) when the folder has no usable calibration
pub fn for_path(path: &Path, scale: MobilityScale) -> Option<Arc<MobilityCalibration>> {
    match scale {
        MobilityScale::Linear => None,
        MobilityScale::Calibrated => match MobilityCalibration::read(path) {
            Ok(calibration) => Some(Arc::new(calibration)),
            Err(e) => {
                log::warn!("{}: ion mobility on the linear scale ({e})", path.display());
                None
            }
        },
    }
}

#[cfg(test)]
mod test {
    use super::*;

    /// The calibration of a timsTOF in 2025 (ModelType 2, one row for all frames)
    const ROW: [f64; 10] = [
        1.0,
        1509.0,
        187.64779662120463,
        78.53435821105788,
        32.72727272727273,
        1.0,
        -0.007497305527082165,
        135.44099330647055,
        13.4184398563705,
        1708.689491089318,
    ];

    #[test]
    fn model_type_2() {
        let model = MobilityModel::new(ROW);
        for (scan, expected) in [(0.0, 1.418389), (750.0, 1.010659), (1500.0, 0.605388)] {
            assert!((model.one_over_k0(scan) - expected).abs() < 1e-6, "scan {scan}: {}", model.one_over_k0(scan));
        }
        // monotone, and continuous at the limits of the model
        let values: Vec<f64> = (-100..1700).map(|s| model.one_over_k0(s as f64)).collect();
        assert!(values.windows(2).all(|w| w[1] < w[0]));
        // fractional scans fall between their neighbours
        let (a, b, f) = (model.one_over_k0(857.0), model.one_over_k0(858.0), model.one_over_k0(857.248));
        assert!(f < a && f > b);
    }

    #[test]
    fn scale_configuration() {
        assert_eq!(serde_json::to_value(MobilityScale::default()).unwrap(), serde_json::json!("calibrated"));
        assert_eq!(serde_json::from_value::<MobilityScale>(serde_json::json!("linear")).unwrap(), MobilityScale::Linear);
    }

    /// A real acquisition: `XNOVA_TDF=<folder>.d`
    #[test]
    #[ignore]
    fn reads_a_real_calibration() {
        let path = PathBuf::from(std::env::var("XNOVA_TDF").expect("set XNOVA_TDF"));
        let calibration = MobilityCalibration::read(&path).unwrap();
        println!(
            "{} models, {} frames, {} precursors; precursor 2: {:?}",
            calibration.models.len(),
            calibration.frame_models.len(),
            calibration.precursors.len(),
            calibration.precursor(2, 1)
        );
        assert!(!calibration.frame_models.is_empty());
        let v = calibration.dominant().one_over_k0(700.0);
        assert!(v > 0.5 && v < 1.6, "{v}");
    }
}
