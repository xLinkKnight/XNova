//! `lfq_long.tsv`: label-free quantification in long form, one row per
//! precursor and run (adapted from Sage Plus, MIT), next to the wide
//! `lfq.tsv`. Each row says whether the run has its own MS2 identification
//! of the precursor (`ms2_confirmed`; otherwise the intensity was found by
//! matching between runs), where that run's peak is (apex, integration
//! bounds and FWHM in the run's own retention times), how well its isotopes
//! match (`spectral_angle`), and a row-level extraction q-value from the
//! paired decoy in the same run (`extraction_q`; a diagnostic, not a
//! calibrated match-between-runs FDR). `precursor_q` is `lfq.tsv`'s q-value.

use super::Runner;
use sage_cloudpath::Url;
use sage_core::database::PeptideIx;
use sage_core::lfq::{Peak, PrecursorId};
use sage_core::ml::retention_alignment::Alignment;
use sage_core::scoring::Feature;
use std::collections::{HashMap, HashSet};

/// PSM q-value for a run's own identification of a precursor
const MS2_CONFIRMED_Q: f32 = 0.01;

impl Runner {
    pub(super) fn write_lfq_long(
        &self,
        areas: &HashMap<(PrecursorId, bool), (Peak, Vec<f64>), fnv::FnvBuildHasher>,
        features: &[Feature],
        filenames: &[String],
        alignments: &[Alignment],
    ) -> anyhow::Result<Url> {
        let peptide_q = self.parameters.quant.lfq_settings.peptide_q_value;
        // (peptide, charge or none for combined charge states, run)
        let mut confirmed: HashSet<(PeptideIx, Option<u8>, usize)> = HashSet::new();
        for f in features
            .iter()
            .filter(|f| f.label == 1 && f.rank == 1 && f.spectrum_q <= MS2_CONFIRMED_Q && f.peptide_q <= peptide_q)
        {
            confirmed.insert((f.peptide_idx, None, f.file_id));
            confirmed.insert((f.peptide_idx, Some(f.charge), f.file_id));
        }

        let mut precursors: Vec<(String, i32, PeptideIx, Option<u8>, &Peak, &Vec<f64>)> = areas
            .iter()
            .filter(|((_, decoy), _)| !decoy)
            .map(|((id, _), (peak, data))| {
                let (ix, charge) = match *id {
                    PrecursorId::Combined(ix) => (ix, None),
                    PrecursorId::Charged((ix, charge)) => (ix, Some(charge)),
                };
                let name = self.database[ix].display(&self.database.modifications).to_string();
                (name, charge.map_or(-1, i32::from), ix, charge, peak, data)
            })
            .collect();
        precursors.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

        let mut wtr = csv::WriterBuilder::new().delimiter(b'\t').from_writer(vec![]);
        wtr.write_record([
            "peptide",
            "charge",
            "proteins",
            "filename",
            "intensity",
            "ms2_confirmed",
            "apex_rt",
            "rt_start",
            "rt_end",
            "fwhm",
            "aligned_apex_rt",
            "spectral_angle",
            "extraction_score",
            "extraction_q",
            "precursor_q",
        ])?;
        let float = |x: f64| ryu::Buffer::new().format(x).to_string();
        let rt = |x: f32| format!("{x:.4}");
        for (name, charge, ix, charge_id, peak, data) in precursors {
            let proteins = self.database[ix].proteins(&self.database.decoy_tag, self.database.generate_decoys);
            for (file, filename) in filenames.iter().enumerate() {
                let Some(run) = peak.files.get(file) else { continue };
                let alignment = &alignments[file];
                let own = |aligned: f32| alignment.inverse(aligned);
                let ms2 = confirmed.contains(&(ix, charge_id, file));
                wtr.write_record([
                    name.clone(),
                    charge.to_string(),
                    proteins.clone(),
                    filename.clone(),
                    float(data.get(file).copied().unwrap_or(0.0)),
                    u8::from(ms2).to_string(),
                    rt(own(run.apex_rt)),
                    rt(own(run.start_rt)),
                    rt(own(run.end_rt)),
                    rt((own(run.half_end_rt) - own(run.half_start_rt)).max(0.0)),
                    format!("{:.5}", run.apex_rt),
                    format!("{:.4}", run.spectral_angle),
                    format!("{:.4}", run.score),
                    ryu::Buffer::new().format(run.q_value).to_string(),
                    ryu::Buffer::new().format(peak.q_value).to_string(),
                ])?;
            }
        }
        let path = self.make_path("lfq_long.tsv");
        sage_cloudpath::write_bytes_sync(&path, wtr.into_inner()?)?;
        Ok(path)
    }
}
