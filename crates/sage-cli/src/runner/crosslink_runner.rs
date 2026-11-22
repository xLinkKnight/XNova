//! Crosslinking (XL-MS) search workflow: read spectra, search for crosslinked
//! peptide pairs, apply target/decoy FDR, and write `crosslinks.sage.tsv`.
//!
//! This is a child module of `runner` so that it can reuse `Runner`'s private
//! spectrum-reading and path helpers.

use super::Runner;
use crate::output::SageResults;
use crate::telemetry;
use log::info;
use rayon::prelude::*;
use sage_core::crosslink::{
    assign_q_values, assign_q_values_in_window, CrosslinkMatch, CrosslinkScorer, CrosslinkSettings, LinkKind,
    PairClass,
};
use sage_core::crosslink_quant::{quantify_lfq, CrosslinkLfq};
use sage_core::mass::Tolerance;
use sage_core::scoring::Scorer;
use sage_core::spectrum::ProcessedSpectrum;
use std::collections::HashMap;
use std::sync::atomic::Ordering;
use std::time::Instant;

/// Isobaric reporter intensities of one MS2 spectrum: (ion injection time,
/// reporter intensities), keyed by (file id, MS2 spectrum id)
/// CSM-level q-value cutoff of the filtered report
const FDR_CUTOFF: f32 = 0.01;
/// Target-target CSMs that pass [`FDR_CUTOFF`], same columns as `crosslinks.sage.tsv`
const FDR_FILTERED: &str = "crosslinks_FDRfiltered.sage.tsv";
/// The proteins the linear search found, which the crosslink search was
/// narrowed to (`linear_proteins`), with their FASTA headers
pub const LINEAR_PROTEINS: &str = "LinearProteins.fasta";

type Reporters = HashMap<(usize, String), (f32, Vec<f32>)>;

/// The matches at 1% FDR in each FDR group (see `sage_core::crosslink`'s
/// `fdr_group`): targets, target-decoy and decoy-decoy matches, as the FDR
/// counted them, so the log shows where the decoys fell
fn log_fdr_groups(matches: &[CrosslinkMatch]) {
    let kinds = ["inter-protein", "intra-protein", "loop-links", "multivalent"];
    // [kind][recorded, guessed charge][TT, TD, DD]
    let mut counts = [[[0usize; 3]; 2]; 4];
    for m in matches.iter().filter(|m| m.spectrum_q <= FDR_CUTOFF) {
        let kind = match (m.kind, m.intra_protein) {
            (LinkKind::Inter, false) => 0,
            (LinkKind::Inter, true) => 1,
            (LinkKind::Loop, _) => 2,
            (LinkKind::Multi, _) => 3,
        };
        let class = match m.class {
            PairClass::TT => 0,
            PairClass::TD => 1,
            PairClass::DD => 2,
        };
        counts[kind][m.charge_guessed as usize][class] += 1;
    }
    let rows: Vec<(String, [usize; 3])> = (0..2)
        .flat_map(|guessed| {
            (0..4).map(move |k| {
                let name = if guessed == 1 { format!("{}, guessed charge", kinds[k]) } else { kinds[k].to_string() };
                (name, counts[k][guessed])
            })
        })
        .filter(|(_, c)| c.iter().sum::<usize>() > 0)
        .collect();
    if rows.is_empty() {
        return;
    }
    info!("crosslink FDR, matches at 1% by group: targets / target-decoy / decoy-decoy");
    let t = crate::progress::thousands;
    for (name, [tt, td, dd]) in rows {
        info!("  {name:<32} {:>7} / {:>5} / {:>5}", t(tt), t(td), t(dd));
    }
}

/// `link_type` column: `crosslink`, `loop` or `multivalent`
fn link_type(kind: LinkKind) -> &'static str {
    match kind {
        LinkKind::Inter => "crosslink",
        LinkKind::Loop => "loop",
        LinkKind::Multi => "multivalent",
    }
}

/// `database_a` / `database_b` columns: `1`, `2`, or `1+2` for a peptide of
/// proteins in both crosslink databases
fn database_name(group: u8) -> &'static str {
    use sage_core::database::{DATABASE_1, DATABASE_2};
    match (group & DATABASE_1 != 0, group & DATABASE_2 != 0) {
        (true, true) => "1+2",
        (true, false) => "1",
        (false, true) => "2",
        (false, false) => "",
    }
}

/// `topology` column: `a1b` (two peptides), `a1a` (two copies of one
/// peptide), `loop`, or a multivalent species' shape (`a2b`, `a1b1c`, …)
fn topology(m: &CrosslinkMatch) -> String {
    match (m.kind, &m.species, m.beta) {
        (LinkKind::Multi, Some(s), _) => s.topology(),
        (LinkKind::Loop, ..) => "loop".into(),
        (_, _, Some(b)) if b == m.alpha => "a1a".into(),
        _ => "a1b".into(),
    }
}

/// `species` column of a multivalent species: its peptides, numbered, and
/// each link by peptide and residue (1-based):
/// `[1] PEPKTIDE [2] KPEPTIDE | 1:K4-2:K1; …`
fn species_text(db: &sage_core::database::IndexedDatabase, s: &sage_core::crosslink::Species) -> String {
    let peptides: Vec<String> = s
        .peptides
        .iter()
        .enumerate()
        .map(|(i, &ix)| format!("[{}] {}", i + 1, db[ix].display(&db.modifications)))
        .collect();
    let residue = |p: u8, site: u16| {
        let seq = &db[s.peptides[p as usize]].sequence;
        let aa = seq.get(site as usize).map_or('?', |&c| c as char);
        format!("{}:{aa}{}", p + 1, site + 1)
    };
    let links: Vec<String> = s.links.iter().map(|&(a, sa, b, sb)| format!("{}-{}", residue(a, sa), residue(b, sb))).collect();
    format!("{} | {}", peptides.join(" "), links.join("; "))
}

/// The `abs_pos_` columns' letters: `a` and `b`, and for a multivalent search
/// one per peptide a species may have (the valency + 1)
fn abs_pos_letters(settings: &CrosslinkSettings) -> Vec<char> {
    let peptides = settings.multivalent.map_or(2, |v| (v as usize + 1).max(2));
    (0..peptides).map(|i| (b'a' + i as u8) as char).collect()
}

/// The `abs_pos_` cells of a match: where its link residues are in their
/// proteins (1-based, the proteins' own numbering). A pair: `abs_pos_a` for
/// `link_a` in `proteins_a`, `abs_pos_b` for `link_b` in `proteins_b`. A
/// loop-link: both in `proteins_a`. A multivalent species: `abs_pos_a`,
/// `abs_pos_b`, `abs_pos_c`, ... for the peptides `[1]`, `[2]`, `[3]`, ... of
/// its `species` column, each with all its link residues. Several proteins
/// are separated by `;`, several residues of one peptide by `,`.
fn absolute_positions(ctx: &super::crosslink_export::ExportContext, m: &CrosslinkMatch, columns: usize) -> Vec<String> {
    let mut cells: Vec<String> = match (m.kind, m.beta) {
        (LinkKind::Multi, _) => super::crosslink_report::members(m)
            .iter()
            .map(|(ix, sites)| ctx.absolute_positions(*ix, sites))
            .collect(),
        (LinkKind::Loop, _) | (_, None) => vec![
            ctx.absolute_positions(m.alpha, &[m.alpha_site]),
            ctx.absolute_positions(m.alpha, &[m.beta_site]),
        ],
        (_, Some(beta)) => vec![
            ctx.absolute_positions(m.alpha, &[m.alpha_site]),
            ctx.absolute_positions(beta, &[m.beta_site]),
        ],
    };
    cells.resize(columns, String::new());
    cells
}

/// Narrow `fasta` to the proteins in `keep` and, with decoys from the FASTA
/// (not generated), their decoys; each protein once (a protein repeated in
/// the FASTA, or in both crosslink databases, whose database groups still say
/// so). Returns `(kept, total)`: the distinct target proteins kept, and in
/// all. A decoy is an accession containing `decoy_tag`, as Sage reads them.
pub(super) fn keep_proteins(
    fasta: &mut sage_core::fasta::Fasta,
    keep: &std::collections::HashSet<std::sync::Arc<str>>,
    decoy_tag: &str,
    generate_decoys: bool,
) -> (usize, usize) {
    let is_decoy = |acc: &str| !generate_decoys && !decoy_tag.is_empty() && acc.contains(decoy_tag);
    let total = fasta.targets.iter().filter(|(acc, _)| !is_decoy(acc)).map(|(acc, _)| acc.clone()).collect::<std::collections::HashSet<_>>().len();
    let mut seen = std::collections::HashSet::new();
    fasta.targets.retain(|(acc, _)| {
        let wanted = if is_decoy(acc) { keep.contains(acc.replacen(decoy_tag, "", 1).as_str()) } else { keep.contains(acc.as_ref()) };
        wanted && seen.insert(acc.clone())
    });
    let kept = fasta.targets.iter().filter(|(acc, _)| !is_decoy(acc)).count();
    (kept, total)
}

/// The crosslink results in numbered parts. Each is announced in the log
/// (`crosslink results 7 of 9: crosslinks.mzid for xiVIEW and PRIDE (108,693
/// matches)`) and ends with its time and rate (`- crosslinks.mzid for xiVIEW
/// and PRIDE: 1.2 s, 90,578 matches/s`); one that runs long has `progress:`
/// lines with its share done and the time left. The Progress tab follows them.
struct Parts {
    n: usize,
    i: usize,
}

struct Part {
    name: String,
    items: usize,
    unit: &'static str,
    start: Instant,
    tracker: crate::progress::Tracker,
}

impl Parts {
    /// Begin the next part, of `items` `unit`s; `counted`: its code counts
    /// them as it goes (otherwise its progress lines give the time only)
    fn start(&mut self, name: &str, items: usize, unit: &'static str, counted: bool) -> Part {
        self.i += 1;
        let label = format!("crosslink results {} of {}: {name}", self.i, self.n.max(self.i));
        let size = if items > 0 { format!(" ({} {unit})", crate::progress::thousands(items)) } else { String::new() };
        info!("{label}{size}");
        let tracker = crate::progress::step(label, if counted { items } else { 0 }, unit);
        Part { name: name.to_string(), items, unit, start: Instant::now(), tracker }
    }
}

impl Part {
    fn counter(&self) -> crate::progress::Counter {
        self.tracker.counter()
    }

    fn add(&self, n: usize) {
        self.tracker.add(n);
    }

    fn detail(&self, text: impl Into<String>) {
        self.tracker.detail(text);
    }

    /// The part is over: its time, and its rate when it had items
    fn done(self) {
        let took = self.start.elapsed();
        let secs = took.as_secs_f64();
        let rate = if self.items > 0 && secs >= 0.001 {
            format!(", {} {}/s", crate::progress::thousands((self.items as f64 / secs).round() as usize), self.unit)
        } else {
            String::new()
        };
        info!("- {}: {}{rate}", self.name, crate::progress::duration(took));
    }
}

/// What a crosslink search keeps from each chunk of spectrum files
#[derive(Default)]
pub(super) struct CrosslinkChunks {
    /// Candidate crosslink-spectrum matches
    matches: Vec<CrosslinkMatch>,
    /// MS1 spectra, for crosslink LFQ and SILAC ratios
    ms1_spectra: Vec<ProcessedSpectrum>,
    /// Isobaric reporter intensities of the MS2 spectra
    reporters: Reporters,
    /// MS2 spectra searched for crosslinks
    searched: usize,
}

impl Runner {
    /// The search's crosslink settings, with its isotope labels
    fn crosslink_settings(&self) -> CrosslinkSettings {
        let mut settings = self
            .parameters
            .crosslink
            .clone()
            .expect("a crosslink search has crosslink settings");
        // isotope labels (SILAC, uniformly labelled proteins) leave a residue free to link
        settings.labels = self.database.isotope_mods();
        settings
    }

    /// The scorer crosslink chains are scored with
    fn crosslink_base_scorer(&self) -> Scorer<'_> {
        Scorer {
            db: &self.database,
            precursor_tol: self.parameters.precursor_tol,
            fragment_tol: self.parameters.fragment_tol,
            min_matched_peaks: self.parameters.min_matched_peaks,
            min_isotope_err: self.parameters.isotope_errors.0,
            max_isotope_err: self.parameters.isotope_errors.1,
            min_precursor_charge: self.parameters.precursor_charge.0,
            max_precursor_charge: self.parameters.precursor_charge.1,
            override_precursor_charge: self.parameters.override_precursor_charge,
            max_fragment_charge: self.parameters.max_fragment_charge,
            chimera: false,
            report_psms: self.parameters.report_psms,
            wide_window: false,
            // a spectral library needs the matched fragments
            annotate_matches: self.parameters.spectral_library.is_some(),
            score_type: self.parameters.score_type,
        }
    }

    /// Whether the crosslink search keeps the MS1 spectra (LFQ, SILAC)
    fn crosslink_needs_ms1(&self) -> bool {
        self.parameters.quant.lfq || self.database.silac.is_some()
    }

    /// Search one chunk's MS2 spectra for crosslinks, adding the candidate
    /// matches and TMT reporters to `found`
    fn search_crosslink_chunk(
        &self,
        xl: &CrosslinkScorer,
        settings: &CrosslinkSettings,
        msn: &[ProcessedSpectrum],
        found: &mut CrosslinkChunks,
    ) {
        if let Some(isobaric) = &self.parameters.quant.tmt {
            // MS3 reporters are keyed by their MS2 precursor scan, like the CSMs
            let level = self.parameters.quant.tmt_settings.level;
            let quant = sage_core::tmt::quantify(msn, isobaric, Tolerance::Ppm(-20.0, 20.0), level);
            found.reporters.extend(
                quant
                    .into_iter()
                    .map(|q| ((q.file_id, q.spec_id), (q.ion_injection_time, q.peaks))),
            );
        }

        let start = Instant::now();
        let searched = msn
            .par_iter()
            .filter(|spec| spec.level == 2 && spec.masses.len() >= self.parameters.min_peaks)
            .count();
        let searching = crate::progress::step(format!("crosslink search, {}", self.files_of(msn)), searched, "spectra");
        let matches: Vec<CrosslinkMatch> = msn
            .par_iter()
            .filter(|spec| spec.level == 2 && spec.masses.len() >= self.parameters.min_peaks)
            .flat_map(|spec| {
                let found = xl.score(spec);
                searching.add(1);
                found
            })
            .collect();
        if let Ok(mut t) = self.timing.lock() {
            t.crosslink += start.elapsed();
        }
        // scans without a precursor are skipped (the linear search, when it
        // runs too, has said so already)
        if !settings.linear_search {
            self.warn_no_precursor(msn);
        }
        info!(
            "- crosslink search: {:8} ms ({} spectra, {} candidate matches)",
            start.elapsed().as_millis(),
            searched,
            matches.len()
        );
        found.searched += searched;
        found.matches.extend(matches);
    }

    pub(super) fn run_crosslink(
        self,
        parallel: usize,
        parquet: bool,
    ) -> anyhow::Result<telemetry::Telemetry> {
        if parquet {
            log::warn!("parquet output is not available for crosslink searches, writing tsv");
        }
        let settings = self.crosslink_settings();
        let batch_size = parallel.max(1);
        info!("=== searching the spectrum files ===");
        let found = {
            let scorer = self.crosslink_base_scorer();
            let kernel = sage_core::crosslink::make_kernel(&settings, &self.database)
                .map_err(anyhow::Error::msg)?;
            let xl = CrosslinkScorer::new(&scorer, &settings).with_kernel(kernel.as_deref());
            if let Some(kernel) = &kernel {
                info!("crosslink chain scoring backend: {}", kernel.name());
            }
            let mut found = CrosslinkChunks::default();
            for (chunk_idx, chunk) in self.parameters.mzml_paths.chunks(batch_size).enumerate() {
                let (ms1, msn) = self.read_processed_spectra(chunk, chunk_idx, batch_size);
                if self.crosslink_needs_ms1() {
                    found.ms1_spectra.extend(ms1);
                }
                self.search_crosslink_chunk(&xl, &settings, &msn, &mut found);
            }
            found
        };
        self.finish_crosslink(settings, found)
    }

    /// A crosslink search with the linear peptide search: each chunk of files
    /// is read once and searched both ways. Reading the files is most of a
    /// timsTOF search, and the two searches used to read every file each.
    /// The linear results are made first, as before: the crosslink spectral
    /// library takes the linear library's scales, and the crosslinks'
    /// entrapment joins the linear search's.
    pub(super) fn run_linear_and_crosslink(
        mut self,
        parallel: usize,
        parquet: bool,
    ) -> anyhow::Result<telemetry::Telemetry> {
        if self.parameters.crosslink.as_ref().is_some_and(|c| c.linear_proteins) {
            return self.run_linear_then_crosslink(parallel, parquet);
        }
        info!("linear peptide search (results.sage.tsv) and crosslink search, reading each spectrum file once");
        if parquet {
            log::warn!("parquet output is not available for crosslink searches, writing tsv");
        }
        let settings = self.crosslink_settings();
        let batch_size = parallel.max(1);
        info!("=== searching the spectrum files: linear and crosslink search ===");
        let (outputs, mut found) = {
            let linear = self.linear_scorer();
            let scorer = self.crosslink_base_scorer();
            let kernel = sage_core::crosslink::make_kernel(&settings, &self.database)
                .map_err(anyhow::Error::msg)?;
            let xl = CrosslinkScorer::new(&scorer, &settings).with_kernel(kernel.as_deref());
            if let Some(kernel) = &kernel {
                info!("crosslink chain scoring backend: {}", kernel.name());
            }
            let mut found = CrosslinkChunks::default();
            let mut results = Vec::new();
            for (chunk_idx, chunk) in self.parameters.mzml_paths.chunks(batch_size).enumerate() {
                let (ms1, msn) = self.read_processed_spectra(chunk, chunk_idx, batch_size);
                let features = self.search_processed_spectra(&linear, &msn);
                self.search_crosslink_chunk(&xl, &settings, &msn, &mut found);
                results.push(self.complete_features(msn, ms1, features));
            }
            (results.into_iter().collect::<SageResults>(), found)
        };
        // the linear search's MS1 spectra are the crosslink search's
        let ms1 = self.finish_linear(outputs, parquet)?;
        if self.crosslink_needs_ms1() {
            found.ms1_spectra = ms1;
        }
        self.finish_crosslink(settings, found)
    }

    /// A crosslink search narrowed to the proteins of the linear search
    /// (`linear_proteins`): the linear search of every file and its results
    /// first, then a database of the proteins it found (LinearProteins.fasta),
    /// then the crosslink search of every file in it. Each file is read
    /// twice: the crosslink database is only known once every file has been
    /// searched linearly.
    fn run_linear_then_crosslink(mut self, parallel: usize, parquet: bool) -> anyhow::Result<telemetry::Telemetry> {
        info!("linear peptide search (results.sage.tsv) of every file first; the proteins it finds make the crosslink search's database");
        if parquet {
            log::warn!("parquet output is not available for crosslink searches, writing tsv");
        }
        let batch_size = parallel.max(1);
        info!("=== searching the spectrum files: linear search ===");
        let outputs = {
            let linear = self.linear_scorer();
            let mut results = Vec::new();
            for (chunk_idx, chunk) in self.parameters.mzml_paths.chunks(batch_size).enumerate() {
                let (ms1, msn) = self.read_processed_spectra(chunk, chunk_idx, batch_size);
                let features = self.search_processed_spectra(&linear, &msn);
                results.push(self.complete_features(msn, ms1, features));
            }
            results.into_iter().collect::<SageResults>()
        };
        let ms1 = self.finish_linear(outputs, parquet)?;
        self.narrow_to_linear_proteins()?;

        // the settings' isotope labels are the new database's
        let settings = self.crosslink_settings();
        info!("=== searching the spectrum files: crosslink search ===");
        let mut found = {
            let scorer = self.crosslink_base_scorer();
            let kernel = sage_core::crosslink::make_kernel(&settings, &self.database).map_err(anyhow::Error::msg)?;
            let xl = CrosslinkScorer::new(&scorer, &settings).with_kernel(kernel.as_deref());
            if let Some(kernel) = &kernel {
                info!("crosslink chain scoring backend: {}", kernel.name());
            }
            let mut found = CrosslinkChunks::default();
            for (chunk_idx, chunk) in self.parameters.mzml_paths.chunks(batch_size).enumerate() {
                let (_, msn) = self.read_processed_spectra(chunk, chunk_idx, batch_size);
                self.search_crosslink_chunk(&xl, &settings, &msn, &mut found);
            }
            found
        };
        // the linear search's MS1 spectra are the crosslink search's
        if self.crosslink_needs_ms1() {
            found.ms1_spectra = ms1;
        }
        self.finish_crosslink(settings, found)
    }

    /// Rebuild the database from the proteins the linear search found, and
    /// write them to LinearProteins.fasta with their FASTA headers. With no
    /// protein found, the whole database stays.
    fn narrow_to_linear_proteins(&mut self) -> anyhow::Result<()> {
        let keep = self.linear_proteins.take().unwrap_or_default();
        let fdr = self.parameters.crosslink.as_ref().map_or(0.01, |c| c.linear_protein_fdr);
        let percent = crate::summary::percent_text(fdr);
        info!("=== crosslink database: the proteins of the linear search ===");
        let mut fasta = super::read_protein_databases(&mut self.parameters.database)?;
        if keep.is_empty() {
            let total = keep_proteins(&mut fasta.clone(), &keep, &self.parameters.database.decoy_tag, self.parameters.database.generate_decoys).1;
            log::warn!(
                "the linear search found no peptide at {percent} FDR: the crosslink search uses the whole database ({} proteins)",
                crate::progress::thousands(total)
            );
            return Ok(());
        }
        let (kept, total) =
            keep_proteins(&mut fasta, &keep, &self.parameters.database.decoy_tag, self.parameters.database.generate_decoys);

        // LinearProteins.fasta: their FASTA records as they were
        let wanted: std::collections::BTreeSet<&str> = keep.iter().map(|a| a.as_ref()).collect();
        let mut written = std::collections::HashSet::new();
        let mut text = String::new();
        for path in self.parameters.database.fasta_files() {
            let source = sage_cloudpath::util::read_text(&sage_cloudpath::to_url(path)?)?;
            text += &super::xlink_proteins::select_records(&source, &wanted, &mut written);
        }
        let path = self.make_path(LINEAR_PROTEINS);
        sage_cloudpath::write_bytes_sync(&path, text.into_bytes())?;
        self.parameters.output_paths.push(path);
        info!(
            "crosslink database: {} of the {} proteins have a peptide at {percent} FDR in the linear search; written to {LINEAR_PROTEINS}",
            crate::progress::thousands(kept),
            crate::progress::thousands(total)
        );

        let start = Instant::now();
        let _step = crate::progress::step("building the crosslink database of the linear search's proteins", 0, "");
        self.database = self.parameters.database.clone().build(fasta);
        info!(
            "crosslink database: {} peptides, {} fragments, built in {}",
            crate::progress::thousands(self.database.peptides.len()),
            crate::progress::thousands(self.database.fragments.len()),
            crate::progress::duration(start.elapsed())
        );
        Ok(())
    }

    /// The crosslink search's results, once every file is searched: FDR,
    /// quantification and the crosslink outputs, in numbered [`Parts`]
    fn finish_crosslink(
        mut self,
        settings: CrosslinkSettings,
        found: CrosslinkChunks,
    ) -> anyhow::Result<telemetry::Telemetry> {
        use super::crosslink_export::{proxl_xml, xiview_csv, ExportContext, SpectrumFile};
        use super::crosslink_report::{csm_report, members, FILE as REPORT, MEDIUM};
        info!("=== crosslink results ===");
        let CrosslinkChunks { mut matches, ms1_spectra, reporters, searched } = found;
        // telemetry counts the spectra searched for crosslinks
        self.spectra_searched.store(searched, Ordering::Relaxed);
        let silac = self.database.silac.is_some();
        let lfq = self.parameters.quant.lfq;
        let tmt = self.parameters.quant.tmt.clone();
        let library = self.parameters.spectral_library.clone();
        let annotator = self.parameters.write_annotator;

        let filenames = self
            .parameters
            .mzml_paths
            .iter()
            .map(|url| {
                sage_cloudpath::filename(url)
                    .unwrap_or_else(|| url.as_str())
                    .to_string()
            })
            .collect::<Vec<_>>();
        // How the viewers find each match's spectrum; spectrum positions are
        // read again only where ids alone are not enough (MGF, mzML without scan=)
        let mut files: Vec<SpectrumFile> = filenames.iter().map(|name| SpectrumFile::from_name(name)).collect();
        let reread: Vec<usize> = files
            .iter()
            .enumerate()
            .filter(|(file_id, file)| {
                let ids = matches.iter().filter(|m| m.file_id == *file_id && m.kind != LinkKind::Multi).map(|m| m.spec_id.as_str());
                file.needs_positions(ids)
            })
            .map(|(file_id, _)| file_id)
            .collect();

        let mut parts = Parts {
            n: 8 + library.is_some() as usize
                + !reread.is_empty() as usize
                + annotator as usize
                + lfq as usize
                + silac as usize,
            i: 0,
        };

        let part = parts.start("crosslink FDR", matches.len(), "matches", false);
        let (csms, links) = if settings.precursor_window {
            // each run's precursor mass error, learned from a first pass
            let (counts, windows) = assign_q_values_in_window(&mut matches);
            if windows.is_empty() {
                info!("precursor mass window: too few confident matches to learn it from; none applied");
            }
            let mut windows: Vec<_> = windows.into_iter().collect();
            windows.sort_by_key(|(file, _)| *file);
            for (file, (lo, hi)) in windows {
                let outside = matches.iter().filter(|m| m.file_id == file && (m.precursor_ppm < lo || m.precursor_ppm > hi)).count();
                let name = filenames.get(file).map(String::as_str).unwrap_or_default();
                info!("precursor mass window of {name}: {lo:+.1} to {hi:+.1} ppm ({outside} matches outside, left out)");
            }
            counts
        } else {
            assign_q_values(&mut matches)
        };
        matches.sort_by(|a, b| b.hyperscore.total_cmp(&a.hyperscore));
        log_fdr_groups(&matches);
        info!("discovered {} target crosslink-spectrum matches at 1% FDR", csms);
        info!("discovered {} unique target crosslinks at 1% FDR", links);
        part.done();

        if let Some(library) = &library {
            let part = parts.start("crosslink spectral library, xlibrary.tsv", matches.len(), "matches", false);
            self.write_crosslink_library(&matches, &settings, library)?;
            part.done();
        }

        if !reread.is_empty() {
            let part = parts.start("spectrum positions for the viewers", reread.len(), "files", true);
            for &file_id in &reread {
                let name = &filenames[file_id];
                part.detail(format!("now {name}"));
                match self
                    .parameters
                    .mzml_paths
                    .get(file_id)
                    .ok_or_else(|| anyhow::anyhow!("no path"))
                    .and_then(|url| Ok(sage_cloudpath::util::read_spectrum_positions(url)?))
                {
                    Ok(positions) => files[file_id].positions = Some(positions),
                    Err(e) => log::warn!(
                        "could not read spectrum positions of {name} ({e}): xiVIEW/ProXL may not find its spectra"
                    ),
                }
                part.add(1);
            }
            part.done();
        }

        // Protein sequences locate each peptide in its proteins
        let fasta = super::read_protein_databases(&mut self.parameters.database.clone())?;
        let wanted: std::collections::HashSet<&str> = matches
            .iter()
            .flat_map(|m| members(m).into_iter().map(|(ix, _)| ix))
            .flat_map(|ix| self.database[ix].proteins.iter().map(|p| p.as_ref()))
            .collect();
        let proteins: HashMap<String, String> = fasta
            .targets
            .iter()
            .filter(|(acc, _)| wanted.contains(acc.as_ref()))
            .map(|(acc, seq)| (acc.to_string(), seq.clone()))
            .collect();
        drop(fasta);
        // the FASTA file, or both crosslink databases
        let fasta_name = self
            .parameters
            .database
            .fasta_files()
            .into_iter()
            .map(|path| {
                sage_cloudpath::to_url(path)
                    .ok()
                    .and_then(|url| sage_cloudpath::filename(&url).map(String::from))
                    .unwrap_or_else(|| path.to_string())
            })
            .collect::<Vec<_>>()
            .join(" + ");
        let static_mods = self.parameters.database.modifications.static_mods.clone();
        let variable_mods = self.parameters.database.modifications.variable_mods.clone();
        let version = self.parameters.version.clone();
        let mut written = Vec::new();
        {
            let mut ctx = ExportContext {
                db: &self.database,
                settings: &settings,
                files: &files,
                proteins: &proteins,
                fasta_name: &fasta_name,
                fragment_tol: self.parameters.fragment_tol,
                ion_kinds: &self.database.ion_kinds,
                static_mods: &static_mods,
                variable_mods: &variable_mods,
                precursor_tol: self.parameters.precursor_tol,
                version: &version,
                progress: Default::default(),
                starts: Default::default(),
            };

            // Only the accepted results: target-target CSMs at 1% FDR whose
            // crosslinks are at 1% FDR too, so that unique crosslinks counted
            // from the file are at 1% (`CrosslinkMatch::accepted`)
            let accepted = matches.iter().filter(|m| m.accepted(FDR_CUTOFF)).cloned().collect::<Vec<_>>();
            let part = parts.start(
                &format!("crosslinks.sage.tsv and {FDR_FILTERED}"),
                matches.len() + accepted.len(),
                "matches",
                true,
            );
            ctx.progress = part.counter();
            let reporters_out = tmt.as_ref().map(|_| &reporters);
            written.push(self.write_crosslinks("crosslinks.sage.tsv", &matches, &filenames, reporters_out, &ctx)?);
            written.push(self.write_crosslinks(FDR_FILTERED, &accepted, &filenames, reporters_out, &ctx)?);
            let unique = accepted.iter().map(sage_core::crosslink::link_key).collect::<std::collections::HashSet<_>>().len();
            info!(
                "wrote {} ({} target crosslink-spectrum matches of {} unique crosslinks, both at 1% FDR)",
                FDR_FILTERED,
                accepted.len(),
                unique
            );
            part.done();

            // the crosslinked proteins at 1%, 5% and 10% FDR, as FASTA files
            let part = parts.start("crosslinked proteins as FASTA files", 0, "", false);
            written.extend(self.write_xlink_proteins(&matches)?);
            part.done();

            // the CSM report: every kind of match, multivalent species too
            let rows = matches.iter().filter(|m| m.class == PairClass::TT && m.spectrum_q <= MEDIUM).count();
            let part = parts.start(REPORT, rows, "matches", true);
            ctx.progress = part.counter();
            let (report, high, medium) = csm_report(&ctx, &matches, &filenames)?;
            let path = self.make_path(REPORT);
            sage_cloudpath::write_bytes_sync(&path, report)?;
            info!("wrote {REPORT}: {high} target crosslink-spectrum matches at 1% FDR (High), {medium} more at 5% (Medium)");
            written.push(path);
            part.done();

            // xiVIEW, ProXL and mzIdentML describe pairs and loop-links only
            let pairs: Vec<CrosslinkMatch> = matches.iter().filter(|m| m.kind != LinkKind::Multi).cloned().collect();
            let multi = matches.len() - pairs.len();
            if multi > 0 {
                info!("{multi} multivalent matches are in crosslinks.sage.tsv and {REPORT} only (not in the xiVIEW, ProXL and mzIdentML files)");
            }

            let part = parts.start("crosslinks_xiView.csv for xiVIEW", pairs.len(), "matches", true);
            ctx.progress = part.counter();
            let (csv, loops) = xiview_csv(&ctx, &pairs)?;
            let path = self.make_path("crosslinks_xiView.csv");
            sage_cloudpath::write_bytes_sync(&path, csv)?;
            if loops > 0 {
                info!("{} loop-link matches are not in crosslinks_xiView.csv: xiVIEW's CSV format cannot describe them", loops);
            }
            written.push(path);
            part.done();

            let proxl = pairs.iter().filter(|m| m.rank == 1 && m.class == PairClass::TT).count();
            let part = parts.start("crosslinks_proXL.xml for ProXL", proxl, "matches", true);
            ctx.progress = part.counter();
            match proxl_xml(&ctx, &pairs) {
                Some(xml) => {
                    let path = self.make_path("crosslinks_proXL.xml");
                    sage_cloudpath::write_bytes_sync(&path, xml.into_bytes())?;
                    written.push(path);
                }
                None => log::warn!("no target crosslink matches: crosslinks_proXL.xml was not written (ProXL needs at least one)"),
            }
            part.done();

            // mzIdentML 1.3.0 for xiVIEW (also a local installation) and PRIDE
            let part = parts.start("crosslinks.mzid for xiVIEW and PRIDE", pairs.len(), "matches", true);
            ctx.progress = part.counter();
            match super::crosslink_mzid::mzid_xml(&ctx, &pairs) {
                Some(xml) => {
                    let path = self.make_path("crosslinks.mzid");
                    sage_cloudpath::write_bytes_sync(&path, xml.into_bytes())?;
                    written.push(path);
                }
                None => log::warn!("no crosslink matches: crosslinks.mzid was not written (mzIdentML needs at least one)"),
            }
            part.done();
        }
        self.parameters.output_paths.extend(written);

        if annotator {
            let part = parts.start("crosslinks.annotator.csv for Annotator", 0, "", false);
            let path = self.write_annotator_crosslinks(&matches, &settings)?;
            self.parameters.output_paths.push(path);
            part.done();
        }
        if tmt.is_some() {
            let quantified = matches
                .iter()
                .filter(|m| reporters.contains_key(&(m.file_id, m.spec_id.clone())))
                .count();
            info!("isobaric reporter ions found for {} of {} crosslink-spectrum matches", quantified, matches.len());
        }

        if lfq {
            let part = parts.start("crosslink LFQ, crosslinks.lfq.tsv", 0, "", false);
            let areas = quantify_lfq(
                &matches,
                &settings,
                &self.parameters.quant.lfq_settings,
                &self.database,
                &ms1_spectra,
                self.parameters.mzml_paths.len(),
                self.parameters.precursor_charge,
            );
            let quantified = areas.iter().filter(|a| !a.decoy && a.peak.q_value <= 0.05).count();
            info!("discovered {} target crosslink MS1 peaks at 5% FDR", quantified);
            let path = self.write_crosslink_lfq(&matches, &areas, &filenames)?;
            self.parameters.output_paths.push(path);
            part.done();
        }
        if silac {
            let part = parts.start("SILAC ratios of the crosslinks", 0, "", false);
            if let Some(path) = self.write_silac_crosslinks(&matches, &settings, &ms1_spectra, &filenames)? {
                self.parameters.output_paths.push(path);
            }
            part.done();
        }

        let part = parts.start("quality control and results.json", 0, "", false);
        // polymer check and diagnostic ions (the linear search, if any, did them already)
        if self.parameters.qc.is_none() {
            self.finish_qc(None, &filenames)?;
        }
        // entrapment of the crosslinks, next to the linear search's
        let crosslink_entrapment = self.entrapment_crosslinks(&matches);
        if !crosslink_entrapment.is_empty() {
            let mut all = self.parameters.qc.as_ref().map(|q| q.entrapment.clone()).unwrap_or_default();
            all.extend(crosslink_entrapment);
            self.write_entrapment(&all)?;
            if let Some(qc) = self.parameters.qc.as_mut() {
                qc.entrapment = all;
            }
        }

        let path = self.make_path("results.json");
        let bytes = serde_json::to_vec_pretty(&self.parameters)?;
        sage_cloudpath::write_bytes_sync(&path, bytes)?;
        self.parameters.output_paths.push(path);
        part.done();
        // the settings, files and checks for people; results.json has them as JSON
        crate::summary::to_log(&self.parameters);

        self.log_time_summary();
        let run_time = (Instant::now() - self.start).as_secs();
        info!("finished in {}s", run_time);

        Ok(telemetry::Telemetry::new(
            self.parameters,
            self.database.peptides.len(),
            self.database.fragments.len(),
            self.spectra_searched.load(Ordering::Relaxed),
            false,
            run_time,
        ))
    }

    fn write_crosslinks(
        &self,
        name: &str,
        matches: &[CrosslinkMatch],
        filenames: &[String],
        reporters: Option<&Reporters>,
        ctx: &super::crosslink_export::ExportContext,
    ) -> anyhow::Result<sage_cloudpath::Url> {
        let path = self.make_path(name);
        let db = &self.database;

        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(vec![]);
        // Isobaric reporter columns follow the standard ones when TMT is set
        let reporter_headers = match (reporters, &self.parameters.quant.tmt) {
            (Some(_), Some(tmt)) => tmt.headers(),
            _ => Vec::new(),
        };
        let mut headers = vec![
            "psm_id",
            "filename",
            "scannr",
            "rank",
            "label",
            "class",
            "link_type",
            "peptide_a",
            "peptide_b",
            "link_a",
            "link_b",
            "proteins_a",
            "proteins_b",
            "intra_protein",
            "expmass",
            "calcmass",
            "charge",
            "isotope_error",
            "precursor_ppm",
            "hyperscore",
            "score_a",
            "score_b",
            "matched_a",
            "matched_b",
            "linked_a",
            "linked_b",
            "delta_next",
            "pairs_scored",
            "stub_peaks",
            "rt",
            "ion_mobility",
            "spectrum_q",
            "crosslink_q",
            "topology",
            "species",
            "topology_peaks",
        ]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
        // the link residues' positions in their proteins: a and b, and for
        // multivalent species one column per peptide (up to the valency + 1)
        let letters = abs_pos_letters(ctx.settings);
        headers.extend(letters.iter().map(|x| format!("abs_pos_{x}")));
        let silac = db.silac.is_some();
        if silac {
            headers.push("silac_channel".into());
        }
        let labelled = db.protein_labels.is_some();
        if labelled {
            headers.push("isotope_label_a".into());
            headers.push("isotope_label_b".into());
        }
        // two crosslink databases: which one each peptide is from
        let two_databases = db.two_databases();
        if two_databases {
            headers.push("database_a".into());
            headers.push("database_b".into());
        }
        if !reporter_headers.is_empty() {
            headers.push("ion_injection_time".into());
            headers.extend(reporter_headers.iter().cloned());
        }
        wtr.write_record(&headers)?;

        for m in matches {
            ctx.progress.add(1);
            let alpha = &db[m.alpha];
            let (peptide_b, proteins_b) = match m.beta {
                Some(ix) => (
                    db[ix].display(&db.modifications).to_string(),
                    db[ix].proteins(&db.decoy_tag, db.generate_decoys),
                ),
                None => (String::new(), String::new()),
            };
            // `link_a`/`link_b` are 1-based residue positions *within each peptide*.
            // For a loop-link both are positions on peptide_a.
            let mut record = vec![
                m.psm_id.to_string(),
                filenames
                    .get(m.file_id)
                    .cloned()
                    .unwrap_or_else(|| m.file_id.to_string()),
                m.spec_id.clone(),
                m.rank.to_string(),
                (if m.class == PairClass::TT { 1 } else { -1 }).to_string(),
                format!("{:?}", m.class),
                link_type(m.kind).to_string(),
                alpha.display(&db.modifications).to_string(),
                peptide_b,
                (m.alpha_site + 1).to_string(),
                (m.beta_site + 1).to_string(),
                alpha.proteins(&db.decoy_tag, db.generate_decoys),
                proteins_b,
                m.intra_protein.to_string(),
                m.expmass.to_string(),
                m.calcmass.to_string(),
                m.charge.to_string(),
                m.isotope_error.to_string(),
                m.precursor_ppm.to_string(),
                m.hyperscore.to_string(),
                m.alpha_score.to_string(),
                m.beta_score.to_string(),
                m.alpha_matched.to_string(),
                m.beta_matched.to_string(),
                m.alpha_linked.to_string(),
                m.beta_linked.to_string(),
                m.delta_next.to_string(),
                m.pairs_scored.to_string(),
                m.stub_peaks.to_string(),
                m.rt.to_string(),
                m.ims.to_string(),
                m.spectrum_q.to_string(),
                m.crosslink_q.to_string(),
                topology(m),
                m.species.as_deref().map(|s| species_text(db, s)).unwrap_or_default(),
                m.species.as_deref().map(|s| s.topology_peaks.to_string()).unwrap_or_default(),
            ];
            record.extend(absolute_positions(ctx, m, letters.len()));
            if silac {
                let peptides: Vec<_> = match &m.species {
                    Some(s) => s.peptides.clone(),
                    None => std::iter::once(m.alpha).chain(m.beta).collect(),
                };
                record.push(super::silac_tables::species_channel(db, &peptides));
            }
            if labelled {
                let name = |ix: Option<sage_core::database::PeptideIx>| {
                    ix.and_then(|ix| db.isotope_label(ix)).map_or("", |l| l.name()).to_string()
                };
                record.push(name(Some(m.alpha)));
                record.push(name(m.beta));
            }
            if two_databases {
                let group = |ix: Option<sage_core::database::PeptideIx>| {
                    ix.and_then(|ix| db.database_group(ix)).map(database_name).unwrap_or_default().to_string()
                };
                record.push(group(Some(m.alpha)));
                record.push(group(m.beta));
            }
            if let Some(reporters) = reporters.filter(|_| !reporter_headers.is_empty()) {
                // Spectra without reporter ions (e.g. no matching MS3) get empty cells
                match reporters.get(&(m.file_id, m.spec_id.clone())) {
                    Some((injection_time, peaks)) => {
                        record.push(injection_time.to_string());
                        record.extend(peaks.iter().map(|x| x.to_string()));
                    }
                    None => record.extend(std::iter::repeat_n(String::new(), reporter_headers.len() + 1)),
                }
            }
            wtr.write_record(&record)?;
        }

        wtr.flush()?;
        let bytes = wtr.into_inner()?;
        sage_cloudpath::write_bytes_sync(&path, bytes)?;
        Ok(path)
    }

    /// `crosslinks.lfq.tsv`: MS1 intensities of each quantified (target)
    /// crosslinked species in every run, named by its best match
    fn write_crosslink_lfq(
        &self,
        matches: &[CrosslinkMatch],
        areas: &[CrosslinkLfq],
        filenames: &[String],
    ) -> anyhow::Result<sage_cloudpath::Url> {
        let path = self.make_path("crosslinks.lfq.tsv");
        let db = &self.database;

        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(vec![]);
        let mut headers = [
            "peptide_a",
            "peptide_b",
            "link_type",
            "link_a",
            "link_b",
            "proteins_a",
            "proteins_b",
            "intra_protein",
            "calcmass",
            "charge",
            "crosslink_q",
            "q_value",
            "score",
            "spectral_angle",
        ]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
        headers.extend(filenames.iter().cloned());
        wtr.write_record(&headers)?;

        for quant in areas.iter().filter(|a| !a.decoy) {
            let m = &matches[quant.csm];
            let alpha = &db[m.alpha];
            let (peptide_b, proteins_b) = match m.beta {
                Some(ix) => (
                    db[ix].display(&db.modifications).to_string(),
                    db[ix].proteins(&db.decoy_tag, db.generate_decoys),
                ),
                None => (String::new(), String::new()),
            };
            let mut record = vec![
                alpha.display(&db.modifications).to_string(),
                peptide_b,
                link_type(m.kind).to_string(),
                (m.alpha_site + 1).to_string(),
                (m.beta_site + 1).to_string(),
                alpha.proteins(&db.decoy_tag, db.generate_decoys),
                proteins_b,
                m.intra_protein.to_string(),
                m.calcmass.to_string(),
                quant.charge.map_or(-1, i32::from).to_string(),
                m.crosslink_q.to_string(),
                quant.peak.q_value.to_string(),
                quant.peak.score.to_string(),
                quant.peak.spectral_angle.to_string(),
            ];
            record.extend(quant.areas.iter().map(|x| x.to_string()));
            wtr.write_record(&record)?;
        }

        wtr.flush()?;
        let bytes = wtr.into_inner()?;
        sage_cloudpath::write_bytes_sync(&path, bytes)?;
        Ok(path)
    }
}

#[cfg(test)]
mod test {
    use super::super::crosslink_export::test::{csm, find, fixture};
    use super::super::crosslink_export::{ExportContext, SpectrumFile};
    use super::*;
    use sage_core::crosslink::Species;
    use sage_core::fasta::Fasta;
    use sage_core::ion_series::Kind;
    use std::collections::HashSet;
    use std::sync::Arc;

    /// abs_pos_*: where the link residues are in their proteins
    #[test]
    fn absolute_positions_follow_the_link_columns() {
        let (db, statics, variables, proteins) = fixture();
        let settings = CrosslinkSettings { multivalent: Some(2), ..Default::default() };
        let files = vec![SpectrumFile::from_name("run01.d")];
        let ctx = ExportContext {
            db: &db,
            settings: &settings,
            files: &files,
            proteins: &proteins,
            fasta_name: "test.fasta",
            fragment_tol: Tolerance::Ppm(-10.0, 10.0),
            ion_kinds: &[Kind::B, Kind::Y],
            static_mods: &statics,
            variable_mods: &variables,
            precursor_tol: Tolerance::Ppm(-10.0, 10.0),
            version: "test",
            progress: Default::default(),
            starts: Default::default(),
        };
        // MAGLEK... is residue 1 on of ONE, CDKEFR 5 on of TWO, STKPLEK 16 on of TWO
        let a = find(&db, "MAGLEKVSTAAM[+15.9949]QR", false);
        let b = find(&db, "C[+57.021465]DKEFR", false);
        let c = find(&db, "STKPLEK", false);
        assert_eq!(abs_pos_letters(&settings), ['a', 'b', 'c']);
        assert_eq!(abs_pos_letters(&CrosslinkSettings::default()), ['a', 'b']);

        // a pair: link_a (K6 of the peptide) is K6 of ONE, link_b (K3) is K7 of TWO
        let pair = csm(a, Some(b), (5, 2), PairClass::TT, 1, 0.01);
        assert_eq!(absolute_positions(&ctx, &pair, 3), ["6", "7", ""]);
        // a loop-link: both residues in peptide_a's protein
        let looped = csm(c, None, (6, 2), PairClass::TT, 1, 0.01);
        assert_eq!(absolute_positions(&ctx, &looped, 2), ["22", "18"]);
        // a multivalent species: one column per peptide of its species column,
        // each with all its link residues
        let mut chain = csm(a, Some(b), (5, 2), PairClass::TT, 1, 0.01);
        chain.kind = LinkKind::Multi;
        chain.species = Some(Box::new(Species {
            peptides: vec![a, b, c],
            links: vec![(0, 5, 1, 2), (2, 2, 0, 7)],
            ..Default::default()
        }));
        assert_eq!(absolute_positions(&ctx, &chain, 3), ["6,8", "7", "18"]);
        // a decoy is placed where its target is
        let decoy = find(&db, "C[+57.021465]FEKDR", true);
        assert_eq!(ctx.absolute_positions(decoy, &[3]), "8");
    }

    /// The crosslink database of the linear search's proteins, and
    /// LinearProteins.fasta: each protein once
    #[test]
    fn linear_proteins_are_kept_once() {
        let text = ">sp|P1|ONE_HUMAN one\nMKTAYIAKQR\n>sp|P2|TWO_HUMAN two\nMSTNPKPQRK\n\
                    >sp|P1|ONE_HUMAN one, again\nMKTAYIAKQR\n>sp|P3|THREE_HUMAN three\nMADEEKLPPG\n";
        let keep: HashSet<Arc<str>> = ["sp|P1|ONE_HUMAN", "sp|P3|THREE_HUMAN"].into_iter().map(Arc::from).collect();
        let mut fasta = Fasta::parse(text.into(), "rev_", true);
        assert_eq!(keep_proteins(&mut fasta, &keep, "rev_", true), (2, 3));
        let accessions: Vec<&str> = fasta.targets.iter().map(|(a, _)| a.as_ref()).collect();
        assert_eq!(accessions, ["sp|P1|ONE_HUMAN", "sp|P3|THREE_HUMAN"]);

        // the FASTA file: each record once, its header as it was
        let wanted: std::collections::BTreeSet<&str> = keep.iter().map(|a| a.as_ref()).collect();
        let mut written = HashSet::new();
        let out = super::super::xlink_proteins::select_records(text, &wanted, &mut written);
        assert_eq!(out, ">sp|P1|ONE_HUMAN one\nMKTAYIAKQR\n>sp|P3|THREE_HUMAN three\nMADEEKLPPG\n");
        // ... also when the protein is in a second FASTA (two crosslink databases)
        let second = super::super::xlink_proteins::select_records(">sp|P3|THREE_HUMAN three\nMADEEKLPPG\n", &wanted, &mut written);
        assert_eq!(second, "");

        // decoys from the FASTA come with their targets
        let text = ">P1\nMKTAYIAKQR\n>rev_P1\nMQKAIYATKR\n>P2\nMSTNPKPQRK\n>rev_P2\nMKRQPKPNTS\n";
        let keep: HashSet<Arc<str>> = ["P1"].into_iter().map(Arc::from).collect();
        let mut fasta = Fasta::parse(text.into(), "rev_", false);
        assert_eq!(keep_proteins(&mut fasta, &keep, "rev_", false), (1, 2));
        let accessions: Vec<&str> = fasta.targets.iter().map(|(a, _)| a.as_ref()).collect();
        assert_eq!(accessions, ["P1", "rev_P1"]);
    }
}
