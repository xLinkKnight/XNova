//! Digestion summary: how the sample was digested, judged from the peptides
//! identified at 1% FDR: how many have missed cleavages, and how many have a
//! terminus the enzyme did not make (semi- and non-enzymatic peptides).
//!
//! Each peptide is located in its proteins and its ends are checked against
//! the enzyme's cleavage sites in that protein (with its restriction rules and
//! every enzyme of a co-digestion). A peptide found in several proteins, or
//! several times, is classified at its most enzymatic occurrence. A decoy is
//! classified as its target, from which it was reversed. Decoys passing the
//! same filter estimate the false targets in each class.

use crate::enzyme::{EnzymeParameters, CLIPPED_AFTER_MET};
use serde::Serialize;
use std::collections::{HashMap, HashSet};

/// How one peptide was cut out of its protein
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Cut {
    /// Enzymatic sites inside the peptide
    pub missed_cleavages: usize,
    /// The N-terminus was made by the enzyme (or is the protein N-terminus)
    pub n_enzymatic: bool,
    /// The C-terminus was made by the enzyme (or is the protein C-terminus)
    pub c_enzymatic: bool,
}

impl Cut {
    /// Most enzymatic first: both ends, one end, none; then fewer missed cleavages
    fn rank(&self) -> (u8, usize) {
        (2 - self.n_enzymatic as u8 - self.c_enzymatic as u8, self.missed_cleavages)
    }
}

/// Classifies peptides by the cleavage sites of their proteins
pub struct Classifier<'a> {
    enzyme: &'a EnzymeParameters,
    proteins: &'a HashMap<String, String>,
    clip_n_term_met: bool,
    /// Cleavage positions of each protein, computed once
    cuts: std::sync::Mutex<HashMap<String, std::sync::Arc<HashSet<usize>>>>,
}

impl<'a> Classifier<'a> {
    /// `proteins`: accession → sequence
    pub fn new(enzyme: &'a EnzymeParameters, proteins: &'a HashMap<String, String>, clip_n_term_met: bool) -> Self {
        Self { enzyme, proteins, clip_n_term_met, cuts: Default::default() }
    }

    fn cuts(&self, accession: &str, sequence: &str) -> std::sync::Arc<HashSet<usize>> {
        if let Some(cuts) = self.cuts.lock().ok().and_then(|c| c.get(accession).cloned()) {
            return cuts;
        }
        let cuts: HashSet<usize> = self
            .enzyme
            .cleavage_sites(sequence)
            .into_iter()
            .flat_map(|site| [site.site.start, site.site.end])
            .collect();
        let cuts = std::sync::Arc::new(cuts);
        if let Ok(mut map) = self.cuts.lock() {
            map.insert(accession.to_string(), cuts.clone());
        }
        cuts
    }

    /// How `peptide` was cut from the first of `proteins` that holds it,
    /// at its most enzymatic occurrence. `None` when no protein holds it.
    /// A non-specific digest cuts everywhere: every peptide is enzymatic.
    pub fn classify(&self, peptide: &str, proteins: &[impl AsRef<str>]) -> Option<Cut> {
        let mut best: Option<Cut> = None;
        for accession in proteins {
            let Some(sequence) = self.proteins.get(accession.as_ref()) else { continue };
            let mut from = 0;
            while let Some(found) = sequence[from..].find(peptide) {
                let start = from + found;
                let end = start + peptide.len();
                let cut = if self.enzyme.enzyme.is_none() {
                    Cut { missed_cleavages: 0, n_enzymatic: true, c_enzymatic: true }
                } else {
                    let cuts = self.cuts(accession.as_ref(), sequence);
                    let residues = sequence.as_bytes();
                    let clipped = self.clip_n_term_met
                        && start == 1
                        && residues[0] == b'M'
                        && CLIPPED_AFTER_MET.contains(&residues[1]);
                    Cut {
                        missed_cleavages: cuts.iter().filter(|&&c| c > start && c < end).count(),
                        n_enzymatic: start == 0 || clipped || cuts.contains(&start),
                        c_enzymatic: end == sequence.len() || cuts.contains(&end),
                    }
                };
                if best.is_none_or(|b| cut.rank() < b.rank()) {
                    best = Some(cut);
                }
                from = start + 1;
            }
        }
        best
    }
}

/// The target a generated decoy was reversed from: the first and last
/// residues stay, the rest is reversed back
pub fn unreverse(decoy: &str) -> String {
    let mut residues: Vec<char> = decoy.chars().collect();
    let n = residues.len();
    if n > 2 {
        residues[1..n - 1].reverse();
    }
    residues.into_iter().collect()
}

/// Peptide counts of one file (or all files), each class estimated as
/// distinct targets minus distinct decoys, floored at zero
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct DigestionRow {
    pub file: String,
    pub target_peptides: usize,
    pub decoy_peptides: usize,
    /// Denominator of the percentages: targets minus decoys
    pub peptides: usize,
    pub missed_cleavages_0: usize,
    pub missed_cleavages_1: usize,
    pub missed_cleavages_2_plus: usize,
    /// N-terminus not made by the enzyme, C-terminus made by it
    pub semi_n: usize,
    /// C-terminus not made by the enzyme, N-terminus made by it
    pub semi_c: usize,
    /// Neither terminus made by the enzyme
    pub non_enzymatic: usize,
    pub missed_cleavage_pct: f64,
    pub semi_n_pct: f64,
    pub semi_c_pct: f64,
    pub non_enzymatic_pct: f64,
}

#[derive(Default, Clone, Copy)]
struct Counts {
    all: usize,
    mc0: usize,
    mc1: usize,
    mc2: usize,
    semi_n: usize,
    semi_c: usize,
    non: usize,
}

impl Counts {
    fn add(&mut self, cut: Option<Cut>) {
        self.all += 1;
        // a peptide no protein holds (e.g. from a peptide list) counts as enzymatic
        let cut = cut.unwrap_or(Cut { missed_cleavages: 0, n_enzymatic: true, c_enzymatic: true });
        match cut.missed_cleavages {
            0 => self.mc0 += 1,
            1 => self.mc1 += 1,
            _ => self.mc2 += 1,
        }
        match (cut.n_enzymatic, cut.c_enzymatic) {
            (false, true) => self.semi_n += 1,
            (true, false) => self.semi_c += 1,
            (false, false) => self.non += 1,
            (true, true) => {}
        }
    }
}

/// One identified peptide: its sequence (no modifications), proteins, and
/// whether it is a decoy
pub struct Identified<'a> {
    pub sequence: &'a str,
    pub proteins: &'a [std::sync::Arc<str>],
    pub decoy: bool,
}

/// The summary of one file's (or all files') distinct identified peptides
pub fn summarize<'a>(file: &str, peptides: impl Iterator<Item = Identified<'a>>, classifier: &Classifier) -> DigestionRow {
    let mut seen: HashSet<(String, bool)> = HashSet::new();
    let (mut targets, mut decoys) = (Counts::default(), Counts::default());
    for p in peptides {
        if !seen.insert((p.sequence.to_string(), p.decoy)) {
            continue;
        }
        if p.decoy {
            // generated decoys are reversed targets; decoys from the FASTA are found as written
            let cut = classifier.classify(p.sequence, p.proteins).or_else(|| {
                let target = unreverse(p.sequence);
                let accessions: Vec<String> = p.proteins.iter().map(|a| a.to_string()).collect();
                let targets: Vec<String> = accessions
                    .iter()
                    .flat_map(|a| [a.clone(), strip_decoy_prefix(a)])
                    .collect();
                classifier.classify(&target, &targets)
            });
            decoys.add(cut);
        } else {
            targets.add(classifier.classify(p.sequence, p.proteins));
        }
    }
    let est = |t: usize, d: usize| t.saturating_sub(d);
    let peptides = est(targets.all, decoys.all);
    let pct = |n: usize| if peptides > 0 { 100.0 * n as f64 / peptides as f64 } else { 0.0 };
    let (mc1, mc2) = (est(targets.mc1, decoys.mc1), est(targets.mc2, decoys.mc2));
    let (semi_n, semi_c, non) =
        (est(targets.semi_n, decoys.semi_n), est(targets.semi_c, decoys.semi_c), est(targets.non, decoys.non));
    DigestionRow {
        file: file.to_string(),
        target_peptides: targets.all,
        decoy_peptides: decoys.all,
        peptides,
        missed_cleavages_0: est(targets.mc0, decoys.mc0),
        missed_cleavages_1: mc1,
        missed_cleavages_2_plus: mc2,
        semi_n,
        semi_c,
        non_enzymatic: non,
        missed_cleavage_pct: pct(mc1 + mc2),
        semi_n_pct: pct(semi_n),
        semi_c_pct: pct(semi_c),
        non_enzymatic_pct: pct(non),
    }
}

/// `rev_P12345` → `P12345` (any prefix ending in `_` before the accession is
/// not stripped: only the usual decoy tags)
fn strip_decoy_prefix(accession: &str) -> String {
    for tag in ["rev_", "REV_", "decoy_", "DECOY_"] {
        if let Some(rest) = accession.strip_prefix(tag) {
            return rest.to_string();
        }
    }
    accession.to_string()
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::enzyme::{Enzyme, EnzymeParameters};
    use std::sync::Arc;

    fn trypsin() -> EnzymeParameters {
        EnzymeParameters { missed_cleavages: 2, min_len: 1, max_len: 50, enzyme: Enzyme::new("KR", "P", true, false) }
    }

    fn proteins() -> HashMap<String, String> {
        HashMap::from([("P1".to_string(), "MSDERAAAKPLLRGGGKTTTR".to_string())])
    }

    #[test]
    fn peptides_are_classified_by_their_protein() {
        let enzyme = trypsin();
        let proteins = proteins();
        let c = Classifier::new(&enzyme, &proteins, false);
        let p1 = ["P1"];
        let cut = |s: &str| c.classify(s, &p1).unwrap();
        // protein N-terminal, enzymatic
        assert_eq!(cut("MSDER"), Cut { missed_cleavages: 0, n_enzymatic: true, c_enzymatic: true });
        // K before P is not a site: AAAKPLLR has no missed cleavage
        assert_eq!(cut("AAAKPLLR"), Cut { missed_cleavages: 0, n_enzymatic: true, c_enzymatic: true });
        assert_eq!(cut("AAAKPLLRGGGK").missed_cleavages, 1);
        // ragged ends
        assert_eq!(cut("AAKPLLR"), Cut { missed_cleavages: 0, n_enzymatic: false, c_enzymatic: true });
        assert_eq!(cut("GGGKTT"), Cut { missed_cleavages: 1, n_enzymatic: true, c_enzymatic: false });
        assert_eq!(cut("GGKTT"), Cut { missed_cleavages: 1, n_enzymatic: false, c_enzymatic: false });
        // protein C-terminal
        assert!(cut("TTTR").c_enzymatic);
        assert_eq!(c.classify("WWWW", &p1), None);
        // without clipping, SDER starts at a non-enzymatic residue 2; with it, at the clipped N-terminus
        assert!(!cut("SDER").n_enzymatic);
        let clipping = Classifier::new(&enzyme, &proteins, true);
        assert!(clipping.classify("SDER", &p1).unwrap().n_enzymatic);
        // a non-specific digest: everything enzymatic
        let nonspecific = EnzymeParameters { missed_cleavages: 0, min_len: 1, max_len: 50, enzyme: None };
        let any = Classifier::new(&nonspecific, &proteins, false);
        assert_eq!(any.classify("GGKTT", &p1).unwrap(), Cut { missed_cleavages: 0, n_enzymatic: true, c_enzymatic: true });
    }

    #[test]
    fn the_summary_subtracts_decoys_in_each_class() {
        let enzyme = trypsin();
        let proteins = proteins();
        let c = Classifier::new(&enzyme, &proteins, false);
        let p1: Vec<Arc<str>> = vec![Arc::from("P1")];
        let rev: Vec<Arc<str>> = vec![Arc::from("rev_P1")];
        let ids = vec![
            Identified { sequence: "MSDER", proteins: &p1, decoy: false },
            Identified { sequence: "MSDER", proteins: &p1, decoy: false }, // counted once
            Identified { sequence: "AAAKPLLRGGGK", proteins: &p1, decoy: false },
            Identified { sequence: "AAKPLLR", proteins: &p1, decoy: false },
            Identified { sequence: "GGKTT", proteins: &p1, decoy: false },
            // the decoy of AAAKPLLRGGGK (one missed cleavage)
            Identified { sequence: "AGGGRLLPKAAK", proteins: &rev, decoy: true },
        ];
        assert_eq!(unreverse("AGGGRLLPKAAK"), "AAAKPLLRGGGK");
        let row = summarize("run1.mzML", ids.into_iter(), &c);
        assert_eq!((row.target_peptides, row.decoy_peptides, row.peptides), (4, 1, 3));
        assert_eq!((row.missed_cleavages_0, row.missed_cleavages_1, row.missed_cleavages_2_plus), (2, 1, 0));
        assert_eq!((row.semi_n, row.semi_c, row.non_enzymatic), (1, 0, 1));
        assert!((row.semi_n_pct - 100.0 / 3.0).abs() < 1e-9);
        assert!((row.missed_cleavage_pct - 100.0 / 3.0).abs() < 1e-9);
    }
}
