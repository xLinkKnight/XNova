//! `xlinkProtein_1pctFDR.fasta`, `xlinkProtein_5pctFDR.fasta` and
//! `xlinkProtein_10pctFDR.fasta`: the proteins of a crosslink search's
//! crosslinks at 1%, 5% and 10% FDR, copied from the search's FASTA file(s)
//! with their header lines as they were. A database of only the crosslinked
//! proteins, e.g. for a second search or for structural modelling, at
//! whichever FDR the results are used.
//!
//! A crosslink match counts at a level when it is target-target and its CSM
//! q-value (`spectrum_q`, as `crosslinks_FDRfiltered.sage.tsv`) or its
//! crosslink q-value (`crosslink_q`, as sage-gui's FDR setting for xiVIEW)
//! is at or below it, so that no protein of the results at that level is
//! missing, whichever of the two is used.

use super::Runner;
use log::{info, warn};
use sage_core::crosslink::{CrosslinkMatch, PairClass};
use sage_core::database::IndexedDatabase;
use std::collections::{BTreeSet, HashSet};

/// The FDR levels and their files, in the results folder
pub const LEVELS: [(f32, &str); 3] = [
    (0.01, "xlinkProtein_1pctFDR.fasta"),
    (0.05, "xlinkProtein_5pctFDR.fasta"),
    (0.10, "xlinkProtein_10pctFDR.fasta"),
];

/// A crosslink match within FDR `level`: target-target, with its CSM or its
/// crosslink q-value at or below `level`
pub fn within(m: &CrosslinkMatch, level: f32) -> bool {
    m.class == PairClass::TT && (m.spectrum_q <= level || m.crosslink_q <= level)
}

/// The proteins (accessions) of `matches`: both peptides of a pair, the
/// peptide of a loop-link and every peptide of a multivalent species, each
/// with every protein it is found in
pub fn crosslinked_proteins<'a, 'm>(db: &'a IndexedDatabase, matches: impl IntoIterator<Item = &'m CrosslinkMatch>) -> BTreeSet<&'a str> {
    matches
        .into_iter()
        .flat_map(|m| {
            let species = m.species.as_ref().map(|s| s.peptides.clone()).unwrap_or_default();
            std::iter::once(m.alpha).chain(m.beta).chain(species)
        })
        .flat_map(|ix| db[ix].proteins.iter().map(|p| p.as_ref()))
        .collect()
}

/// The records of FASTA `text` whose accession (the header's first word, as
/// Sage reads it) is in `wanted`, in the file's order, their header and
/// sequence lines as they were. `written` keeps a protein found in several
/// files to its first record.
pub fn select_records(text: &str, wanted: &BTreeSet<&str>, written: &mut HashSet<String>) -> String {
    let mut out = String::new();
    let mut keep = false;
    for line in text.lines() {
        let line = line.trim_end();
        if line.trim().is_empty() {
            continue;
        }
        if let Some(header) = line.trim_start().strip_prefix('>') {
            let accession = header.split_ascii_whitespace().next().unwrap_or("");
            keep = wanted.contains(accession) && written.insert(accession.to_string());
        }
        if keep {
            out.push_str(line.trim_start());
            out.push('\n');
        }
    }
    out
}

impl Runner {
    /// Write the crosslinked proteins at each FDR level (a level without
    /// crosslinks writes no file); returns the files written
    pub(super) fn write_xlink_proteins(&self, matches: &[CrosslinkMatch]) -> anyhow::Result<Vec<sage_cloudpath::Url>> {
        let (widest, _) = LEVELS[LEVELS.len() - 1];
        if !matches.iter().any(|m| within(m, widest)) {
            info!("no crosslinks at 10% FDR: no xlinkProtein FASTA file was written");
            return Ok(Vec::new());
        }
        // the FASTA file, or both crosslink databases, read once
        let mut texts = Vec::new();
        for path in self.parameters.database.fasta_files() {
            texts.push(sage_cloudpath::util::read_text(&sage_cloudpath::to_url(path)?)?);
        }
        let mut paths = Vec::new();
        let mut counts = Vec::new();
        for (level, name) in LEVELS {
            let wanted = crosslinked_proteins(&self.database, matches.iter().filter(|m| within(m, level)));
            let percent = (level * 100.0).round();
            if wanted.is_empty() {
                info!("no crosslinks at {percent}% FDR: {name} was not written");
                continue;
            }
            let mut written = HashSet::new();
            let fasta: String = texts.iter().map(|text| select_records(text, &wanted, &mut written)).collect();
            let missing: Vec<&str> = wanted.iter().copied().filter(|a| !written.contains(*a)).collect();
            if !missing.is_empty() {
                warn!("{} crosslinked protein(s) not found in the FASTA file(s), left out of {name}: {}", missing.len(), missing.join(", "));
            }
            let path = self.make_path(name);
            sage_cloudpath::write_bytes_sync(&path, fasta.into_bytes())?;
            counts.push(format!("{} at {percent}%", written.len()));
            paths.push(path);
        }
        info!("wrote the crosslinked proteins as FASTA files (xlinkProtein_*pctFDR.fasta): {} FDR", counts.join(", "));
        Ok(paths)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::runner::crosslink_export::test::{csm, find, fixture};
    use sage_core::crosslink::{LinkKind, Species};

    #[test]
    fn records_are_copied_whole() {
        let text = ">sp|P1|ONE_HUMAN First protein OS=Homo sapiens GN=ONE\r\n\
                    MKTAYIAKQR\r\n\
                    QISFVKSHFS\r\n\
                    \r\n\
                    >sp|P2|TWO_HUMAN Second protein\r\n\
                    MSTNPKPQRK\r\n\
                    >sp|P3|THREE_HUMAN Third\n\
                    MADEEKLPPG\n";
        let wanted: BTreeSet<&str> = ["sp|P1|ONE_HUMAN", "sp|P3|THREE_HUMAN"].into();
        let mut written = HashSet::new();
        let out = select_records(text, &wanted, &mut written);
        assert_eq!(
            out,
            ">sp|P1|ONE_HUMAN First protein OS=Homo sapiens GN=ONE\nMKTAYIAKQR\nQISFVKSHFS\n\
             >sp|P3|THREE_HUMAN Third\nMADEEKLPPG\n"
        );
        // reads back as the same proteins and sequences
        let back = sage_core::fasta::Fasta::parse(out, "rev_", true);
        let accessions: Vec<&str> = back.targets.iter().map(|(a, _)| a.as_ref()).collect();
        assert_eq!(accessions, ["sp|P1|ONE_HUMAN", "sp|P3|THREE_HUMAN"]);
        assert_eq!(back.targets[0].1, "MKTAYIAKQRQISFVKSHFS");
        // a protein in a second database too: once
        assert_eq!(select_records(">sp|P3|THREE_HUMAN again\nMADEEKLPPG\n", &wanted, &mut written), "");
    }

    #[test]
    fn proteins_of_pairs_loops_and_species() {
        let (db, ..) = fixture();
        let one = find(&db, "VSTAAMQR", false);
        let two = find(&db, "QLMNR", false);
        let pair = csm(one, Some(two), (0, 0), PairClass::TT, 1, 0.001);
        let looped = csm(one, None, (0, 3), PairClass::TT, 1, 0.001);
        assert_eq!(crosslinked_proteins(&db, [&pair]).into_iter().collect::<Vec<_>>(), ["sp|P11111|ONE_HUMAN", "sp|P22222|TWO_HUMAN"]);
        assert_eq!(crosslinked_proteins(&db, [&looped]).into_iter().collect::<Vec<_>>(), ["sp|P11111|ONE_HUMAN"]);
        // a multivalent species: every peptide's protein
        let mut multi = looped;
        multi.kind = LinkKind::Multi;
        multi.species = Some(Box::new(Species { peptides: vec![one, two], links: vec![(0, 0, 1, 0)], matched: vec![3, 3], ..Default::default() }));
        assert_eq!(crosslinked_proteins(&db, [&multi]).len(), 2);
        assert!(crosslinked_proteins(&db, []).is_empty());
    }

    /// The levels nest: a match within 1% is within 5% and 10%; either
    /// q-value counts; decoys never do
    #[test]
    fn matches_within_each_level() {
        let (db, ..) = fixture();
        let one = find(&db, "VSTAAMQR", false);
        let two = find(&db, "QLMNR", false);
        let at = |spectrum_q: f32, crosslink_q: f32, class: PairClass| {
            let mut m = csm(one, Some(two), (0, 0), class, 1, crosslink_q);
            m.spectrum_q = spectrum_q;
            LEVELS.iter().filter(|(level, _)| within(&m, *level)).count()
        };
        assert_eq!(at(0.004, 0.2, PairClass::TT), 3, "CSM q-value within 1%");
        assert_eq!(at(0.5, 0.004, PairClass::TT), 3, "crosslink q-value within 1%");
        assert_eq!(at(0.03, 0.2, PairClass::TT), 2, "5% and 10%");
        assert_eq!(at(0.2, 0.08, PairClass::TT), 1, "10% only");
        assert_eq!(at(0.2, 0.2, PairClass::TT), 0);
        assert_eq!(at(0.001, 0.001, PairClass::TD), 0, "decoys never");
        assert_eq!(at(0.001, 0.001, PairClass::DD), 0);
        let names: Vec<&str> = LEVELS.iter().map(|(_, n)| *n).collect();
        assert_eq!(names, ["xlinkProtein_1pctFDR.fasta", "xlinkProtein_5pctFDR.fasta", "xlinkProtein_10pctFDR.fasta"]);
    }
}
