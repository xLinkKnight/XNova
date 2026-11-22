//! `crosslinks_CSMs.csv`: the accepted crosslink-spectrum matches as a table
//! to read, in the manner of MS Annika's and xiSEARCH's CSM tables.
//!
//! One row per target CSM at 5% FDR or better, best first, with its
//! confidence (High: q <= 1%, Medium: q <= 5%). Each peptide of a match has
//! columns of its own: A and B for a pair, A for a loop-link, and A, B, C, ...
//! for a multivalent species, as many as the largest species in the file
//! needs. Each crosslink is also given in its proteins (`ONE_HUMAN:K12-
//! TWO_HUMAN:K3`), from the protein sequences of the FASTA.

use super::crosslink_export::{mz, ExportContext};
use sage_core::crosslink::{rank_score, CrosslinkMatch, LinkKind, PairClass};
use sage_core::database::PeptideIx;

/// The report's file name, in the results folder
pub const FILE: &str = "crosslinks_CSMs.csv";
/// CSM q-values of the confidence classes: High at or below 1%, Medium at or below 5%
pub const HIGH: f32 = 0.01;
pub const MEDIUM: f32 = 0.05;

/// The peptides of a match, each with its 0-based link residues: one peptide
/// with both residues for a loop-link, two for a pair, and every peptide of a
/// multivalent species in its order (copies of one peptide separately)
pub fn members(m: &CrosslinkMatch) -> Vec<(PeptideIx, Vec<u16>)> {
    match (&m.species, m.beta) {
        (Some(s), _) => {
            let mut out: Vec<(PeptideIx, Vec<u16>)> = s.peptides.iter().map(|&p| (p, Vec::new())).collect();
            for &(a, sa, b, sb) in &s.links {
                if let Some(x) = out.get_mut(a as usize) {
                    x.1.push(sa);
                }
                if let Some(x) = out.get_mut(b as usize) {
                    x.1.push(sb);
                }
            }
            for (_, sites) in &mut out {
                sites.sort_unstable();
                sites.dedup();
            }
            out
        }
        (None, Some(b)) => vec![(m.alpha, vec![m.alpha_site]), (b, vec![m.beta_site])],
        (None, None) => {
            let mut sites = vec![m.alpha_site, m.beta_site];
            sites.sort_unstable();
            vec![(m.alpha, sites)]
        }
    }
}

/// The links of a match as `(peptide, residue, peptide, residue)`, indices
/// into [`members`] and 0-based residues
fn links(m: &CrosslinkMatch) -> Vec<(usize, u16, usize, u16)> {
    match (&m.species, m.beta) {
        (Some(s), _) => s.links.iter().map(|&(a, sa, b, sb)| (a as usize, sa, b as usize, sb)).collect(),
        (None, Some(_)) => vec![(0, m.alpha_site, 1, m.beta_site)],
        (None, None) => vec![(0, m.alpha_site.min(m.beta_site), 0, m.alpha_site.max(m.beta_site))],
    }
}

/// Each peptide's own score and matched fragment ions
fn peptide_scores(m: &CrosslinkMatch) -> Vec<(Option<f64>, u16)> {
    match (&m.species, m.beta) {
        (Some(s), _) => (0..s.peptides.len()).map(|i| (s.scores.get(i).copied(), s.matched.get(i).copied().unwrap_or(0))).collect(),
        (None, Some(_)) => vec![(Some(m.alpha_score), m.alpha_matched), (Some(m.beta_score), m.beta_matched)],
        (None, None) => vec![(Some(m.hyperscore), m.alpha_matched)],
    }
}

/// `Crosslink type` column
fn crosslink_type(m: &CrosslinkMatch) -> &'static str {
    match m.kind {
        LinkKind::Loop => "Loop-link",
        LinkKind::Multi => "Multivalent",
        LinkKind::Inter if m.beta == Some(m.alpha) => "Homomeric",
        LinkKind::Inter if m.intra_protein => "Intra-protein",
        LinkKind::Inter => "Inter-protein",
    }
}

/// `Topology` column, as `crosslinks.sage.tsv`'s
fn topology(m: &CrosslinkMatch) -> String {
    match (m.kind, &m.species, m.beta) {
        (LinkKind::Multi, Some(s), _) => s.topology(),
        (LinkKind::Loop, ..) => "loop".into(),
        (_, _, Some(b)) if b == m.alpha => "a1a".into(),
        _ => "a1b".into(),
    }
}

/// The reagent's name, or the custom linker's mass
fn crosslinker(ctx: &ExportContext) -> String {
    match ctx.settings.reagent {
        Some(reagent) => reagent.name().into(),
        None => format!("custom ({:.4} Da)", ctx.settings.linker_mass()),
    }
}

/// A protein's short name: the entry name of a UniProt header
/// (`sp|P04050|RPB1_YEAST` → `RPB1_YEAST`), else the accession as it is
pub fn short_name(accession: &str) -> &str {
    let parts: Vec<&str> = accession.split('|').collect();
    match parts.as_slice() {
        [db, _, entry, ..] if matches!(*db, "sp" | "tr") && !entry.is_empty() => entry,
        _ => accession,
    }
}

fn fixed(x: f64, digits: usize) -> String {
    if x.is_finite() {
        format!("{x:.digits$}")
    } else {
        String::new()
    }
}

/// The report, and how many CSMs it holds at each confidence: `(bytes, high, medium)`
pub fn csm_report(ctx: &ExportContext, matches: &[CrosslinkMatch], filenames: &[String]) -> anyhow::Result<(Vec<u8>, usize, usize)> {
    let mut rows: Vec<&CrosslinkMatch> = matches.iter().filter(|m| m.class == PairClass::TT && m.spectrum_q <= MEDIUM).collect();
    rows.sort_by(|a, b| {
        a.spectrum_q
            .total_cmp(&b.spectrum_q)
            .then(rank_score(b).total_cmp(&rank_score(a)))
            .then(b.hyperscore.total_cmp(&a.hyperscore))
            .then(a.psm_id.cmp(&b.psm_id))
    });
    let blocks = rows.iter().map(|m| members(m).len()).max().unwrap_or(2).max(2);
    let db = ctx.db;
    let linker = crosslinker(ctx);

    let mut wtr = csv::WriterBuilder::new().from_writer(vec![]);
    let mut header: Vec<String> = [
        "CSM ID",
        "Confidence",
        "CSM q-value",
        "Crosslink q-value",
        "Crosslink type",
        "Topology",
        "Peptides",
        "Sequence",
        "Crosslinks",
        "Score",
        "Hyperscore",
        "Stub peaks",
        "Topology peaks",
        "Crosslinker",
        "Spectrum file",
        "Scan",
        "RT [min]",
        "Ion mobility [1/K0]",
        "Charge",
        "Precursor m/z",
        "Experimental mass [Da]",
        "Calculated mass [Da]",
        "Mass error [ppm]",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    for i in 0..blocks {
        let x = (b'A' + i as u8) as char;
        for name in ["Peptide", "Modified peptide", "Link sites", "Proteins", "Peptide start", "Protein link sites", "Score", "Matched ions"] {
            header.push(format!("{name} {x}"));
        }
    }
    wtr.write_record(&header)?;

    let (mut high, mut medium) = (0, 0);
    for m in rows {
        ctx.progress.add(1);
        if m.spectrum_q <= HIGH {
            high += 1;
        } else {
            medium += 1;
        }
        let peptides = members(m);
        let scores = peptide_scores(m);
        // each peptide's proteins and 1-based start in each
        let starts: Vec<_> = peptides.iter().map(|(ix, _)| ctx.starts_of(*ix)).collect();
        let residue = |p: usize, site: u16| db[peptides[p].0].sequence.get(site as usize).map_or('?', |&c| c as char);
        // a link residue in the first protein of its peptide
        let in_protein = |p: usize, site: u16| match starts[p].first() {
            Some((name, Some(start))) => format!("{}:{}{}", short_name(name), residue(p, site), start + site as usize),
            Some((name, None)) => format!("{}:{}?", short_name(name), residue(p, site)),
            None => format!("?:{}?", residue(p, site)),
        };
        let crosslinks = links(m)
            .iter()
            .map(|&(a, sa, b, sb)| format!("{}-{}", in_protein(a, sa), in_protein(b, sb)))
            .collect::<Vec<_>>()
            .join("; ");
        let sequence = peptides
            .iter()
            .map(|(ix, _)| String::from_utf8_lossy(&db[*ix].sequence).to_string())
            .collect::<Vec<_>>()
            .join("-");
        let file = filenames.get(m.file_id).cloned().unwrap_or_else(|| m.file_id.to_string());
        let mut record = vec![
            m.psm_id.to_string(),
            (if m.spectrum_q <= HIGH { "High" } else { "Medium" }).to_string(),
            m.spectrum_q.to_string(),
            m.crosslink_q.to_string(),
            crosslink_type(m).to_string(),
            topology(m),
            peptides.len().to_string(),
            sequence,
            crosslinks,
            fixed(rank_score(m), 2),
            fixed(m.hyperscore, 2),
            if ctx.settings.is_cleavable() { m.stub_peaks.to_string() } else { String::new() },
            m.species.as_deref().map(|s| s.topology_peaks.to_string()).unwrap_or_default(),
            linker.clone(),
            file,
            m.spec_id.clone(),
            fixed(m.rt as f64, 3),
            if m.ims > 0.0 { fixed(m.ims as f64, 4) } else { String::new() },
            m.charge.to_string(),
            fixed(mz(m.expmass, m.charge) as f64, 5),
            fixed(m.expmass as f64, 4),
            fixed(m.calcmass as f64, 4),
            fixed(m.precursor_ppm as f64, 2),
        ];
        for i in 0..blocks {
            let Some((ix, sites)) = peptides.get(i) else {
                record.extend(std::iter::repeat_n(String::new(), 8));
                continue;
            };
            let p = &db[*ix];
            let site_text = |offset: Option<usize>| {
                sites
                    .iter()
                    .map(|&s| format!("{}{}", residue(i, s), offset.map_or("?".to_string(), |o| (o + s as usize).to_string())))
                    .collect::<Vec<_>>()
                    .join(",")
            };
            let (score, matched) = scores.get(i).copied().unwrap_or((None, 0));
            record.push(String::from_utf8_lossy(&p.sequence).to_string());
            record.push(p.display(&db.modifications).to_string());
            record.push(site_text(Some(1)));
            record.push(starts[i].iter().map(|(name, _)| name.as_str()).collect::<Vec<_>>().join(";"));
            record.push(starts[i].iter().map(|(_, s)| s.map_or("?".to_string(), |s| s.to_string())).collect::<Vec<_>>().join(";"));
            record.push(starts[i].iter().map(|(_, s)| site_text(*s)).collect::<Vec<_>>().join(";"));
            record.push(score.map(|s| fixed(s, 2)).unwrap_or_default());
            record.push(matched.to_string());
        }
        wtr.write_record(&record)?;
    }
    wtr.flush()?;
    Ok((wtr.into_inner()?, high, medium))
}

#[cfg(test)]
mod test {
    use super::super::crosslink_export::test::{csm, find, fixture};
    use super::super::crosslink_export::SpectrumFile;
    use super::*;
    use sage_core::crosslink::{CrosslinkSettings, Reagent, Species};
    use sage_core::ion_series::Kind;
    use sage_core::mass::Tolerance;

    fn read(bytes: Vec<u8>) -> Vec<std::collections::HashMap<String, String>> {
        let mut rdr = csv::Reader::from_reader(bytes.as_slice());
        let header: Vec<String> = rdr.headers().unwrap().iter().map(String::from).collect();
        rdr.records()
            .map(|r| header.iter().cloned().zip(r.unwrap().iter().map(String::from)).collect())
            .collect()
    }

    #[test]
    fn report_of_pairs_loops_and_species() {
        let (db, statics, variables, proteins) = fixture();
        let settings = CrosslinkSettings { reagent: Some(Reagent::Dsso), ..Default::default() };
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
        // MAGLEK... starts at residue 1 of P11111 (ONE), CDKEFR at 5 of P22222 (TWO), STKPLEK at 16 of TWO
        let a = find(&db, "MAGLEKVSTAAM[+15.9949]QR", false);
        let b = find(&db, "C[+57.021465]DKEFR", false);
        let c = find(&db, "STKPLEK", false);
        let decoy = find(&db, "C[+57.021465]FEKDR", true);
        let mut pair = csm(a, Some(b), (5, 2), PairClass::TT, 1, 0.002);
        pair.spectrum_q = 0.004;
        let mut looped = csm(c, None, (6, 2), PairClass::TT, 1, 0.003);
        looped.spectrum_q = 0.002;
        looped.hyperscore = 33.0;
        looped.ims = 1.0123;
        // a1b1c: MAGLEK...'s K6 to CDKEFR's K3, and STKPLEK's K3 to MAGLEK...'s S8
        let mut chain = csm(a, Some(b), (5, 2), PairClass::TT, 1, 0.03);
        chain.kind = LinkKind::Multi;
        chain.spectrum_q = 0.03;
        chain.species = Some(Box::new(Species {
            peptides: vec![a, b, c],
            links: vec![(0, 5, 1, 2), (2, 2, 0, 7)],
            matched: vec![9, 6, 4],
            scores: vec![30.0, 20.0, 12.5],
            ..Default::default()
        }));
        let mut weak = csm(a, Some(b), (5, 2), PairClass::TT, 1, 0.2);
        weak.spectrum_q = 0.2;
        let mut false_pair = csm(a, Some(decoy), (5, 3), PairClass::TD, 1, 0.001);
        false_pair.spectrum_q = 0.001;
        let matches = vec![weak, chain, false_pair, pair, looped];
        let (bytes, high, medium) = csm_report(&ctx, &matches, &["run01.d".to_string()]).unwrap();
        let text = String::from_utf8(bytes.clone()).unwrap();
        let rows = read(bytes);
        assert_eq!((rows.len(), high, medium), (3, 2, 1), "targets at 5%, no decoys: {text}");
        // three peptide blocks: the species has three
        assert!(text.lines().next().unwrap().ends_with("Peptide C,Modified peptide C,Link sites C,Proteins C,Peptide start C,Protein link sites C,Score C,Matched ions C"));

        // best first: the loop-link (q 0.002), the pair (0.004), the species (0.03)
        let looped = &rows[0];
        assert_eq!(looped["Crosslink type"], "Loop-link");
        assert_eq!(looped["Confidence"], "High");
        assert_eq!(looped["Peptides"], "1");
        assert_eq!(looped["Link sites A"], "K3,K7");
        assert_eq!(looped["Peptide start A"], "16");
        assert_eq!(looped["Protein link sites A"], "K18,K22");
        assert_eq!(looped["Crosslinks"], "TWO_HUMAN:K18-TWO_HUMAN:K22");
        assert_eq!(looped["Score"], "33.00");
        assert_eq!(looped["Ion mobility [1/K0]"], "1.0123");
        assert_eq!(looped["Peptide B"], "");

        let pair = &rows[1];
        assert_eq!(pair["Crosslink type"], "Inter-protein");
        assert_eq!(pair["Topology"], "a1b");
        assert_eq!(pair["Sequence"], "MAGLEKVSTAAMQR-CDKEFR");
        assert_eq!(pair["Modified peptide A"], "MAGLEKVSTAAM[+15.9949]QR");
        assert_eq!(pair["Crosslinks"], "ONE_HUMAN:K6-TWO_HUMAN:K7");
        assert_eq!(pair["Proteins A"], "sp|P11111|ONE_HUMAN");
        assert_eq!((pair["Link sites A"].as_str(), pair["Protein link sites B"].as_str()), ("K6", "K7"));
        assert_eq!(pair["Score"], "20.00", "the weaker peptide's");
        assert_eq!((pair["Score A"].as_str(), pair["Score B"].as_str(), pair["Matched ions B"].as_str()), ("30.00", "20.00", "6"));
        assert_eq!(pair["Crosslinker"], "DSSO");
        assert_eq!(pair["Spectrum file"], "run01.d");
        assert_eq!(pair["Peptide C"], "");

        let chain = &rows[2];
        assert_eq!(chain["Confidence"], "Medium");
        assert_eq!(chain["Crosslink type"], "Multivalent");
        assert_eq!(chain["Topology"], "a1b1c");
        assert_eq!(chain["Peptides"], "3");
        assert_eq!(chain["Sequence"], "MAGLEKVSTAAMQR-CDKEFR-STKPLEK");
        assert_eq!(chain["Crosslinks"], "ONE_HUMAN:K6-TWO_HUMAN:K7; TWO_HUMAN:K18-ONE_HUMAN:S8");
        assert_eq!((chain["Peptide C"].as_str(), chain["Link sites C"].as_str(), chain["Protein link sites C"].as_str()), ("STKPLEK", "K3", "K18"));
        assert_eq!((chain["Score C"].as_str(), chain["Matched ions C"].as_str()), ("12.50", "4"));
        assert_eq!((chain["Link sites A"].as_str(), chain["Protein link sites A"].as_str()), ("K6,S8", "K6,S8"), "a peptide with two links");

        // nothing at 5%: a header with two peptide blocks
        let (bytes, high, medium) = csm_report(&ctx, &[], &[]).unwrap();
        assert_eq!((high, medium), (0, 0));
        let text = String::from_utf8(bytes).unwrap();
        assert_eq!(text.lines().count(), 1);
        assert!(text.contains("Matched ions B") && !text.contains("Peptide C"));
        assert_eq!(short_name("sp|P04050|RPB1_YEAST"), "RPB1_YEAST");
        assert_eq!(short_name("custom_protein"), "custom_protein");
    }
}
