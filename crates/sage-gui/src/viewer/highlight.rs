//! The spectrum and the peptide map, linked: the residues of an ion picked
//! in the spectrum are framed in the map, and the peaks of a glyph picked in
//! the map (a hook, a stub triangle) are picked out in the spectrum.

use super::model::{Annotation, Kind, Series, StubIonKind};
use std::collections::BTreeSet;
use std::ops::Range;

/// A glyph of the peptide map that stands for ions
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Glyph {
    /// A hook at `residue`: the backbone ions ending after it (N-terminal:
    /// a, b, c) or starting at it (C-terminal: x, y, z)
    Hook { peptide: usize, n_term: bool, residue: usize },
    /// A nested hook: those ions with linker stub `stub` in place of the
    /// partner (MS-cleavable linkers)
    StubHook { peptide: usize, n_term: bool, residue: usize, stub: usize },
    /// A triangle by a bridge: the whole peptide with stub `stub` (its
    /// doublet peaks); in a multivalent species, `link` is the linker that
    /// broke, and the triangle also stands for the parts freed when only it broke
    Stub { peptide: usize, stub: usize, link: Option<usize> },
    /// Multivalent species: a peptide with two or more links freed whole,
    /// with a stub on each (`stubs`, a bit per stub: its triplet peaks)
    Stubs { peptide: usize, stubs: u8 },
}

/// A bit per stub
pub fn stub_mask(stubs: &[usize]) -> u8 {
    stubs.iter().fold(0, |mask, &j| mask | (1u8 << j.min(7)))
}

/// What the pointer is on, or what was clicked
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Focus {
    /// A peak of the spectrum (index into [`Annotation::peaks`])
    Peak(usize),
    Glyph(Glyph),
}

/// One interpretation of a peak: the glyph that stands for it in the map,
/// and the residues the ion holds
#[derive(Clone, Debug, PartialEq)]
pub struct Ion {
    pub peak: usize,
    pub text: String,
    pub theoretical_mz: f64,
    pub glyph: Option<Glyph>,
    /// (peptide, residues), 0-based
    pub residues: Vec<(usize, Range<usize>)>,
}

/// The residues a fragment of `peptide` made of `range` holds: its own, and
/// when it holds a link residue (the linker intact), every peptide linked
/// through it, whole (as the engine computes the fragment's mass)
fn held(a: &Annotation, peptide: usize, range: Range<usize>) -> Vec<(usize, Range<usize>)> {
    let len = |k: usize| a.coverage.get(k).map_or(0, |c| c.residues.len());
    let mut out = vec![(peptide, range.clone())];
    let Some(link) = a.link.as_ref().filter(|l| l.kind != Kind::Loop) else { return out };
    let mut seen = vec![peptide];
    let mut next = Vec::new();
    for &(pa, ra, pb, rb) in &link.links {
        for (own, residue, other) in [(pa, ra, pb), (pb, rb, pa)] {
            if own == peptide && range.contains(&residue) && !seen.contains(&other) {
                seen.push(other);
                next.push(other);
            }
        }
    }
    // and the peptides linked to those
    while let Some(k) = next.pop() {
        out.push((k, 0..len(k)));
        for &(pa, _, pb, _) in &link.links {
            for (from, to) in [(pa, pb), (pb, pa)] {
                if from == k && !seen.contains(&to) {
                    seen.push(to);
                    next.push(to);
                }
            }
        }
    }
    out
}

/// Every interpretation of every peak of `a`, with its glyph and residues
pub fn ions(a: &Annotation) -> Vec<Ion> {
    let len = |k: usize| a.coverage.get(k).map_or(0, |c| c.residues.len());
    let mut out = Vec::new();
    for (p, peak) in a.peaks.iter().enumerate() {
        for l in &peak.labels {
            let ion = |glyph, residues| Ion { peak: p, text: l.text.clone(), theoretical_mz: l.theoretical_mz, glyph, residues };
            let n = len(l.peptide);
            match (l.series, l.backbone, l.internal) {
                // stub ions, below (they say where they lie)
                (Series::Stub(_), ..) => {}
                (_, Some((n_term, r)), _) => {
                    let range = if n_term { 0..(r + 1).min(n) } else { r.min(n)..n };
                    let glyph = Glyph::Hook { peptide: l.peptide, n_term, residue: r };
                    out.push(ion(Some(glyph), held(a, l.peptide, range)));
                }
                // internal ions: no glyph of their own, but residues to frame
                (_, _, Some((first, last))) => out.push(ion(None, held(a, l.peptide, first.min(n)..(last + 1).min(n)))),
                (Series::Precursor, ..) => out.push(ion(None, (0..a.coverage.len()).map(|k| (k, 0..len(k))).collect())),
                // immonium and other ions: nothing in the map
                _ => out.push(ion(None, Vec::new())),
            }
        }
    }
    for s in &a.stub_ions {
        let n = len(s.peptide);
        let (glyphs, residues): (Vec<Glyph>, Vec<(usize, Range<usize>)>) = match s.kind {
            StubIonKind::NTerm(r) | StubIonKind::CTerm(r) => {
                let n_term = matches!(s.kind, StubIonKind::NTerm(_));
                let range = if n_term { 0..(r + 1).min(n) } else { r.min(n)..n };
                let glyphs = s.stubs.iter().map(|&stub| Glyph::StubHook { peptide: s.peptide, n_term, residue: r, stub }).collect();
                (glyphs, vec![(s.peptide, range)])
            }
            StubIonKind::Peptide => {
                let glyph = match a.link.as_ref().filter(|l| l.kind == Kind::Multi) {
                    // a species: a peptide with one link, by the bridge of that link;
                    // one with two or more, after its row (as the map draws them)
                    Some(link) if link.ends_of(s.peptide).len() == 1 => {
                        let e = link.links.iter().position(|&(pa, _, pb, _)| pa == s.peptide || pb == s.peptide);
                        Glyph::Stub { peptide: s.peptide, stub: s.stub, link: e }
                    }
                    Some(_) => Glyph::Stubs { peptide: s.peptide, stubs: stub_mask(&s.stubs) },
                    None => Glyph::Stub { peptide: s.peptide, stub: s.stub, link: None },
                };
                (vec![glyph], vec![(s.peptide, 0..n)])
            }
            StubIonKind::Part { link, peptides } => (
                vec![Glyph::Stub { peptide: s.peptide, stub: s.stub, link: Some(link as usize) }],
                (0..8).filter(|k| peptides & (1 << k) != 0).map(|k| (k, 0..len(k))).collect(),
            ),
        };
        for glyph in glyphs {
            out.push(Ion { peak: s.peak, text: s.text.clone(), theoretical_mz: s.theoretical_mz, glyph: Some(glyph), residues: residues.clone() });
        }
    }
    out
}

/// What to pick out: the peaks in the spectrum, the residues and glyphs in the map
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Highlight {
    pub peaks: BTreeSet<usize>,
    pub residues: Vec<(usize, Range<usize>)>,
    pub glyphs: BTreeSet<Glyph>,
    /// The ions picked (indices into [`ions`])
    pub ions: Vec<usize>,
}

impl Highlight {
    /// A peak: its interpretations' residues and glyphs; a glyph: the peaks
    /// of its ions, and their residues
    pub fn of(ions: &[Ion], focus: Focus) -> Highlight {
        let mut h = Highlight::default();
        for (i, ion) in ions.iter().enumerate() {
            let picked = match focus {
                Focus::Peak(p) => ion.peak == p,
                Focus::Glyph(g) => ion.glyph == Some(g),
            };
            if picked {
                h.ions.push(i);
                h.peaks.insert(ion.peak);
                h.glyphs.extend(ion.glyph);
                h.residues.extend(ion.residues.iter().cloned());
            }
        }
        match focus {
            Focus::Peak(p) => h.peaks.insert(p),
            Focus::Glyph(g) => h.glyphs.insert(g),
        };
        h
    }

    pub fn is_empty(&self) -> bool {
        self.peaks.is_empty() && self.glyphs.is_empty()
    }

    /// The residue ranges of `peptide`, joined where they overlap or touch
    pub fn ranges(&self, peptide: usize) -> Vec<Range<usize>> {
        let mut ranges: Vec<Range<usize>> =
            self.residues.iter().filter(|(k, r)| *k == peptide && !r.is_empty()).map(|(_, r)| r.clone()).collect();
        ranges.sort_by_key(|r| r.start);
        let mut out: Vec<Range<usize>> = Vec::new();
        for r in ranges {
            match out.last_mut() {
                Some(last) if r.start <= last.end => last.end = last.end.max(r.end),
                _ => out.push(r),
            }
        }
        out
    }

    pub fn has_residue(&self, peptide: usize, residue: usize) -> bool {
        self.residues.iter().any(|(k, r)| *k == peptide && r.contains(&residue))
    }

    /// One line on what is picked out: the ions with their m/z, and the residues
    pub fn describe(&self, a: &Annotation, ions: &[Ion]) -> String {
        if self.ions.is_empty() {
            return match self.peaks.first() {
                Some(&p) => format!("m/z {:.4}: not annotated", a.peaks[p].mz),
                None => String::new(),
            };
        }
        let mut names: Vec<String> = Vec::new();
        for &i in &self.ions {
            let ion = &ions[i];
            let observed = a.peaks[ion.peak].mz;
            let ppm = (observed - ion.theoretical_mz) / ion.theoretical_mz * 1e6;
            let name = format!("{} {observed:.4} ({ppm:+.1} ppm)", ion.text);
            if !names.contains(&name) {
                names.push(name);
            }
        }
        let mut parts: Vec<String> = Vec::new();
        for k in 0..a.coverage.len() {
            for r in self.ranges(k) {
                let residues = &a.coverage[k].residues;
                let what = if r.len() == residues.len() {
                    "whole".to_string()
                } else {
                    format!("{}{}–{}{}", residues[r.start], r.start + 1, residues[r.end - 1], r.end)
                };
                let peptide = if a.coverage.len() > 1 { format!("peptide {} ", k + 1) } else { String::new() };
                parts.push(format!("{peptide}{what}"));
            }
        }
        let peaks = if self.peaks.len() > 1 { format!("{} peaks: ", self.peaks.len()) } else { String::new() };
        if parts.is_empty() {
            format!("{peaks}{}", names.join(", "))
        } else {
            format!("{peaks}{} · {}", names.join(", "), parts.join(", "))
        }
    }
}

#[cfg(test)]
mod test {
    use super::super::model;
    use super::*;

    /// The ions of a cleaved DSSO crosslink, both ways: a peak frames its
    /// residues, a glyph picks out its peaks
    #[test]
    fn peaks_and_glyphs_point_at_each_other() {
        let (_, a) = model::test::cleaved_dsso();
        let ions = ions(&a);
        // y5+S of TDYHFPNEWKLIR: from W9 on, the S stub in place of the partner
        let y5 = a.peaks.iter().position(|p| p.labels[0].text == "y5+S").unwrap();
        let h = Highlight::of(&ions, Focus::Peak(y5));
        assert_eq!(h.residues, vec![(1, 8..13)]);
        assert_eq!(h.glyphs.iter().copied().collect::<Vec<_>>(), vec![Glyph::StubHook { peptide: 1, n_term: false, residue: 8, stub: 2 }]);
        assert!(h.has_residue(1, 9) && !h.has_residue(1, 7) && !h.has_residue(0, 3));
        assert!(h.describe(&a, &ions).starts_with("y5+S "), "{}", h.describe(&a, &ions));
        assert!(h.describe(&a, &ions).ends_with("peptide 2 W9–R13"), "{}", h.describe(&a, &ions));

        // the reverse: the nested hook picks out that peak alone
        let g = Highlight::of(&ions, Focus::Glyph(Glyph::StubHook { peptide: 1, n_term: false, residue: 8, stub: 2 }));
        assert_eq!(g.peaks.iter().copied().collect::<Vec<_>>(), vec![y5]);
        // a triangle: the first peptide with the A stub, whole
        let a_stub = Glyph::Stub { peptide: 0, stub: 0, link: None };
        let t = Highlight::of(&ions, Focus::Glyph(a_stub));
        assert_eq!(t.peaks.len(), 1);
        assert_eq!(a.peaks[*t.peaks.first().unwrap()].labels[0].text, "P+A");
        assert_eq!(t.ranges(0), vec![0..12]);
        assert!(t.ranges(1).is_empty());
    }

    /// A fragment holding the link residue (the linker intact) carries the
    /// whole partner; one short of it does not
    #[test]
    fn intact_linker_fragments_hold_the_partner() {
        let (_, a) = model::test::cleaved_dsso();
        // GLEK|VSTAAMQR, linked at K4 (index 3) to TDYHFPNEWK|LIR at K10 (index 9)
        assert_eq!(held(&a, 0, 0..4), vec![(0, 0..4), (1, 0..13)]);
        assert_eq!(held(&a, 0, 0..3), vec![(0, 0..3)]);
        assert_eq!(held(&a, 1, 9..13), vec![(1, 9..13), (0, 0..12)]);
        assert_eq!(held(&a, 1, 10..13), vec![(1, 10..13)]);
    }

    /// A multivalent species: each triangle picks out its peaks; the parts
    /// freed when one linker broke frame every peptide they hold
    #[test]
    fn multivalent_glyphs_point_at_their_peaks() {
        let (link, a) = model::test::cleaved_multi();
        let ions = ions(&a);
        let stub = |name: &str| link.stubs.iter().position(|s| s.name == name).unwrap();
        let (sa, ss) = (stub("A"), stub("S"));
        // peptide 2 (two links) freed whole with A and S: its triplet, framed whole
        let set = Highlight::of(&ions, Focus::Glyph(Glyph::Stubs { peptide: 1, stubs: stub_mask(&[sa, ss]) }));
        assert!(!set.peaks.is_empty());
        assert!(set.peaks.iter().all(|&p| a.peaks[p].labels.iter().any(|l| l.text.contains("A+S"))));
        assert_eq!(set.ranges(1), vec![0..26]);
        // peptide 1 (one link, XL1) with A: by XL1's bridge
        let one = Highlight::of(&ions, Focus::Glyph(Glyph::Stub { peptide: 0, stub: sa, link: Some(0) }));
        assert_eq!(one.peaks.len(), 1);
        assert_eq!(one.ranges(0), vec![0..12]);
        // {2,3} + A, freed when XL1 broke: by XL1's end on peptide 2, peptides 2 and 3 whole
        let part = ions.iter().find(|i| i.text.starts_with("{2,3}")).expect("the part freed by XL1");
        assert_eq!(part.glyph, Some(Glyph::Stub { peptide: 1, stub: sa, link: Some(0) }));
        let h = Highlight::of(&ions, Focus::Peak(part.peak));
        assert_eq!((h.ranges(0), h.ranges(1), h.ranges(2)), (vec![], vec![0..26], vec![0..13]));
    }

    /// Internal ions frame the residues they hold; immonium ions nothing
    #[test]
    fn internal_ions_frame_their_residues() {
        use model::{Coverage, IonLabel, PeakView};
        let label = |text: &str, series, internal| IonLabel {
            text: text.into(),
            series,
            peptide: 0,
            theoretical_mz: 500.0,
            backbone: None,
            internal,
        };
        let a = Annotation {
            proforma: "PEPTIDEK".into(),
            peaks: vec![
                PeakView { mz: 500.0, intensity: 10.0, labels: vec![label("m3:5", Series::Internal, Some((2, 4)))] },
                PeakView { mz: 600.0, intensity: 10.0, labels: vec![label("I", Series::Immonium, None)] },
            ],
            coverage: vec![Coverage { residues: "PEPTIDEK".chars().map(String::from).collect(), n_term: vec![false; 8], c_term: vec![false; 8] }],
            explained_peaks: 2,
            explained_intensity: 1.0,
            max_intensity: 10.0,
            link: None,
            stub_ions: Vec::new(),
        };
        let ions = ions(&a);
        let h = Highlight::of(&ions, Focus::Peak(0));
        assert_eq!(h.ranges(0), vec![2..5]);
        assert!(h.glyphs.is_empty());
        assert!(h.describe(&a, &ions).ends_with("P3–I5"), "{}", h.describe(&a, &ions));
        assert!(Highlight::of(&ions, Focus::Peak(1)).residues.is_empty());
    }

    #[test]
    fn ranges_are_joined() {
        let h = Highlight { residues: vec![(0, 4..6), (0, 0..2), (0, 1..3), (1, 0..5), (0, 3..4)], ..Default::default() };
        assert_eq!(h.ranges(0), vec![0..6]);
        let h = Highlight { residues: vec![(0, 0..2), (0, 3..5)], ..Default::default() };
        assert_eq!(h.ranges(0), vec![0..2, 3..5]);
    }
}
