//! Where a modification may sit, beyond its residue (named modifications,
//! adapted from Sage Plus, MIT): a position rule on the peptide, and a
//! PROSITE-style motif matched against the protein around the site.
//!
//! In a configuration, both go on the object form of a modification:
//!
//! ```json
//! "variable_mods": {
//!   "K": [{"mass": 114.042927, "name": "UNIMOD:121", "position": "not_c_term"}],
//!   "N": [{"mass": 0.984016, "name": "UNIMOD:7", "motif": "N-{P}-[ST]"}],
//!   "S": [{"mass": 79.966331, "name": "UNIMOD:21", "neutral_losses": [97.976896]}]
//! }
//! ```
//!
//! or in the key: `internal_residue:K`, `not_n_term:K`, `not_c_term:K`,
//! `motif:N-{P}-[ST]` (one entry per residue the motif's site allows). The
//! terminal sites have names too (`first_residue:K` is `^K`, `last_residue:K`
//! `$K`, `protein_first:M` `[M`, `protein_last:K` `]K`, `peptide_n_term` `^`,
//! and so on), next to the symbols, so upstream configurations load as before.
//!
//! Decoys are reversed after modification, so they keep their target's sites.
//! A motif is matched against the first protein of a peptide (with the
//! residues around it); a peptide shared by proteins that differ around the
//! site follows that protein.

use crate::mass::VALID_AA;
use crate::modification::{ModEntry, ModId, ModificationSpecificity};
use serde::{de, Deserialize, Deserializer, Serialize, Serializer};
use std::collections::HashMap;

/// Residues kept on each side of a peptide for motif matching
pub const MOTIF_FLANK: usize = 8;

/// Where on the peptide a modified residue may be
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SitePosition {
    #[default]
    Anywhere,
    /// Not the peptide's first residue
    NotNTerm,
    /// Not the peptide's last residue (e.g. a modified lysine that blocks
    /// trypsin: GlyGly, acetyl)
    NotCTerm,
    /// Neither
    Internal,
}

impl SitePosition {
    pub fn is_anywhere(&self) -> bool {
        *self == SitePosition::Anywhere
    }

    /// May the residue at `idx` of a peptide of `len` residues carry it?
    pub fn allows(self, idx: usize, len: usize) -> bool {
        match self {
            SitePosition::Anywhere => true,
            SitePosition::NotNTerm => idx != 0,
            SitePosition::NotCTerm => idx + 1 != len,
            SitePosition::Internal => idx != 0 && idx + 1 != len,
        }
    }
}

/// One position of a motif
#[derive(Clone, Debug, PartialEq, Eq)]
enum Element {
    OneOf(Vec<u8>),
    NoneOf(Vec<u8>),
}

impl Element {
    fn matches(&self, residue: u8) -> bool {
        match self {
            Element::OneOf(set) => set.contains(&residue),
            Element::NoneOf(set) => !set.contains(&residue),
        }
    }

    fn residues(&self) -> Vec<u8> {
        VALID_AA.iter().copied().filter(|&r| r.is_ascii_uppercase() && self.matches(r)).collect()
    }
}

/// A PROSITE-style motif: elements separated by `-`, each a residue (`N`),
/// any residue (`x`), one of (`[ST]`) or none of (`{P}`), optionally
/// repeated (`x(2)`). The modified residue is the first element, or the one
/// marked with `*` (`R-x-x-[ST]*`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Motif {
    pattern: String,
    elements: Vec<Element>,
    /// Index of the modified residue in `elements`
    site: usize,
}

impl Motif {
    pub fn parse(pattern: &str) -> Result<Motif, String> {
        let bad = |why: &str| format!("motif {pattern:?}: {why}");
        let residues = |s: &str| -> Result<Vec<u8>, String> {
            let set: Vec<u8> = s.bytes().map(|b| b.to_ascii_uppercase()).collect();
            match set.is_empty() || set.iter().any(|b| !VALID_AA.contains(b)) {
                true => Err(bad(&format!("{s:?} is not a set of residues"))),
                false => Ok(set),
            }
        };
        let mut elements = Vec::new();
        let mut site = None;
        for part in pattern.trim().split('-') {
            let mut part = part.trim();
            let marked = part.ends_with('*');
            if marked {
                part = part[..part.len() - 1].trim_end();
            }
            // repetition: x(2)
            let (body, repeat) = match part.find('(') {
                Some(open) if part.ends_with(')') => {
                    let n: usize = part[open + 1..part.len() - 1]
                        .trim()
                        .parse()
                        .map_err(|_| bad(&format!("{part:?} has a bad repetition")))?;
                    (&part[..open], n)
                }
                Some(_) => return Err(bad(&format!("{part:?} has a bad repetition"))),
                None => (part, 1),
            };
            if repeat == 0 || repeat > 2 * MOTIF_FLANK {
                return Err(bad("a repetition must be 1 to 16"));
            }
            let element = if body.eq_ignore_ascii_case("x") {
                Element::NoneOf(Vec::new())
            } else if let Some(inner) = body.strip_prefix('[').and_then(|b| b.strip_suffix(']')) {
                Element::OneOf(residues(inner)?)
            } else if let Some(inner) = body.strip_prefix('{').and_then(|b| b.strip_suffix('}')) {
                Element::NoneOf(residues(inner)?)
            } else if body.len() == 1 {
                Element::OneOf(residues(body)?)
            } else {
                return Err(bad(&format!("{body:?} is not a residue, x, [..] or {{..}}")));
            };
            if marked {
                if site.is_some() || repeat != 1 {
                    return Err(bad("mark one single element as the site with *"));
                }
                site = Some(elements.len());
            }
            elements.extend(std::iter::repeat(element).take(repeat));
        }
        if elements.is_empty() {
            return Err(bad("empty"));
        }
        let site = site.unwrap_or(0);
        if site > MOTIF_FLANK || elements.len() - 1 - site > MOTIF_FLANK {
            return Err(bad(&format!("at most {MOTIF_FLANK} positions on each side of the site")));
        }
        if elements[site].residues().is_empty() {
            return Err(bad("no residue can be the site"));
        }
        Ok(Motif { pattern: pattern.trim().to_string(), elements, site })
    }

    pub fn pattern(&self) -> &str {
        &self.pattern
    }

    /// Residues the modified position can be
    pub fn site_residues(&self) -> Vec<u8> {
        self.elements[self.site].residues()
    }

    /// Does the motif match `context` with its site at `at`? Positions
    /// beyond `context` (past a protein terminus) do not match.
    pub fn matches(&self, context: &[u8], at: usize) -> bool {
        let Some(start) = at.checked_sub(self.site) else { return false };
        if start + self.elements.len() > context.len() {
            return false;
        }
        self.elements.iter().zip(&context[start..]).all(|(e, &r)| e.matches(r))
    }
}

impl Serialize for Motif {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(&self.pattern)
    }
}

impl<'de> Deserialize<'de> for Motif {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let s = String::deserialize(deserializer)?;
        Motif::parse(&s).map_err(de::Error::custom)
    }
}

/// Position rule and motif of one configured modification
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SiteRule {
    pub position: SitePosition,
    pub motif: Option<Motif>,
    /// `site_mode: library`: only where the site library lists this
    /// modification (the peptide's [`SiteContext::library`])
    pub library: Option<ModId>,
}

impl SiteRule {
    pub fn is_any(&self) -> bool {
        self.position.is_anywhere() && self.motif.is_none() && self.library.is_none()
    }

    /// The rule on the peptide alone, without the site library (where site
    /// localization moves a modification)
    pub fn without_library(&self) -> SiteRule {
        SiteRule { library: None, ..self.clone() }
    }

    /// May residue `idx` of `peptide` carry the modification? `context`: the
    /// peptide with up to [`MOTIF_FLANK`] protein residues on each side, and
    /// where the peptide starts in it (else the peptide alone)
    pub fn allows(&self, idx: usize, peptide: &[u8], context: Option<&SiteContext>) -> bool {
        if !self.position.allows(idx, peptide.len()) {
            return false;
        }
        if let Some(id) = self.library {
            if !context.is_some_and(|c| c.library.contains(&(idx, id))) {
                return false;
            }
        }
        match (&self.motif, context) {
            (None, _) => true,
            (Some(motif), Some(c)) => motif.matches(&c.sequence, c.offset + idx),
            (Some(motif), None) => motif.matches(peptide, idx),
        }
    }
}

/// A peptide in its protein, for motifs that reach past the peptide
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SiteContext {
    /// Up to [`MOTIF_FLANK`] residues before, the peptide, and up to
    /// [`MOTIF_FLANK`] after
    pub sequence: Vec<u8>,
    /// Where the peptide starts in `sequence`
    pub offset: usize,
    /// Residues of the peptide (by index) at which the site library lists a
    /// modification, on any of its proteins
    pub library: Vec<(usize, ModId)>,
}

impl SiteContext {
    /// The context of `peptide` found in `protein` (its first occurrence)
    pub fn find(protein: &[u8], peptide: &[u8]) -> Option<SiteContext> {
        let at = protein.windows(peptide.len().max(1)).position(|w| w == peptide)?;
        let start = at.saturating_sub(MOTIF_FLANK);
        let end = (at + peptide.len() + MOTIF_FLANK).min(protein.len());
        Some(SiteContext { sequence: protein[start..end].to_vec(), offset: at - start, library: Vec::new() })
    }
}

/// A modification key: the residue/terminal symbol or its name, or a rule
/// and residue (`internal_residue:K`, `motif:N-{P}-[ST]`)
pub fn parse_site_key(key: &str) -> Result<Vec<(ModificationSpecificity, SiteRule)>, String> {
    let key = key.trim();
    if let Some(pattern) = key.strip_prefix("motif:") {
        let motif = Motif::parse(pattern)?;
        return Ok(motif
            .site_residues()
            .into_iter()
            .map(|r| {
                (
                    ModificationSpecificity::Residue(r),
                    SiteRule { position: SitePosition::Anywhere, motif: Some(motif.clone()), library: None },
                )
            })
            .collect());
    }
    let positional = [
        ("internal_residue:", SitePosition::Internal),
        ("not_n_term:", SitePosition::NotNTerm),
        ("not_c_term:", SitePosition::NotCTerm),
    ];
    for (prefix, position) in positional {
        if let Some(residue) = key.strip_prefix(prefix) {
            let specificity: ModificationSpecificity = residue.parse().map_err(|e| format!("{e}"))?;
            if !matches!(specificity, ModificationSpecificity::Residue(_)) {
                return Err(format!("{key}: give one residue after {prefix}"));
            }
            return Ok(vec![(specificity, SiteRule { position, motif: None, library: None })]);
        }
    }
    let specificity: ModificationSpecificity = key.parse().map_err(|e| format!("{e}"))?;
    Ok(vec![(specificity, SiteRule::default())])
}

/// `database.static_mods`, with named and rule keys
pub fn named_static_mods<'de, D>(deserializer: D) -> Result<Option<HashMap<ModificationSpecificity, ModEntry>>, D::Error>
where
    D: Deserializer<'de>,
{
    let Some(raw) = Option::<HashMap<String, ModEntry>>::deserialize(deserializer)? else { return Ok(None) };
    let mut out = HashMap::new();
    let mut keys: Vec<_> = raw.into_iter().collect();
    keys.sort_by(|a, b| a.0.cmp(&b.0));
    for (key, entry) in keys {
        for (specificity, rule) in parse_site_key(&key).map_err(de::Error::custom)? {
            let entry = entry.clone().with_rule(&rule).map_err(de::Error::custom)?;
            if out.insert(specificity, entry).is_some() {
                return Err(de::Error::custom(format!("static modification on {specificity} is given twice ({key})")));
            }
        }
    }
    Ok(Some(out))
}

/// `database.variable_mods`, with named and rule keys
pub fn named_variable_mods<'de, D>(
    deserializer: D,
) -> Result<Option<HashMap<ModificationSpecificity, Vec<ModEntry>>>, D::Error>
where
    D: Deserializer<'de>,
{
    let Some(raw) = Option::<HashMap<String, Vec<ModEntry>>>::deserialize(deserializer)? else { return Ok(None) };
    let mut out: HashMap<ModificationSpecificity, Vec<ModEntry>> = HashMap::new();
    let mut keys: Vec<_> = raw.into_iter().collect();
    keys.sort_by(|a, b| a.0.cmp(&b.0));
    for (key, entries) in keys {
        for (specificity, rule) in parse_site_key(&key).map_err(de::Error::custom)? {
            for entry in &entries {
                out.entry(specificity).or_default().push(entry.clone().with_rule(&rule).map_err(de::Error::custom)?);
            }
        }
    }
    Ok(Some(out))
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn motifs() {
        let sequon = Motif::parse("N-{P}-[ST]").unwrap();
        assert_eq!(sequon.site_residues(), b"N");
        assert!(sequon.matches(b"AANGTA", 2));
        assert!(!sequon.matches(b"AANPTA", 2));
        assert!(!sequon.matches(b"AANGAA", 2));
        assert!(!sequon.matches(b"AANG", 2), "past the end");
        let basophilic = Motif::parse("R-x-x-[ST]*").unwrap();
        assert_eq!(basophilic.site_residues(), b"ST");
        assert!(basophilic.matches(b"RAAS", 3));
        assert!(!basophilic.matches(b"RAAS", 2));
        assert!(!basophilic.matches(b"AAS", 2), "before the start");
        let repeated = Motif::parse("[ST]-x(2)-E").unwrap();
        assert!(repeated.matches(b"SAAE", 0));
        for bad in ["", "N-", "B?", "x(0)", "N*-S*", "x(2)*-S", "[]", "x(20)"] {
            assert!(Motif::parse(bad).is_err(), "{bad}");
        }
    }

    #[test]
    fn positions() {
        use SitePosition::*;
        assert!(Anywhere.allows(0, 5) && Anywhere.allows(4, 5));
        assert!(!NotCTerm.allows(4, 5) && NotCTerm.allows(0, 5));
        assert!(!NotNTerm.allows(0, 5) && NotNTerm.allows(4, 5));
        assert!(!Internal.allows(0, 5) && !Internal.allows(4, 5) && Internal.allows(2, 5));
    }

    #[test]
    fn keys() {
        use ModificationSpecificity::*;
        assert_eq!(parse_site_key("K").unwrap(), vec![(Residue(b'K'), SiteRule::default())]);
        assert_eq!(parse_site_key("first_residue:Q").unwrap()[0].0, PeptideN(Some(b'Q')));
        assert_eq!(parse_site_key("protein_last:K").unwrap()[0].0, ProteinC(Some(b'K')));
        assert_eq!(parse_site_key("peptide_n_term").unwrap()[0].0, PeptideN(None));
        let internal = parse_site_key("internal_residue:K").unwrap();
        assert_eq!(internal[0], (Residue(b'K'), SiteRule { position: SitePosition::Internal, motif: None, library: None }));
        let motif = parse_site_key("motif:[ST]-P").unwrap();
        assert_eq!(motif.iter().map(|(s, _)| *s).collect::<Vec<_>>(), vec![Residue(b'S'), Residue(b'T')]);
        assert!(parse_site_key("internal_residue:^").is_err());
        assert!(parse_site_key("elsewhere:K").is_err());
    }

    #[test]
    fn context_reaches_past_the_peptide() {
        let protein = b"MAAAKNGTRLLLR";
        let c = SiteContext::find(protein, b"NGTR").unwrap();
        let rule = SiteRule {
            position: SitePosition::Anywhere,
            motif: Some(Motif::parse("N-{P}-[ST]").unwrap()),
            library: None,
        };
        assert!(rule.allows(0, b"NGTR", Some(&c)));
        // the site needs the residues after the peptide: AAKN|GT
        let c = SiteContext::find(protein, b"AAAKN").unwrap();
        assert!(rule.allows(4, b"AAAKN", Some(&c)));
        assert!(!rule.allows(4, b"AAAKN", None), "the peptide alone does not hold the motif");
    }
}
