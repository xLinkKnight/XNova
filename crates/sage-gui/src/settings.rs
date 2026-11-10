//! The GUI's model of a Sage search, converted to and from the JSON
//! configuration file that `sage.exe` reads.
//!
//! Keys the GUI does not edit (e.g. `bruker_config`, LFQ peak scoring) are kept
//! when a configuration is loaded and written back unchanged.

use serde_json::{json, Map, Value};

/// Residues Sage accepts in modification targets and enzyme rules
pub const VALID_AA: &str = "ACDEFGHIKLMNPQRSTVWYUO";

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum TolUnit {
    Ppm,
    Da,
}

impl TolUnit {
    pub fn key(self) -> &'static str {
        match self {
            TolUnit::Ppm => "ppm",
            TolUnit::Da => "da",
        }
    }
}

#[derive(Clone, Debug)]
pub struct Tolerance {
    pub unit: TolUnit,
    pub lo: f64,
    pub hi: f64,
}

impl Tolerance {
    fn ppm(lo: f64, hi: f64) -> Self {
        Self {
            unit: TolUnit::Ppm,
            lo,
            hi,
        }
    }

    fn to_json(&self) -> Value {
        let mut map = Map::new();
        map.insert(self.unit.key().into(), json!([self.lo, self.hi]));
        Value::Object(map)
    }

    fn from_json(v: &Value) -> Option<Self> {
        let (key, range) = v.as_object()?.iter().next()?;
        let unit = match key.as_str() {
            "ppm" => TolUnit::Ppm,
            "da" => TolUnit::Da,
            _ => return None,
        };
        let range = range.as_array()?;
        Some(Self {
            unit,
            lo: range.first()?.as_f64()?,
            hi: range.get(1)?.as_f64()?,
        })
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub enum ModKind {
    Fixed,
    #[default]
    Variable,
}

#[derive(Clone, Debug, Default)]
pub struct Modification {
    pub target: String,
    pub mass: f64,
    /// Optional ProForma name, e.g. `UNIMOD:35`
    pub name: String,
    pub kind: ModKind,
    /// Variable mods only: maximum occurrences per peptide, 0 = no limit
    pub max_count: u32,
    /// Where on the peptide the residue may be
    pub position: SitePosition,
    /// PROSITE-style motif around the site, e.g. `N-{P}-[ST]`; empty: none
    pub motif: String,
    /// Neutral losses of fragments containing the site, Da, comma-separated
    pub losses: String,
    /// Variable mods only: indexed, at the site library's sites, or a mass offset
    pub search: ModSearch,
}

pub use sage_core::site_rules::SitePosition;

/// How a variable modification is searched
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub enum ModSearch {
    /// Every modified form indexed
    #[default]
    Indexed,
    /// Indexed at the site library's sites only (`site_mode: library`)
    Library,
    /// Tried at search time (`search_mode: mass_offset`)
    Offset,
}

/// The search choices offered for a variable modification
pub const MOD_SEARCHES: [(ModSearch, &str); 3] =
    [(ModSearch::Indexed, "indexed"), (ModSearch::Library, "library sites"), (ModSearch::Offset, "mass offset")];

pub fn search_label(s: ModSearch) -> &'static str {
    MOD_SEARCHES.iter().find(|(x, _)| *x == s).map_or("indexed", |(_, l)| *l)
}

/// The position choices offered for a modification
pub const SITE_POSITIONS: [(SitePosition, &str); 4] = [
    (SitePosition::Anywhere, "anywhere"),
    (SitePosition::NotNTerm, "not N-term"),
    (SitePosition::NotCTerm, "not C-term"),
    (SitePosition::Internal, "internal"),
];

pub fn position_label(p: SitePosition) -> &'static str {
    SITE_POSITIONS.iter().find(|(x, _)| *x == p).map_or("anywhere", |(_, l)| *l)
}

/// The neutral losses written in `losses` (Da), or what is wrong with them
pub fn parse_losses(losses: &str) -> Result<Vec<f64>, String> {
    losses
        .split([',', ';', ' '])
        .map(str::trim)
        .filter(|x| !x.is_empty())
        .map(|x| match x.parse::<f64>() {
            Ok(v) if v.is_finite() && v > 0.0 => Ok(v),
            _ => Err(format!("`{x}` is not a positive mass")),
        })
        .collect()
}

/// What is wrong with a modification's position, motif or neutral losses
pub fn rule_problem(m: &Modification) -> Option<String> {
    let target = m.target.trim();
    let residue = target.chars().last().filter(|c| c.is_ascii_uppercase());
    if m.position != SitePosition::Anywhere && (target.len() != 1 || residue.is_none()) {
        return Some(format!("`{target}`: a position (not N-/C-term, internal) needs a single residue as the target"));
    }
    if !m.motif.trim().is_empty() {
        match sage_core::site_rules::Motif::parse(m.motif.trim()) {
            Err(e) => return Some(format!("`{target}`: {e}")),
            Ok(motif) => {
                let Some(r) = residue else {
                    return Some(format!("`{target}`: a motif needs a residue as the target"));
                };
                if !motif.site_residues().contains(&(r as u8)) {
                    return Some(format!("`{target}`: the site of motif {} is never {r}", motif.pattern()));
                }
            }
        }
    }
    if m.kind == ModKind::Variable {
        match m.search {
            ModSearch::Offset if !m.motif.trim().is_empty() => {
                return Some(format!("`{target}`: a mass offset is placed on the peptide alone and cannot have a motif"));
            }
            ModSearch::Offset if m.mass == 0.0 => return Some(format!("`{target}`: a mass offset needs a mass")),
            ModSearch::Library if target.len() != 1 || residue.is_none() => {
                return Some(format!("`{target}`: library sites need a single residue as the target"));
            }
            _ => {}
        }
    }
    parse_losses(&m.losses).err().map(|e| format!("`{target}`: neutral loss {e}"))
}

/// A common modification offered in the "Add" menu
pub struct Preset {
    pub label: &'static str,
    pub target: &'static str,
    pub mass: f64,
    pub name: &'static str,
    pub kind: ModKind,
}

pub const PRESETS: &[Preset] = &[
    Preset { label: "Carbamidomethyl (C)", target: "C", mass: 57.021464, name: "UNIMOD:4", kind: ModKind::Fixed },
    Preset { label: "Propionamide (C)", target: "C", mass: 71.037114, name: "UNIMOD:24", kind: ModKind::Fixed },
    Preset { label: "Oxidation (M)", target: "M", mass: 15.994915, name: "UNIMOD:35", kind: ModKind::Variable },
    Preset { label: "Acetyl (protein N-term)", target: "[", mass: 42.010565, name: "UNIMOD:1", kind: ModKind::Variable },
    Preset { label: "Phospho (S)", target: "S", mass: 79.966331, name: "UNIMOD:21", kind: ModKind::Variable },
    Preset { label: "Phospho (T)", target: "T", mass: 79.966331, name: "UNIMOD:21", kind: ModKind::Variable },
    Preset { label: "Phospho (Y)", target: "Y", mass: 79.966331, name: "UNIMOD:21", kind: ModKind::Variable },
    Preset { label: "Deamidated (N)", target: "N", mass: 0.984016, name: "UNIMOD:7", kind: ModKind::Variable },
    Preset { label: "Deamidated (Q)", target: "Q", mass: 0.984016, name: "UNIMOD:7", kind: ModKind::Variable },
    Preset { label: "Gln->pyro-Glu (peptide N-term Q)", target: "^Q", mass: -17.026549, name: "UNIMOD:28", kind: ModKind::Variable },
    Preset { label: "Glu->pyro-Glu (peptide N-term E)", target: "^E", mass: -18.010565, name: "UNIMOD:27", kind: ModKind::Variable },
    Preset { label: "GlyGly / ubiquitin remnant (K)", target: "K", mass: 114.042927, name: "UNIMOD:121", kind: ModKind::Variable },
    Preset { label: "Methyl (K)", target: "K", mass: 14.01565, name: "UNIMOD:34", kind: ModKind::Variable },
    Preset { label: "Methyl (R)", target: "R", mass: 14.01565, name: "UNIMOD:34", kind: ModKind::Variable },
    Preset { label: "Dimethyl (K)", target: "K", mass: 28.0313, name: "UNIMOD:36", kind: ModKind::Variable },
    Preset { label: "Acetyl (K)", target: "K", mass: 42.010565, name: "UNIMOD:1", kind: ModKind::Variable },
    Preset { label: "TMT 6/10/11-plex (K)", target: "K", mass: 229.162932, name: "UNIMOD:737", kind: ModKind::Fixed },
    Preset { label: "TMT 6/10/11-plex (peptide N-term)", target: "^", mass: 229.162932, name: "UNIMOD:737", kind: ModKind::Fixed },
    Preset { label: "TMTpro 16/18-plex (K)", target: "K", mass: 304.207146, name: "UNIMOD:2016", kind: ModKind::Fixed },
    Preset { label: "TMTpro 16/18-plex (peptide N-term)", target: "^", mass: 304.207146, name: "UNIMOD:2016", kind: ModKind::Fixed },
    Preset { label: "DSS/BS3 mono-link, hydrolyzed (K)", target: "K", mass: 156.078644, name: "", kind: ModKind::Variable },
];

impl From<&Preset> for Modification {
    fn from(p: &Preset) -> Self {
        Self {
            target: p.target.into(),
            mass: p.mass,
            name: p.name.into(),
            kind: p.kind,
            max_count: 0,
            // phospho: the phosphate leaves S and T fragments as H3PO4
            losses: if p.name == "UNIMOD:21" && matches!(p.target, "S" | "T") { "97.976896".into() } else { String::new() },
            ..Default::default()
        }
    }
}

/// Check a modification target such as `C`, `^`, `^Q`, `$`, `[`, `]M`
pub fn check_target(target: &str) -> Result<(), String> {
    let mut chars = target.chars();
    let first = chars.next().ok_or("target is empty")?;
    let rest: Vec<char> = chars.collect();
    let aa_ok = |c: char| VALID_AA.contains(c);
    match first {
        '^' | '$' | '[' | ']' => match rest[..] {
            [] => Ok(()),
            [c] if aa_ok(c) => Ok(()),
            _ => Err(format!("`{target}`: a terminus may be followed by at most one amino acid")),
        },
        c if rest.is_empty() && aa_ok(c) => Ok(()),
        _ => Err(format!(
            "`{target}`: use one amino acid ({VALID_AA}) or a terminus (^ $ [ ]) with an optional amino acid"
        )),
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Reagent {
    Dss,
    Bs3,
    Dsg,
    Bs2g,
    Dsso,
    Dsbu,
    Edc,
    Dmtmm,
    Clixlink,
    Custom,
}

impl Reagent {
    pub const ALL: [Reagent; 10] = [
        Reagent::Dss,
        Reagent::Bs3,
        Reagent::Dsg,
        Reagent::Bs2g,
        Reagent::Dsso,
        Reagent::Dsbu,
        Reagent::Edc,
        Reagent::Dmtmm,
        Reagent::Clixlink,
        Reagent::Custom,
    ];

    pub fn key(self) -> &'static str {
        match self {
            Reagent::Dss => "dss",
            Reagent::Bs3 => "bs3",
            Reagent::Dsg => "dsg",
            Reagent::Bs2g => "bs2g",
            Reagent::Dsso => "dsso",
            Reagent::Dsbu => "dsbu",
            Reagent::Edc => "edc",
            Reagent::Dmtmm => "dmtmm",
            Reagent::Clixlink => "clixlink",
            Reagent::Custom => "custom",
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Reagent::Dss => "DSS (138.068 Da, non-cleavable)",
            Reagent::Bs3 => "BS3 (138.068 Da, non-cleavable)",
            Reagent::Dsg => "DSG (96.021 Da, non-cleavable)",
            Reagent::Bs2g => "BS2G (96.021 Da, non-cleavable)",
            Reagent::Dsso => "DSSO (158.004 Da, MS-cleavable)",
            Reagent::Dsbu => "DSBU (196.085 Da, MS-cleavable)",
            Reagent::Edc => "EDC (zero-length, -18.011 Da)",
            Reagent::Dmtmm => "DMTMM (zero-length, -18.011 Da)",
            Reagent::Clixlink => "cliXlink (267.057 Da, MS-cleavable, asymmetric)",
            Reagent::Custom => "Custom mass",
        }
    }

    pub fn is_cleavable(self) -> bool {
        self.stubs().is_some()
    }

    /// Linker mass of a preset (Da), as in `sage_core::crosslink::Reagent::mass`
    pub fn mass(self) -> Option<f64> {
        match self {
            Reagent::Dss | Reagent::Bs3 => Some(138.06808),
            Reagent::Dsg | Reagent::Bs2g => Some(96.02113),
            Reagent::Dsso => Some(158.00376),
            Reagent::Dsbu => Some(196.08479),
            Reagent::Edc | Reagent::Dmtmm => Some(-18.010565),
            Reagent::Clixlink => Some(267.05653),
            Reagent::Custom => None,
        }
    }

    /// Preset `(short, long, alt)` stub masses of an MS-cleavable reagent
    /// (for an asymmetric one, its first arm)
    pub fn stubs(self) -> Option<(f64, f64, &'static [f64])> {
        match self {
            Reagent::Dsso | Reagent::Clixlink => Some((54.01056, 103.9932, &[85.98264])),
            Reagent::Dsbu => Some((85.05276, 111.03203, &[])),
            _ => None,
        }
    }

    /// Asymmetric cleavable reagents: the stubs of the second arm
    pub fn second_arm(self) -> Option<(f64, f64, &'static [f64])> {
        match self {
            // the alkyne-tagged arm: alkene / sulfenic acid / thiol
            Reagent::Clixlink => Some((163.06333, 213.04596, &[195.0354])),
            _ => None,
        }
    }

    /// Default `(site1, site2)` for reagents with two different reactive ends
    pub fn sites(self) -> Option<(Site, Site)> {
        match self {
            Reagent::Edc | Reagent::Dmtmm => Some((
                Site { residues: "K".into(), protein_nterm: true, protein_cterm: false },
                Site { residues: "DE".into(), protein_nterm: false, protein_cterm: true },
            )),
            _ => None,
        }
    }
}

/// What one end of a heterobifunctional linker reacts with
#[derive(Clone, Debug, PartialEq)]
pub struct Site {
    pub residues: String,
    pub protein_nterm: bool,
    pub protein_cterm: bool,
}

impl Site {
    /// Sage's notation: residue letters, `[` protein N-terminus, `]` protein C-terminus
    pub fn to_spec(&self) -> String {
        let mut s = self.residues.clone();
        if self.protein_nterm {
            s.push('[');
        }
        if self.protein_cterm {
            s.push(']');
        }
        s
    }

    pub fn from_spec(spec: &str) -> Self {
        Site {
            residues: spec.chars().filter(char::is_ascii_alphabetic).collect(),
            protein_nterm: spec.contains('['),
            protein_cterm: spec.contains(']'),
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Backend {
    Cpu,
    Reflect,
    Gpu,
    Auto,
}

impl Backend {
    pub const ALL: [Backend; 4] = [Backend::Cpu, Backend::Reflect, Backend::Gpu, Backend::Auto];

    pub fn key(self) -> &'static str {
        match self {
            Backend::Cpu => "cpu",
            Backend::Reflect => "reflect",
            Backend::Gpu => "gpu",
            Backend::Auto => "auto",
        }
    }
}

#[derive(Clone, Debug)]
pub struct Crosslink {
    pub enabled: bool,
    pub reagent: Reagent,
    pub custom_mass: f64,
    pub custom_cleavable: bool,
    /// Fragment (stub) masses of an MS-cleavable linker, preset or custom
    pub stub_short: f64,
    pub stub_long: f64,
    /// Extra stub masses, comma separated
    pub stub_alt: String,
    /// Custom cleavable linkers: the two arms leave different stubs
    pub asymmetric: bool,
    /// Asymmetric linkers: stubs of the second arm
    pub stub2_short: f64,
    pub stub2_long: f64,
    pub stub2_alt: String,
    /// Homobifunctional: residues both ends react with
    pub residues: String,
    pub protein_nterm: bool,
    /// Heterobifunctional: each end reacts with a different group
    pub two_sites: bool,
    pub site1: Site,
    pub site2: Site,
    pub allow_cterm_residue: bool,
    pub loop_links: bool,
    pub candidates: u32,
    pub min_chain_peaks: u32,
    pub require_complete_doublets: bool,
    pub backend: Backend,
    pub chain_survivors: u32,
    /// Also run a linear peptide search (results.sage.tsv)
    pub linear_search: bool,
    /// Non-cleavable: matched linker-carrying ions required per chain
    pub min_linked_peaks: u32,
    /// Also search multivalent species (several linkers), up to `valency` linkers
    pub multivalent: bool,
    pub valency: u8,
}

impl Default for Crosslink {
    fn default() -> Self {
        let (site1, site2) = Reagent::Edc.sites().unwrap();
        Self {
            enabled: false,
            reagent: Reagent::Dss,
            custom_mass: 138.06808,
            custom_cleavable: false,
            stub_short: 54.01056,
            stub_long: 103.9932,
            stub_alt: String::new(),
            asymmetric: false,
            stub2_short: 163.06333,
            stub2_long: 213.04596,
            stub2_alt: String::new(),
            residues: "K".into(),
            two_sites: false,
            site1,
            site2,
            protein_nterm: true,
            allow_cterm_residue: false,
            loop_links: true,
            candidates: 30,
            min_chain_peaks: 3,
            require_complete_doublets: false,
            backend: Backend::Cpu,
            chain_survivors: 20000,
            linear_search: true,
            min_linked_peaks: 2,
            multivalent: false,
            valency: 2,
        }
    }
}

impl Crosslink {
    /// Is the selected linker MS-cleavable (a cleavable preset, or custom with stubs)?
    pub fn is_cleavable(&self) -> bool {
        self.reagent.is_cleavable() || (self.reagent == Reagent::Custom && self.custom_cleavable)
    }

    /// Mass the linker adds to the peptide pair (Da)
    pub fn linker_mass(&self) -> f64 {
        self.reagent.mass().unwrap_or(self.custom_mass)
    }

    /// Does the linker leave different stubs on its two arms?
    pub fn has_second_arm(&self) -> bool {
        match self.reagent {
            Reagent::Custom => self.custom_cleavable && self.asymmetric,
            r => r.second_arm().is_some(),
        }
    }

    /// Load the selected preset's stub masses (no change for non-cleavable / custom)
    pub fn reset_stubs(&mut self) {
        let list = |alt: &[f64]| alt.iter().map(f64::to_string).collect::<Vec<_>>().join(", ");
        if let Some((short, long, alt)) = self.reagent.stubs() {
            self.stub_short = short;
            self.stub_long = long;
            self.stub_alt = list(alt);
        }
        if let Some((short, long, alt)) = self.reagent.second_arm() {
            self.stub2_short = short;
            self.stub2_long = long;
            self.stub2_alt = list(alt);
        }
    }

    /// Have the stub masses of a cleavable preset been changed from its defaults?
    pub fn stubs_edited(&self) -> bool {
        let Some(first) = self.reagent.stubs() else {
            return false;
        };
        let same = |a: f64, b: f64| (a - b).abs() < 1e-6;
        let arm_same = |(short, long, alt): (f64, f64, &[f64]), (s, l, a): (f64, f64, &str)| {
            same(s, short)
                && same(l, long)
                && parse_list(a).is_ok_and(|v| v.len() == alt.len() && v.iter().zip(alt).all(|(x, y)| same(*x, *y)))
        };
        let first_same = arm_same(first, (self.stub_short, self.stub_long, &self.stub_alt));
        let second_same = self
            .reagent
            .second_arm()
            .is_none_or(|second| arm_same(second, (self.stub2_short, self.stub2_long, &self.stub2_alt)));
        !(first_same && second_same)
    }

    fn stubs_json(&self) -> Value {
        let mut stubs = json!({
            "short": self.stub_short,
            "long": self.stub_long,
            "alt": parse_list(&self.stub_alt).unwrap_or_default(),
        });
        if self.has_second_arm() {
            stubs["other"] = json!({
                "short": self.stub2_short,
                "long": self.stub2_long,
                "alt": parse_list(&self.stub2_alt).unwrap_or_default(),
            });
        }
        stubs
    }

    fn to_json(&self) -> Value {
        let mut m = Map::new();
        if self.reagent == Reagent::Custom {
            m.insert("mass".into(), json!(self.custom_mass));
            if self.custom_cleavable {
                m.insert("cleavable".into(), self.stubs_json());
            }
        } else {
            m.insert("reagent".into(), json!(self.reagent.key()));
            // Preset stubs are implied by the reagent; write them only when edited
            if self.stubs_edited() {
                m.insert("cleavable".into(), self.stubs_json());
            }
        }
        if self.two_sites {
            m.insert("site1".into(), json!(self.site1.to_spec()));
            m.insert("site2".into(), json!(self.site2.to_spec()));
        } else {
            m.insert("residues".into(), json!(self.residues));
            m.insert("protein_nterm".into(), json!(self.protein_nterm));
        }
        m.insert("allow_cterm_residue".into(), json!(self.allow_cterm_residue));
        m.insert("loop_links".into(), json!(self.loop_links));
        m.insert("candidates".into(), json!(self.candidates));
        m.insert("min_chain_peaks".into(), json!(self.min_chain_peaks));
        m.insert(
            "require_complete_doublets".into(),
            json!(self.require_complete_doublets),
        );
        m.insert("backend".into(), json!(self.backend.key()));
        m.insert("chain_survivors".into(), json!(self.chain_survivors));
        // always written: Sage turns the linear search on when it is absent
        m.insert("linear_search".into(), json!(self.linear_search));
        // Sage's default (2) is left out, so configs still run in an older sage.exe
        if self.min_linked_peaks != 2 {
            m.insert("min_linked_peaks".into(), json!(self.min_linked_peaks));
        }
        // only when on (off is Sage's default)
        if self.multivalent {
            m.insert("multivalent".into(), json!(self.valency));
        }
        Value::Object(m)
    }

    fn from_json(v: &Map<String, Value>) -> Self {
        let mut c = Crosslink {
            enabled: true,
            ..Default::default()
        };
        if let Some(r) = v.get("reagent").and_then(Value::as_str) {
            if let Some(r) = Reagent::ALL.iter().find(|x| x.key() == r) {
                c.reagent = *r;
            }
        } else if let Some(mass) = v.get("mass").and_then(Value::as_f64) {
            c.reagent = Reagent::Custom;
            c.custom_mass = mass;
        }
        c.reset_stubs();
        if let Some(stubs) = v.get("cleavable").and_then(Value::as_object) {
            c.custom_cleavable = c.reagent == Reagent::Custom;
            c.stub_short = f64_or(stubs.get("short"), c.stub_short);
            c.stub_long = f64_or(stubs.get("long"), c.stub_long);
            let list = |v: Option<&Value>| {
                v.and_then(Value::as_array)
                    .map(|a| a.iter().filter_map(Value::as_f64).map(|x| x.to_string()).collect::<Vec<_>>().join(", "))
                    .unwrap_or_default()
            };
            // Sage's default `alt` is empty, not the preset's
            c.stub_alt = list(stubs.get("alt"));
            if let Some(other) = stubs.get("other").and_then(Value::as_object) {
                c.asymmetric = true;
                c.stub2_short = f64_or(other.get("short"), c.stub2_short);
                c.stub2_long = f64_or(other.get("long"), c.stub2_long);
                c.stub2_alt = list(other.get("alt"));
            }
        }
        if let Some(r) = v.get("residues").and_then(Value::as_str) {
            c.residues = r.into();
        }
        c.protein_nterm = bool_or(v.get("protein_nterm"), c.protein_nterm);
        if let (Some(s1), Some(s2)) = (
            v.get("site1").and_then(Value::as_str),
            v.get("site2").and_then(Value::as_str),
        ) {
            c.two_sites = true;
            c.site1 = Site::from_spec(s1);
            c.site2 = Site::from_spec(s2);
        } else if let Some((s1, s2)) = c.reagent.sites() {
            c.two_sites = true;
            c.site1 = s1;
            c.site2 = s2;
        }
        c.allow_cterm_residue = bool_or(v.get("allow_cterm_residue"), c.allow_cterm_residue);
        c.loop_links = bool_or(v.get("loop_links"), c.loop_links);
        c.candidates = u32_or(v.get("candidates"), c.candidates);
        c.min_chain_peaks = u32_or(v.get("min_chain_peaks"), c.min_chain_peaks);
        c.require_complete_doublets =
            bool_or(v.get("require_complete_doublets"), c.require_complete_doublets);
        if let Some(b) = v.get("backend").and_then(Value::as_str) {
            if let Some(b) = Backend::ALL.iter().find(|x| x.key() == b) {
                c.backend = *b;
            }
        }
        c.chain_survivors = u32_or(v.get("chain_survivors"), c.chain_survivors);
        c.linear_search = bool_or(v.get("linear_search"), c.linear_search);
        c.min_linked_peaks = u32_or(v.get("min_linked_peaks"), c.min_linked_peaks);
        if let Some(n) = v.get("multivalent").and_then(Value::as_u64) {
            c.multivalent = true;
            c.valency = (n as u8).clamp(2, sage_core::crosslink::MAX_VALENCY);
        }
        c
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Tmt {
    None,
    Tmt6,
    Tmt10,
    Tmt11,
    Tmt16,
    Tmt18,
}

impl Tmt {
    pub const ALL: [Tmt; 6] = [Tmt::None, Tmt::Tmt6, Tmt::Tmt10, Tmt::Tmt11, Tmt::Tmt16, Tmt::Tmt18];

    /// Name used in the Sage JSON, `None` for no TMT quantification
    pub fn key(self) -> Option<&'static str> {
        match self {
            Tmt::None => None,
            Tmt::Tmt6 => Some("Tmt6"),
            Tmt::Tmt10 => Some("Tmt10"),
            Tmt::Tmt11 => Some("Tmt11"),
            Tmt::Tmt16 => Some("Tmt16"),
            Tmt::Tmt18 => Some("Tmt18"),
        }
    }

    pub fn label(self) -> &'static str {
        self.key().unwrap_or("None")
    }
}

/// SILAC: no labels, light + heavy, or light + medium + heavy; or dimethyl
/// labelling on K and peptide N-termini (light + medium, or all three)
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Silac {
    Off,
    Duplex,
    Triplex,
    DimethylDuplex,
    DimethylTriplex,
}

impl Silac {
    pub const ALL: [Silac; 5] =
        [Silac::Off, Silac::Duplex, Silac::Triplex, Silac::DimethylDuplex, Silac::DimethylTriplex];

    pub fn label(self) -> &'static str {
        match self {
            Silac::Off => "Off",
            Silac::Duplex => "SILAC duplex (light, heavy)",
            Silac::Triplex => "SILAC triplex (light, medium, heavy)",
            Silac::DimethylDuplex => "Dimethyl duplex (light, medium)",
            Silac::DimethylTriplex => "Dimethyl triplex (light, medium, heavy)",
        }
    }

    /// Dimethyl labelling (K and N-termini; light carries a label too)
    pub fn is_dimethyl(self) -> bool {
        matches!(self, Silac::DimethylDuplex | Silac::DimethylTriplex)
    }

    /// SILAC on residues, with editable K / R labels
    pub fn is_residue_silac(self) -> bool {
        matches!(self, Silac::Duplex | Silac::Triplex)
    }

    /// Label states every peptide is searched in
    pub fn states(self) -> usize {
        match self {
            Silac::Off => 1,
            Silac::Duplex | Silac::DimethylDuplex => 2,
            Silac::Triplex | Silac::DimethylTriplex => 3,
        }
    }
}

/// The light dimethyl label (K and N-terminus) Sage adds as a static
/// modification: not shown in the modifications table
fn is_light_dimethyl(m: &Modification) -> bool {
    m.kind == ModKind::Fixed
        && (m.target == "K" || m.target == "^")
        && (m.mass - f64::from(sage_core::silac::DIMETHYL_LIGHT)).abs() < 1e-3
}

/// SILAC label presets: (name, residue, mass)
pub const SILAC_LABELS: [(&str, char, f64); 5] = [
    ("Lys4 (2H4)", 'K', 4.025107),
    ("Lys6 (13C6)", 'K', 6.020129),
    ("Lys8 (13C6 15N2)", 'K', 8.014199),
    ("Arg6 (13C6)", 'R', 6.020129),
    ("Arg10 (13C6 15N4)", 'R', 10.008269),
];

pub const ION_KINDS: [&str; 7] = ["a", "b", "c", "x", "y", "z", "z_radical"];

/// A protease Sage knows by name (`sage_core::enzyme::protease`)
pub struct Protease {
    pub name: &'static str,
    pub label: &'static str,
    pub cleave_at: &'static str,
    pub restrict: &'static str,
    pub c_terminal: bool,
    pub help: &'static str,
}

pub const PROTEASES: &[Protease] = &[
    Protease { name: "trypsin", label: "Trypsin", cleave_at: "KR", restrict: "P", c_terminal: true,
        help: "Cleaves after K and R, not before P. The standard enzyme for bottom-up proteomics." },
    Protease { name: "trypsin/p", label: "Trypsin/P", cleave_at: "KR", restrict: "", c_terminal: true,
        help: "Cleaves after K and R, also before P." },
    Protease { name: "lys-c", label: "Lys-C", cleave_at: "K", restrict: "P", c_terminal: true,
        help: "Cleaves after K, not before P. Often combined with trypsin for more complete digestion." },
    Protease { name: "lys-n", label: "Lys-N", cleave_at: "K", restrict: "", c_terminal: false,
        help: "Cleaves before K (N-terminal side)." },
    Protease { name: "arg-c", label: "Arg-C", cleave_at: "R", restrict: "P", c_terminal: true,
        help: "Cleaves after R, not before P." },
    Protease { name: "asp-n", label: "Asp-N", cleave_at: "D", restrict: "", c_terminal: false,
        help: "Cleaves before D (N-terminal side). Combined with trypsin for sequential digests, e.g. in XL-MS." },
    Protease { name: "glu-c", label: "Glu-C", cleave_at: "E", restrict: "P", c_terminal: true,
        help: "Cleaves after E, not before P (ammonium bicarbonate buffer). In phosphate buffer it also \
               cuts after D: add a custom rule cleave at D." },
    Protease { name: "chymotrypsin", label: "Chymotrypsin", cleave_at: "FYWL", restrict: "P", c_terminal: true,
        help: "Cleaves after F, Y, W and L, not before P (low specificity)." },
    Protease { name: "chymotrypsin-high", label: "Chymotrypsin (high specificity)", cleave_at: "FYW", restrict: "P", c_terminal: true,
        help: "Cleaves after F, Y and W only, not before P." },
    Protease { name: "proalanase", label: "ProAlanase", cleave_at: "PA", restrict: "", c_terminal: true,
        help: "Cleaves after P and A. Useful for proline-rich proteins and histones." },
];

/// A hand-written cleavage rule
#[derive(Clone, Debug, PartialEq)]
pub struct EnzymeRule {
    pub cleave_at: String,
    pub restrict: String,
    pub c_terminal: bool,
}

/// Lower-case alphanumerics only, so `Asp-N` matches `aspn`
fn normalize(name: &str) -> String {
    name.chars()
        .filter(char::is_ascii_alphanumeric)
        .collect::<String>()
        .to_ascii_lowercase()
}

/// Sort a rule into a named protease when one matches it exactly
fn classify(rule: EnzymeRule, names: &mut Vec<String>, custom: &mut Vec<EnzymeRule>) {
    let preset = PROTEASES.iter().find(|p| {
        p.cleave_at == rule.cleave_at && p.restrict == rule.restrict && p.c_terminal == rule.c_terminal
    });
    match preset {
        Some(p) if !names.iter().any(|n| n == p.name) => names.push(p.name.into()),
        Some(_) => {}
        None => custom.push(rule),
    }
}

/// `(named proteases, custom rules, non-specific)` from a Sage `enzyme` block
fn parse_enzymes(enzyme: Option<&Map<String, Value>>) -> (Vec<String>, Vec<EnzymeRule>, bool) {
    let (mut names, mut custom) = (Vec::new(), Vec::new());
    let Some(e) = enzyme else {
        // Sage's default enzyme is trypsin
        return (vec!["trypsin".into()], custom, false);
    };
    if let Some(list) = e.get("enzymes").and_then(Value::as_array) {
        for item in list {
            match item {
                Value::String(name) => {
                    if let Some(p) = PROTEASES.iter().find(|p| normalize(p.name) == normalize(name)) {
                        if !names.iter().any(|n| n == p.name) {
                            names.push(p.name.into());
                        }
                    }
                }
                Value::Object(o) => classify(
                    EnzymeRule {
                        cleave_at: str_or(o.get("cleave_at"), ""),
                        restrict: str_or(o.get("restrict"), ""),
                        c_terminal: bool_or(o.get("c_terminal"), true),
                    },
                    &mut names,
                    &mut custom,
                ),
                _ => {}
            }
        }
        return (names, custom, false);
    }
    // Single enzyme, Sage's defaults for missing keys
    let rule = EnzymeRule {
        cleave_at: str_or(e.get("cleave_at"), "KR"),
        restrict: str_or(e.get("restrict"), ""),
        c_terminal: bool_or(e.get("c_terminal"), true),
    };
    if rule.cleave_at.is_empty() {
        return (names, custom, true);
    }
    classify(rule, &mut names, &mut custom);
    (names, custom, false)
}

#[derive(Clone, Debug)]
pub struct Settings {
    // Files
    pub fasta: String,
    /// Two-database crosslink search: crosslinks only between a protein of
    /// Database 1 and one of Database 2; set, they replace `fasta`
    pub database_1: String,
    pub database_2: String,
    pub spectra: Vec<String>,
    pub output_directory: String,

    // Database & enzyme
    pub missed_cleavages: u32,
    pub min_len: u32,
    pub max_len: u32,
    /// Named proteases used together (co-digestion), see [`PROTEASES`]
    pub enzymes: Vec<String>,
    /// Further hand-written cleavage rules
    pub custom_enzymes: Vec<EnzymeRule>,
    /// Non-specific search: every peptide bond
    pub nonspecific: bool,
    pub semi_enzymatic: bool,
    /// `None` lets Sage choose (500 Da, or 350 Da for crosslink searches)
    pub peptide_min_mass: Option<f64>,
    pub peptide_max_mass: f64,
    pub ion_kinds: [bool; 7],
    pub min_ion_index: u32,
    pub max_variable_mods: u32,
    pub max_combinations: Option<u32>,
    pub generate_decoys: bool,
    pub decoy_tag: String,
    pub bucket_size: u32,
    pub prefilter: bool,
    pub prefilter_chunk_size: u32,
    pub prefilter_low_memory: bool,
    /// Also search proteins without their initiator Met (MetAP)
    pub clip_n_term_met: bool,

    // Modifications
    pub mods: Vec<Modification>,
    /// Known modification sites for "library sites" modifications (empty: none)
    pub site_library: String,
    /// With a site_q column: only sites at or below it
    pub site_library_q: f64,

    // Spectra & tolerances
    pub precursor_tol: Tolerance,
    pub fragment_tol: Tolerance,
    pub precursor_charge: (u8, u8),
    pub override_precursor_charge: bool,
    pub isotope_errors: (i8, i8),
    pub deisotope: bool,
    pub chimera: bool,
    pub wide_window: bool,
    pub min_peaks: u32,
    /// `None` lets Sage choose (150, or 300 for cleavable crosslinkers)
    pub max_peaks: Option<u32>,
    pub max_fragment_charge: Option<u8>,
    pub min_matched_peaks: u32,
    pub report_psms: u32,

    // Scoring, quantification & output
    pub openms_score: bool,
    pub predict_rt: bool,
    pub protein_grouping: bool,
    pub protein_grouping_peptide_fdr: f64,
    pub write_pin: bool,
    pub write_report: bool,
    /// `results.annotator.csv` / `crosslinks.annotator.csv` for Annotator
    pub write_annotator: bool,
    /// Search MS2 spectra for diagnostic ions (diagnostic_ions.tsv)
    pub diagnostic_ions: bool,
    /// A custom ion list from a loaded configuration, written back as it was
    pub diagnostic_ions_list: Option<Value>,
    pub annotate_matches: bool,
    pub lfq: bool,
    /// MS1 mass tolerance (ppm) for LFQ
    pub lfq_ppm: f64,
    /// Only peptides / crosslinks at or below this q-value are quantified
    pub lfq_q: f64,
    pub lfq_combine_charges: bool,
    /// Write lfq_long.tsv (one row per precursor and run)
    pub lfq_long_form: bool,
    /// Minimum isotope-envelope spectral angle of an integrated peak
    pub lfq_spectral_angle: f64,
    pub tmt: Tmt,
    pub tmt_level: u8,
    pub tmt_sn: bool,
    pub silac: Silac,
    /// Label masses (K, R) of the medium (triplex) and heavy channels; 0: unlabelled
    pub silac_medium: (f64, f64),
    pub silac_heavy: (f64, f64),
    /// Crosslinks between peptides of different channels
    pub silac_mixed: bool,
    /// Proteins uniformly labelled with 15N and/or 13C (accession → label)
    pub labelled_proteins: std::collections::BTreeMap<String, sage_core::labelled::Isotope>,
    /// Also search the labelled proteins unlabelled
    pub include_unlabelled: bool,

    pub crosslink: Crosslink,

    /// timsTOF 1/K0 on timsrust's linear scale instead of the instrument's
    /// calibration (`bruker_config.ion_mobility_scale`)
    pub linear_mobility: bool,
    /// Align runs with a straight line (`retention_time_alignment: linear`)
    /// instead of the nonlinear warp
    pub linear_alignment: bool,
    /// Entrapment proteins: accessions containing this tag (empty = off)
    pub entrapment_tag: String,
    /// Localize variable modifications (site probabilities, FLR, sites.tsv)
    pub localization: bool,
    /// Residues used as decoy sites for the false localization rate
    pub localization_decoys: String,

    /// Keys from a loaded file the GUI does not edit, written back unchanged
    extra: Map<String, Value>,
    extra_database: Map<String, Value>,
    extra_enzyme: Map<String, Value>,
    extra_quant: Map<String, Value>,
    extra_lfq: Map<String, Value>,
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            fasta: String::new(),
            database_1: String::new(),
            database_2: String::new(),
            spectra: Vec::new(),
            output_directory: String::new(),
            missed_cleavages: 2,
            min_len: 5,
            max_len: 50,
            enzymes: vec!["trypsin".into()],
            custom_enzymes: Vec::new(),
            nonspecific: false,
            semi_enzymatic: false,
            peptide_min_mass: None,
            peptide_max_mass: 5000.0,
            ion_kinds: [false, true, false, false, true, false, false],
            min_ion_index: 2,
            max_variable_mods: 2,
            max_combinations: None,
            generate_decoys: true,
            decoy_tag: "rev_".into(),
            bucket_size: 8192,
            prefilter: false,
            prefilter_chunk_size: 0,
            prefilter_low_memory: true,
            clip_n_term_met: false,
            mods: vec![(&PRESETS[0]).into(), (&PRESETS[2]).into()],
            site_library: String::new(),
            site_library_q: 0.01,
            precursor_tol: Tolerance::ppm(-10.0, 10.0),
            fragment_tol: Tolerance::ppm(-10.0, 10.0),
            precursor_charge: (2, 4),
            override_precursor_charge: false,
            isotope_errors: (0, 0),
            deisotope: true,
            chimera: false,
            wide_window: false,
            min_peaks: 15,
            max_peaks: None,
            max_fragment_charge: None,
            min_matched_peaks: 4,
            report_psms: 1,
            openms_score: false,
            predict_rt: true,
            protein_grouping: true,
            protein_grouping_peptide_fdr: 0.01,
            write_pin: false,
            write_report: true,
            write_annotator: false,
            diagnostic_ions: false,
            diagnostic_ions_list: None,
            annotate_matches: false,
            lfq: false,
            lfq_ppm: 5.0,
            lfq_q: 0.01,
            lfq_combine_charges: true,
            lfq_long_form: true,
            lfq_spectral_angle: 0.7,
            tmt: Tmt::None,
            tmt_level: 3,
            tmt_sn: false,
            silac: Silac::Off,
            silac_medium: (4.025107, 6.020129),
            silac_heavy: (8.014199, 10.008269),
            silac_mixed: false,
            labelled_proteins: Default::default(),
            include_unlabelled: true,
            crosslink: Crosslink::default(),
            linear_mobility: false,
            linear_alignment: false,
            entrapment_tag: String::new(),
            localization: false,
            localization_decoys: "A".into(),
            extra: Map::new(),
            extra_database: Map::new(),
            extra_enzyme: Map::new(),
            extra_quant: Map::new(),
            extra_lfq: Map::new(),
        }
    }
}

const TOP_KEYS: &[&str] = &[
    "database", "precursor_tol", "fragment_tol", "precursor_charge", "override_precursor_charge",
    "isotope_errors", "deisotope", "chimera", "wide_window", "min_peaks", "max_peaks",
    "max_fragment_charge", "min_matched_peaks", "report_psms", "score_type", "predict_rt",
    "protein_grouping", "protein_grouping_peptide_fdr", "write_pin", "write_report",
    "write_annotator", "annotate_matches", "diagnostic_ions", "retention_time_alignment", "entrapment", "localization", "quant", "crosslink", "mzml_paths", "output_directory",
];
const DATABASE_KEYS: &[&str] = &[
    "fasta", "enzyme", "peptide_min_mass", "peptide_max_mass", "ion_kinds", "min_ion_index",
    "static_mods", "variable_mods", "max_variable_mods", "max_combinations", "generate_decoys",
    "decoy_tag", "bucket_size", "prefilter", "prefilter_chunk_size", "prefilter_low_memory", "silac",
    "labelled_proteins", "crosslink_databases", "clip_n_term_met", "site_library",
];
const ENZYME_KEYS: &[&str] = &[
    "missed_cleavages", "min_len", "max_len", "cleave_at", "restrict", "c_terminal", "enzymes",
    "semi_enzymatic",
];
const QUANT_KEYS: &[&str] = &["lfq", "tmt", "tmt_settings", "lfq_settings"];
const LFQ_KEYS: &[&str] = &["ppm_tolerance", "peptide_q_value", "combine_charge_states", "spectral_angle", "long_form"];

impl Settings {
    /// Two-database crosslink search: Database 1 or 2 set (they replace the
    /// FASTA file)
    pub fn two_databases(&self) -> bool {
        !self.database_1.trim().is_empty() || !self.database_2.trim().is_empty()
    }

    /// The FASTA files searched (set ones only): the FASTA file, or Database
    /// 1 and 2
    pub fn fasta_files(&self) -> Vec<&str> {
        let files = if self.two_databases() {
            vec![self.database_1.trim(), self.database_2.trim()]
        } else {
            vec![self.fasta.trim()]
        };
        files.into_iter().filter(|f| !f.is_empty()).collect()
    }
}

/// The proteins of a FASTA file: accession (the header's first word, as Sage
/// names it) and the rest of the header; `None` if it cannot be read
pub fn fasta_proteins(path: &str) -> Option<Vec<(String, String)>> {
    let text = std::fs::read_to_string(path).ok()?;
    Some(
        text.lines()
            .filter_map(|l| l.strip_prefix('>'))
            .map(|header| {
                let (accession, rest) = header.trim().split_once(char::is_whitespace).unwrap_or((header.trim(), ""));
                // UniProt: the name, up to `OS=`
                let name = rest[..rest.find(" OS=").unwrap_or(rest.len())].trim();
                (accession.to_string(), name.to_string())
            })
            .collect(),
    )
}

fn fasta_accessions(path: &str) -> Option<Vec<String>> {
    Some(fasta_proteins(path)?.into_iter().map(|(a, _)| a).collect())
}

/// The GUI's SILAC settings from `database.silac`: (mode, medium, heavy, mixed)
fn parse_silac(v: Option<&Value>, d: &Settings) -> (Silac, (f64, f64), (f64, f64), bool) {
    let Some(silac) = v.and_then(Value::as_object) else {
        return (Silac::Off, d.silac_medium, d.silac_heavy, false);
    };
    let labels = |c: &Value| {
        let l = c.get("labels");
        let mass = |r: &str| l.and_then(|l| l.get(r)).and_then(Value::as_f64).unwrap_or(0.0);
        (mass("K"), mass("R"))
    };
    let channels: Vec<&Value> = silac.get("channels").and_then(Value::as_array).map(|a| a.iter().collect()).unwrap_or_default();
    let mixed = silac.get("mixed_crosslinks").and_then(Value::as_bool).unwrap_or(false);
    // dimethyl: the light channel is labelled too
    if silac.get("light").and_then(Value::as_object).is_some_and(|l| !l.is_empty()) {
        let mode = if channels.len() >= 2 { Silac::DimethylTriplex } else { Silac::DimethylDuplex };
        return (mode, d.silac_medium, d.silac_heavy, mixed);
    }
    match channels[..] {
        [heavy] => (Silac::Duplex, d.silac_medium, labels(heavy), mixed),
        [medium, heavy, ..] => (Silac::Triplex, labels(medium), labels(heavy), mixed),
        [] => (Silac::Off, d.silac_medium, d.silac_heavy, false),
    }
}

fn f64_or(v: Option<&Value>, default: f64) -> f64 {
    v.and_then(Value::as_f64).unwrap_or(default)
}

fn u32_or(v: Option<&Value>, default: u32) -> u32 {
    v.and_then(Value::as_u64).map(|x| x as u32).unwrap_or(default)
}

fn bool_or(v: Option<&Value>, default: bool) -> bool {
    v.and_then(Value::as_bool).unwrap_or(default)
}

fn str_or(v: Option<&Value>, default: &str) -> String {
    v.and_then(Value::as_str).unwrap_or(default).to_string()
}

fn pair<T: TryFrom<i64>>(v: Option<&Value>, default: (T, T)) -> (T, T) {
    let get = |i: usize| {
        v.and_then(|v| v.get(i))
            .and_then(Value::as_i64)
            .and_then(|x| T::try_from(x).ok())
    };
    match (get(0), get(1)) {
        (Some(a), Some(b)) => (a, b),
        _ => default,
    }
}

/// Keys of `obj` that are not in `known`
fn unknown_keys(obj: Option<&Map<String, Value>>, known: &[&str]) -> Map<String, Value> {
    obj.map(|o| {
        o.iter()
            .filter(|(k, _)| !known.contains(&k.as_str()))
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect()
    })
    .unwrap_or_default()
}

/// Warning shown while the multivalent crosslink search is on
pub fn multivalent_warning(valency: u8) -> String {
    format!(
        "Multivalent crosslinks (up to {valency} linkers per species): the search space is very large. \
         Every combination of up to {} peptides and every arrangement of their links is scored, so \
         searches take far longer and use more memory. Only recommended with a small protein subset \
         (a purified complex of a few dozen proteins), not a proteome.",
        valency + 1
    )
}

/// Parse a comma separated list of numbers
pub fn parse_list(s: &str) -> Result<Vec<f64>, String> {
    s.split(',')
        .map(str::trim)
        .filter(|x| !x.is_empty())
        .map(|x| x.parse::<f64>().map_err(|_| format!("`{x}` is not a number")))
        .collect()
}

fn mod_entry(m: &Modification) -> Value {
    let mut e = Map::new();
    e.insert("mass".into(), json!(m.mass));
    if m.kind == ModKind::Variable && m.max_count > 0 {
        e.insert("max_count".into(), json!(m.max_count));
    }
    if !m.name.trim().is_empty() {
        e.insert("name".into(), json!(m.name.trim()));
    }
    if m.position != SitePosition::Anywhere {
        e.insert("position".into(), serde_json::to_value(m.position).unwrap_or_default());
    }
    if !m.motif.trim().is_empty() {
        e.insert("motif".into(), json!(m.motif.trim()));
    }
    if let Ok(losses) = parse_losses(&m.losses) {
        if !losses.is_empty() {
            e.insert("neutral_losses".into(), json!(losses));
        }
    }
    if m.kind == ModKind::Variable {
        match m.search {
            ModSearch::Indexed => {}
            ModSearch::Library => {
                e.insert("site_mode".into(), json!("library"));
            }
            ModSearch::Offset => {
                e.insert("search_mode".into(), json!("mass_offset"));
            }
        }
    }
    if e.len() == 1 {
        json!(m.mass)
    } else {
        Value::Object(e)
    }
}

fn parse_mod_entry(target: &str, v: &Value, kind: ModKind) -> Option<Modification> {
    let (mass, name, max_count) = match v {
        Value::Number(n) => (n.as_f64()?, String::new(), 0),
        Value::Object(o) => (
            o.get("mass")?.as_f64()?,
            str_or(o.get("name"), ""),
            u32_or(o.get("max_count"), 0),
        ),
        _ => return None,
    };
    let field = |k: &str| v.as_object().and_then(|o| o.get(k));
    Some(Modification {
        target: target.into(),
        mass,
        name,
        kind,
        max_count,
        position: field("position").and_then(|p| serde_json::from_value(p.clone()).ok()).unwrap_or_default(),
        motif: str_or(field("motif"), ""),
        losses: field("neutral_losses")
            .and_then(Value::as_array)
            .map(|l| l.iter().filter_map(Value::as_f64).map(|x| x.to_string()).collect::<Vec<_>>().join(", "))
            .unwrap_or_default(),
        search: if field("search_mode").and_then(Value::as_str) == Some("mass_offset") {
            ModSearch::Offset
        } else if field("site_mode").and_then(Value::as_str) == Some("library") {
            ModSearch::Library
        } else {
            ModSearch::Indexed
        },
    })
}

/// A configuration key's modifications: named and rule keys
/// (`internal_residue:K`, `motif:N-{P}-[ST]`) become a target with the rule
fn mods_of_key(key: &str, v: &Value, kind: ModKind) -> Vec<Modification> {
    let Some(m) = parse_mod_entry(key, v, kind) else { return Vec::new() };
    match sage_core::site_rules::parse_site_key(key) {
        Ok(sites) => sites
            .into_iter()
            .map(|(specificity, rule)| {
                let mut m = m.clone();
                m.target = specificity.to_string();
                if !rule.position.is_anywhere() {
                    m.position = rule.position;
                }
                if let Some(motif) = rule.motif {
                    m.motif = motif.pattern().to_string();
                }
                m
            })
            .collect(),
        // kept as written; check() reports it
        Err(_) => vec![m],
    }
}

impl Settings {
    pub fn is_thermo_raw(path: &str) -> bool {
        path.to_ascii_lowercase().ends_with(".raw") && std::path::Path::new(path).is_file()
    }

    /// Build the Sage JSON configuration. `spectra` replaces the file list
    /// (used to substitute converted mzML files for Thermo .raw files).
    pub fn to_json(&self, spectra: &[String]) -> Value {
        let mut enzyme = self.extra_enzyme.clone();
        enzyme.insert("missed_cleavages".into(), json!(self.missed_cleavages));
        enzyme.insert("min_len".into(), json!(self.min_len));
        enzyme.insert("max_len".into(), json!(self.max_len));
        match (self.nonspecific, &self.enzymes[..], &self.custom_enzymes[..]) {
            (true, _, _) => {
                enzyme.insert("cleave_at".into(), json!(""));
            }
            // A single hand-written rule keeps the classic form, which also
            // allows `$` (no cleavage)
            (false, [], [rule]) => {
                enzyme.insert("cleave_at".into(), json!(rule.cleave_at));
                enzyme.insert("restrict".into(), json!(rule.restrict));
                enzyme.insert("c_terminal".into(), json!(rule.c_terminal));
            }
            (false, names, rules) => {
                let list = names
                    .iter()
                    .map(|n| json!(n))
                    .chain(rules.iter().map(|r| {
                        json!({
                            "cleave_at": r.cleave_at,
                            "restrict": r.restrict,
                            "c_terminal": r.c_terminal,
                        })
                    }))
                    .collect::<Vec<_>>();
                enzyme.insert("enzymes".into(), Value::Array(list));
            }
        }
        enzyme.insert("semi_enzymatic".into(), json!(self.semi_enzymatic));

        let mut static_mods = Map::new();
        let mut variable_mods = Map::new();
        for m in &self.mods {
            match m.kind {
                ModKind::Fixed => {
                    static_mods.insert(m.target.clone(), mod_entry(m));
                }
                ModKind::Variable => {
                    let list = variable_mods
                        .entry(m.target.clone())
                        .or_insert_with(|| json!([]));
                    list.as_array_mut().unwrap().push(mod_entry(m));
                }
            }
        }

        let ion_kinds: Vec<&str> = ION_KINDS
            .iter()
            .zip(self.ion_kinds)
            .filter(|(_, on)| *on)
            .map(|(k, _)| *k)
            .collect();

        let mut db = self.extra_database.clone();
        // two crosslink databases replace the FASTA file
        if self.two_databases() {
            db.insert(
                "crosslink_databases".into(),
                json!({"database_1": self.database_1.trim(), "database_2": self.database_2.trim()}),
            );
        } else {
            db.insert("fasta".into(), json!(self.fasta));
        }
        db.insert("enzyme".into(), Value::Object(enzyme));
        if let Some(mass) = self.peptide_min_mass {
            db.insert("peptide_min_mass".into(), json!(mass));
        }
        db.insert("peptide_max_mass".into(), json!(self.peptide_max_mass));
        db.insert("ion_kinds".into(), json!(ion_kinds));
        db.insert("min_ion_index".into(), json!(self.min_ion_index));
        db.insert("static_mods".into(), Value::Object(static_mods));
        db.insert("variable_mods".into(), Value::Object(variable_mods));
        db.insert("max_variable_mods".into(), json!(self.max_variable_mods));
        if let Some(n) = self.max_combinations {
            db.insert("max_combinations".into(), json!(n));
        }
        db.insert("generate_decoys".into(), json!(self.generate_decoys));
        // only when on (off is Sage's default), so configs still run in an older sage.exe
        if self.clip_n_term_met {
            db.insert("clip_n_term_met".into(), json!(true));
        }
        if !self.site_library.trim().is_empty() {
            db.insert(
                "site_library".into(),
                json!({"path": self.site_library.trim(), "max_q": self.site_library_q}),
            );
        }
        db.insert("decoy_tag".into(), json!(self.decoy_tag));
        db.insert("bucket_size".into(), json!(self.bucket_size));
        db.insert("prefilter".into(), json!(self.prefilter));
        if self.prefilter {
            db.insert("prefilter_chunk_size".into(), json!(self.prefilter_chunk_size));
            db.insert("prefilter_low_memory".into(), json!(self.prefilter_low_memory));
        }
        if let Some(silac) = self.silac_json() {
            db.insert("silac".into(), silac);
        }
        if !self.labelled_proteins.is_empty() {
            let labelled = sage_core::labelled::LabelledProteins {
                proteins: self.labelled_proteins.clone(),
                include_unlabelled: self.include_unlabelled,
            };
            db.insert("labelled_proteins".into(), serde_json::to_value(labelled).unwrap_or_default());
        }

        let mut quant = self.extra_quant.clone();
        quant.insert("lfq".into(), json!(self.lfq));
        let mut lfq = self.extra_lfq.clone();
        lfq.insert("ppm_tolerance".into(), json!(self.lfq_ppm));
        lfq.insert("peptide_q_value".into(), json!(self.lfq_q));
        lfq.insert("combine_charge_states".into(), json!(self.lfq_combine_charges));
        lfq.insert("spectral_angle".into(), json!(self.lfq_spectral_angle));
        // on by default in sage.exe: written only when off
        if !self.lfq_long_form {
            lfq.insert("long_form".into(), json!(false));
        }
        quant.insert("lfq_settings".into(), Value::Object(lfq));
        if let Some(tmt) = self.tmt.key() {
            quant.insert("tmt".into(), json!(tmt));
            quant.insert(
                "tmt_settings".into(),
                json!({ "level": self.tmt_level, "sn": self.tmt_sn }),
            );
        }

        let mut top = self.extra.clone();
        // the one key of bruker_config the GUI edits; the rest stays as loaded
        let bruker = top.entry("bruker_config").or_insert_with(|| json!({}));
        if let Some(bruker) = bruker.as_object_mut() {
            if self.linear_mobility {
                bruker.insert("ion_mobility_scale".into(), json!("linear"));
            } else {
                bruker.remove("ion_mobility_scale");
            }
        }
        if top.get("bruker_config").and_then(Value::as_object).is_some_and(|b| b.is_empty()) {
            top.remove("bruker_config");
        }
        top.insert("database".into(), Value::Object(db));
        top.insert("precursor_tol".into(), self.precursor_tol.to_json());
        top.insert("fragment_tol".into(), self.fragment_tol.to_json());
        top.insert(
            "precursor_charge".into(),
            json!([self.precursor_charge.0, self.precursor_charge.1]),
        );
        top.insert(
            "override_precursor_charge".into(),
            json!(self.override_precursor_charge),
        );
        top.insert(
            "isotope_errors".into(),
            json!([self.isotope_errors.0, self.isotope_errors.1]),
        );
        top.insert("deisotope".into(), json!(self.deisotope));
        top.insert("chimera".into(), json!(self.chimera));
        top.insert("wide_window".into(), json!(self.wide_window));
        top.insert("min_peaks".into(), json!(self.min_peaks));
        if let Some(n) = self.max_peaks {
            top.insert("max_peaks".into(), json!(n));
        }
        if let Some(z) = self.max_fragment_charge {
            top.insert("max_fragment_charge".into(), json!(z));
        }
        top.insert("min_matched_peaks".into(), json!(self.min_matched_peaks));
        top.insert("report_psms".into(), json!(self.report_psms));
        top.insert(
            "score_type".into(),
            json!(if self.openms_score { "OpenMSHyperScore" } else { "SageHyperScore" }),
        );
        top.insert("predict_rt".into(), json!(self.predict_rt));
        top.insert("protein_grouping".into(), json!(self.protein_grouping));
        top.insert(
            "protein_grouping_peptide_fdr".into(),
            json!(self.protein_grouping_peptide_fdr),
        );
        top.insert("write_pin".into(), json!(self.write_pin));
        top.insert("write_report".into(), json!(self.write_report));
        top.insert("write_annotator".into(), json!(self.write_annotator));
        if self.localization {
            let decoys: String = self.localization_decoys.chars().filter(|c| c.is_ascii_alphabetic()).collect();
            let decoys = if decoys.is_empty() { "A".to_string() } else { decoys.to_ascii_uppercase() };
            top.insert("localization".into(), json!({"decoy_residues": decoys}));
        }
        if !self.entrapment_tag.trim().is_empty() {
            top.insert("entrapment".into(), json!({"tag": self.entrapment_tag.trim()}));
        }
        // only when not the default (nonlinear)
        if self.linear_alignment {
            top.insert("retention_time_alignment".into(), json!("linear"));
        }
        // only when on (off is Sage's default)
        if self.diagnostic_ions {
            top.insert("diagnostic_ions".into(), self.diagnostic_ions_list.clone().unwrap_or(json!(true)));
        }
        top.insert("annotate_matches".into(), json!(self.annotate_matches));
        top.insert("quant".into(), Value::Object(quant));
        if self.crosslink.enabled {
            top.insert("crosslink".into(), self.crosslink.to_json());
        }
        top.insert("mzml_paths".into(), json!(spectra));
        top.insert("output_directory".into(), json!(self.output_directory));
        Value::Object(top)
    }

    /// Read a Sage JSON configuration. Missing keys take Sage's defaults.
    pub fn from_json(v: &Value) -> Result<Self, String> {
        let top = v.as_object().ok_or("the file is not a JSON object")?;
        let db = top.get("database").and_then(Value::as_object);
        let enzyme = db.and_then(|d| d.get("enzyme")).and_then(Value::as_object);
        let quant = top.get("quant").and_then(Value::as_object);
        let d = Settings::default();
        let dbv = |k: &str| db.and_then(|d| d.get(k));
        let ev = |k: &str| enzyme.and_then(|e| e.get(k));
        let lfq_settings = quant.and_then(|q| q.get("lfq_settings")).and_then(Value::as_object);
        let lv = |k: &str| lfq_settings.and_then(|l| l.get(k));
        let silac = parse_silac(dbv("silac"), &d);
        let labelled: Option<sage_core::labelled::LabelledProteins> =
            dbv("labelled_proteins").and_then(|v| serde_json::from_value(v.clone()).ok());

        let mut mods = Vec::new();
        if let Some(fixed) = dbv("static_mods").and_then(Value::as_object) {
            for (target, entry) in fixed {
                mods.extend(mods_of_key(target, entry, ModKind::Fixed));
            }
        }
        // Sage adds dimethyl's light label as a fixed modification (a search's
        // own configuration lists it): it is part of the dimethyl setting
        if silac.0.is_dimethyl() {
            mods.retain(|m| !is_light_dimethyl(m));
        }
        if let Some(var) = dbv("variable_mods").and_then(Value::as_object) {
            for (target, entries) in var {
                match entries {
                    Value::Array(list) => mods.extend(
                        list.iter()
                            .flat_map(|e| mods_of_key(target, e, ModKind::Variable)),
                    ),
                    e => mods.extend(mods_of_key(target, e, ModKind::Variable)),
                }
            }
        }

        let mut ion_kinds = d.ion_kinds;
        if let Some(kinds) = dbv("ion_kinds").and_then(Value::as_array) {
            ion_kinds = [false; 7];
            for k in kinds.iter().filter_map(Value::as_str) {
                if let Some(i) = ION_KINDS.iter().position(|x| *x == k) {
                    ion_kinds[i] = true;
                }
            }
        }

        let (enzymes, custom_enzymes, nonspecific) = parse_enzymes(enzyme);

        let tmt = quant
            .and_then(|q| q.get("tmt"))
            .and_then(Value::as_str)
            .and_then(|t| Tmt::ALL.iter().find(|x| x.key() == Some(t)).copied())
            .unwrap_or(Tmt::None);
        let tmt_settings = quant
            .and_then(|q| q.get("tmt_settings"))
            .and_then(Value::as_object);

        Ok(Settings {
            fasta: str_or(dbv("fasta"), ""),
            database_1: str_or(dbv("crosslink_databases").and_then(|d| d.get("database_1")), ""),
            database_2: str_or(dbv("crosslink_databases").and_then(|d| d.get("database_2")), ""),
            spectra: top
                .get("mzml_paths")
                .and_then(Value::as_array)
                .map(|a| a.iter().filter_map(Value::as_str).map(String::from).collect())
                .unwrap_or_default(),
            output_directory: str_or(top.get("output_directory"), ""),
            missed_cleavages: u32_or(ev("missed_cleavages"), 0),
            min_len: u32_or(ev("min_len"), d.min_len),
            max_len: u32_or(ev("max_len"), d.max_len),
            enzymes,
            custom_enzymes,
            nonspecific,
            semi_enzymatic: bool_or(ev("semi_enzymatic"), d.semi_enzymatic),
            peptide_min_mass: dbv("peptide_min_mass").and_then(Value::as_f64),
            peptide_max_mass: f64_or(dbv("peptide_max_mass"), d.peptide_max_mass),
            ion_kinds,
            min_ion_index: u32_or(dbv("min_ion_index"), d.min_ion_index),
            max_variable_mods: u32_or(dbv("max_variable_mods"), d.max_variable_mods),
            max_combinations: dbv("max_combinations").and_then(Value::as_u64).map(|x| x as u32),
            generate_decoys: bool_or(dbv("generate_decoys"), d.generate_decoys),
            decoy_tag: str_or(dbv("decoy_tag"), &d.decoy_tag),
            bucket_size: u32_or(dbv("bucket_size"), d.bucket_size),
            prefilter: bool_or(dbv("prefilter"), d.prefilter),
            prefilter_chunk_size: u32_or(dbv("prefilter_chunk_size"), d.prefilter_chunk_size),
            prefilter_low_memory: bool_or(dbv("prefilter_low_memory"), d.prefilter_low_memory),
            clip_n_term_met: bool_or(dbv("clip_n_term_met"), d.clip_n_term_met),
            mods,
            site_library: str_or(dbv("site_library").and_then(|l| l.get("path")), ""),
            site_library_q: dbv("site_library")
                .and_then(|l| l.get("max_q"))
                .and_then(Value::as_f64)
                .unwrap_or(d.site_library_q),
            precursor_tol: top
                .get("precursor_tol")
                .and_then(Tolerance::from_json)
                .unwrap_or(d.precursor_tol),
            fragment_tol: top
                .get("fragment_tol")
                .and_then(Tolerance::from_json)
                .unwrap_or(d.fragment_tol),
            precursor_charge: pair(top.get("precursor_charge"), d.precursor_charge),
            override_precursor_charge: bool_or(
                top.get("override_precursor_charge"),
                d.override_precursor_charge,
            ),
            isotope_errors: pair(top.get("isotope_errors"), d.isotope_errors),
            deisotope: bool_or(top.get("deisotope"), d.deisotope),
            chimera: bool_or(top.get("chimera"), d.chimera),
            wide_window: bool_or(top.get("wide_window"), d.wide_window),
            min_peaks: u32_or(top.get("min_peaks"), d.min_peaks),
            max_peaks: top.get("max_peaks").and_then(Value::as_u64).map(|x| x as u32),
            max_fragment_charge: top
                .get("max_fragment_charge")
                .and_then(Value::as_u64)
                .map(|x| x as u8),
            min_matched_peaks: u32_or(top.get("min_matched_peaks"), d.min_matched_peaks),
            report_psms: u32_or(top.get("report_psms"), d.report_psms),
            openms_score: top.get("score_type").and_then(Value::as_str)
                == Some("OpenMSHyperScore"),
            predict_rt: bool_or(top.get("predict_rt"), d.predict_rt),
            protein_grouping: bool_or(top.get("protein_grouping"), d.protein_grouping),
            protein_grouping_peptide_fdr: f64_or(
                top.get("protein_grouping_peptide_fdr"),
                d.protein_grouping_peptide_fdr,
            ),
            write_pin: bool_or(top.get("write_pin"), d.write_pin),
            write_report: bool_or(top.get("write_report"), d.write_report),
            // absent from a file means off, as in Sage
            write_annotator: bool_or(top.get("write_annotator"), false),
            diagnostic_ions: match top.get("diagnostic_ions") {
                Some(Value::Bool(on)) => *on,
                Some(Value::Array(_)) => true,
                _ => false,
            },
            diagnostic_ions_list: top.get("diagnostic_ions").filter(|v| v.is_array()).cloned(),
            annotate_matches: bool_or(top.get("annotate_matches"), d.annotate_matches),
            lfq: bool_or(quant.and_then(|q| q.get("lfq")), false),
            lfq_ppm: f64_or(lv("ppm_tolerance"), d.lfq_ppm),
            lfq_q: f64_or(lv("peptide_q_value"), d.lfq_q),
            lfq_combine_charges: bool_or(lv("combine_charge_states"), d.lfq_combine_charges),
            lfq_long_form: bool_or(lv("long_form"), d.lfq_long_form),
            lfq_spectral_angle: f64_or(lv("spectral_angle"), d.lfq_spectral_angle),
            tmt,
            tmt_level: tmt_settings
                .and_then(|t| t.get("level"))
                .and_then(Value::as_u64)
                .map(|x| x as u8)
                .unwrap_or(d.tmt_level),
            tmt_sn: bool_or(tmt_settings.and_then(|t| t.get("sn")), d.tmt_sn),
            silac: silac.0,
            silac_medium: silac.1,
            silac_heavy: silac.2,
            silac_mixed: silac.3,
            labelled_proteins: labelled.as_ref().map(|l| l.proteins.clone()).unwrap_or_default(),
            include_unlabelled: labelled.as_ref().is_none_or(|l| l.include_unlabelled),
            crosslink: top
                .get("crosslink")
                .and_then(Value::as_object)
                .map(Crosslink::from_json)
                .unwrap_or_default(),
            linear_alignment: top.get("retention_time_alignment").and_then(Value::as_str) == Some("linear"),
            localization: top.get("localization").is_some_and(|l| l.is_object()),
            localization_decoys: top
                .get("localization")
                .and_then(|l| l.get("decoy_residues"))
                .and_then(Value::as_str)
                .unwrap_or("A")
                .to_string(),
            entrapment_tag: top
                .get("entrapment")
                .and_then(|e| e.get("tag"))
                .and_then(Value::as_str)
                .unwrap_or("")
                .to_string(),
            linear_mobility: top
                .get("bruker_config")
                .and_then(|b| b.get("ion_mobility_scale"))
                .and_then(Value::as_str)
                == Some("linear"),
            extra: unknown_keys(Some(top), TOP_KEYS),
            extra_database: unknown_keys(db, DATABASE_KEYS),
            extra_enzyme: unknown_keys(enzyme, ENZYME_KEYS),
            extra_quant: unknown_keys(quant, QUANT_KEYS),
            extra_lfq: unknown_keys(lfq_settings, LFQ_KEYS),
        })
    }

    /// `database.silac`: the labelled channels, `None` when off
    pub fn silac_json(&self) -> Option<Value> {
        let channel = |name: &str, (k, r): (f64, f64)| {
            let mut labels = Map::new();
            for (residue, mass) in [("K", k), ("R", r)] {
                if mass > 0.0 {
                    labels.insert(residue.into(), json!(mass));
                }
            }
            json!({ "name": name, "labels": labels })
        };
        let channels = match self.silac {
            Silac::Off => return None,
            Silac::Duplex => vec![channel("heavy", self.silac_heavy)],
            Silac::Triplex => vec![channel("medium", self.silac_medium), channel("heavy", self.silac_heavy)],
            Silac::DimethylDuplex | Silac::DimethylTriplex => {
                // sage_core::silac::SilacSettings::dimethyl, with the masses as written
                let sites = |mass: f64| json!({ "K": mass, "^": mass });
                let mut channels = vec![json!({ "name": "medium", "labels": sites(32.056407) })];
                if self.silac == Silac::DimethylTriplex {
                    channels.push(json!({ "name": "heavy", "labels": sites(36.07567) }));
                }
                return Some(json!({
                    "light": sites(28.0313),
                    "channels": channels,
                    "mixed_crosslinks": self.silac_mixed,
                }));
            }
        };
        Some(json!({ "channels": channels, "mixed_crosslinks": self.silac_mixed }))
    }

    /// Problems that stop a search (`errors`) and ones worth a look (`warnings`)
    pub fn check(&self) -> (Vec<String>, Vec<String>) {
        let mut errors = Vec::new();
        let mut warnings = Vec::new();
        if !self.labelled_proteins.is_empty() {
            // every labelled protein should be in the FASTA
            let headers: Option<Vec<String>> = self
                .fasta_files()
                .into_iter()
                .map(fasta_accessions)
                .collect::<Option<Vec<_>>>()
                .map(|files| files.concat());
            if let Some(accessions) = headers {
                let missing: Vec<&String> = self
                    .labelled_proteins
                    .keys()
                    .filter(|k| !accessions.iter().any(|a| a == *k || a.split('|').skip(1).any(|p| p == k.as_str())))
                    .collect();
                if !missing.is_empty() {
                    errors.push(format!(
                        "Isotope-labelled proteins not in the FASTA: {}",
                        missing.iter().map(|s| s.as_str()).collect::<Vec<_>>().join(", ")
                    ));
                }
            }
            warnings.push(format!(
                "{} protein(s) are searched isotope-labelled{}. Incompletely labelled proteins have a lighter \
                 isotope envelope: allow isotope errors (e.g. -1 to 0) if the instrument picks a lighter peak.",
                self.labelled_proteins.len(),
                if self.include_unlabelled { " and unlabelled" } else { " only (not unlabelled)" }
            ));
        }
        if self.silac.is_dimethyl() {
            let clash: Vec<&str> = self
                .mods
                .iter()
                .filter(|m| m.kind == ModKind::Fixed && (m.target == "K" || m.target == "^") && !is_light_dimethyl(m))
                .map(|m| m.target.as_str())
                .collect();
            if !clash.is_empty() {
                errors.push(format!(
                    "Dimethyl labelling puts the light label (+28.0313) on K and the peptide N-terminus as a fixed \
                     modification: remove the other fixed modification on {}.",
                    clash.join(" and ")
                ));
            }
            warnings.push(format!(
                "Dimethyl: every peptide is searched in {} label states (light +28.0313, medium +32.0564{} on K and \
                 the peptide N-terminus), so the search takes about {} times as long. Ratios to light are written to \
                 silac.sage.tsv and silac_proteins.sage.tsv; they need MS1 spectra.",
                self.silac.states(),
                if self.silac == Silac::DimethylTriplex { ", heavy +36.0757" } else { "" },
                self.silac.states()
            ));
            if self.crosslink.enabled {
                warnings.push(
                    "Dimethyl with crosslinking: every lysine is searched dimethylated, so crosslinks through \
                     lysine are not found. Dimethyl labelling suits linkers on other residues."
                        .into(),
                );
            }
        }
        if self.silac.is_residue_silac() {
            let mut channels = vec![("heavy", self.silac_heavy)];
            if self.silac == Silac::Triplex {
                channels.push(("medium", self.silac_medium));
                if self.silac_medium == self.silac_heavy {
                    errors.push("SILAC: the medium and heavy channels have the same labels.".into());
                }
            }
            for (name, (k, r)) in channels {
                if k <= 0.0 && r <= 0.0 {
                    errors.push(format!("SILAC: the {name} channel labels neither K nor R."));
                }
            }
            let n = self.silac.states();
            warnings.push(format!(
                "SILAC: every peptide is searched in {n} label states, so the search takes about {n} times as long. \
                 Ratios to light are written to silac.sage.tsv and silac_proteins.sage.tsv{}; they need MS1 spectra.",
                if self.crosslink.enabled { " (crosslinks: crosslinks.silac.tsv)" } else { "" }
            ));
        }

        if self.two_databases() {
            for (name, path) in [("Database 1", &self.database_1), ("Database 2", &self.database_2)] {
                if path.trim().is_empty() {
                    errors.push(format!("Two-database crosslink search: select {name} (Files tab), or clear both."));
                } else if !std::path::Path::new(path.trim()).is_file() {
                    errors.push(format!("{name} not found: {}", path.trim()));
                }
            }
            if !self.crosslink.enabled {
                warnings.push(
                    "Database 1 and Database 2 are set but crosslinking is off: they are searched as one FASTA."
                        .into(),
                );
            } else {
                warnings.push(
                    "Two-database crosslink search: only crosslinks between a Database 1 protein and a Database 2 \
                     protein are searched (loop-links and crosslinks within one database are not). Linear \
                     peptides are searched in both."
                        .into(),
                );
            }
            if !self.database_1.trim().is_empty() && self.database_1.trim() == self.database_2.trim() {
                warnings.push(
                    "Database 1 and Database 2 are the same file: every protein is in both, so every crosslink \
                     is searched."
                        .into(),
                );
            }
        } else if self.fasta.trim().is_empty() {
            errors.push("Select a FASTA database (Files tab).".into());
        } else if !std::path::Path::new(&self.fasta).is_file() {
            errors.push(format!("FASTA file not found: {}", self.fasta));
        }
        if self.spectra.is_empty() {
            errors.push("Add at least one spectrum file or Bruker .d folder (Files tab).".into());
        }
        for f in &self.spectra {
            if !std::path::Path::new(f).exists() {
                errors.push(format!("Spectrum file not found: {f}"));
            }
        }
        if self.output_directory.trim().is_empty() {
            errors.push("Select an output folder for the results (Files tab).".into());
        }

        if self.min_len > self.max_len {
            errors.push("Minimum peptide length is larger than the maximum.".into());
        }
        if !self.nonspecific {
            if self.enzymes.is_empty() && self.custom_enzymes.is_empty() {
                errors.push("Select at least one enzyme, or choose a non-specific search (Database tab).".into());
            }
            let single = self.enzymes.is_empty() && self.custom_enzymes.len() == 1;
            for rule in &self.custom_enzymes {
                let no_cleavage = single && rule.cleave_at == "$";
                if !no_cleavage
                    && (rule.cleave_at.is_empty() || !rule.cleave_at.chars().all(|c| VALID_AA.contains(c)))
                {
                    errors.push(format!(
                        "Custom enzyme: `cleave at` must be amino acid letters from {VALID_AA} \
                         (or `$` for no cleavage, as the only enzyme), got `{}`.",
                        rule.cleave_at
                    ));
                }
                if !rule.restrict.chars().all(|c| VALID_AA.contains(c)) {
                    errors.push(format!("Custom enzyme: `restrict` may only contain {VALID_AA}."));
                }
            }
        }
        if !self.ion_kinds.iter().any(|x| *x) {
            errors.push("Select at least one fragment ion type.".into());
        }

        let variable = |search: ModSearch| -> Vec<&Modification> {
            self.mods.iter().filter(|m| m.kind == ModKind::Variable && m.search == search).collect()
        };
        let (offsets, library_mods) = (variable(ModSearch::Offset), variable(ModSearch::Library));
        let library = self.site_library.trim();
        if !library_mods.is_empty() && library.is_empty() {
            errors.push("Modifications searched at library sites need a site library (Modifications tab).".into());
        }
        if !library.is_empty() {
            if sage_cloudpath::try_parse_url(library).is_none_or(|u| u.scheme() == "file")
                && !std::path::Path::new(library).is_file()
            {
                errors.push(format!("Site library not found: {library}"));
            }
            if library_mods.is_empty() {
                warnings.push(
                    "A site library is set, but no variable modification is searched at library sites: it is not used."
                        .into(),
                );
            }
        }
        if !library_mods.is_empty() && !self.generate_decoys {
            errors.push(
                "Library sites need generated decoys (Database & enzyme tab): decoys must carry their targets' sites."
                    .into(),
            );
        }
        if !offsets.is_empty() {
            if self.crosslink.enabled {
                errors.push("Mass offsets are for linear searches: turn crosslinking off, or search the modification indexed.".into());
            }
            if self.silac != Silac::Off || !self.labelled_proteins.is_empty() {
                errors.push("Mass offsets cannot be combined with SILAC, dimethyl or labelled proteins.".into());
            }
            for o in &offsets {
                let indexed = self.mods.iter().any(|m| {
                    m.kind == ModKind::Variable
                        && m.search != ModSearch::Offset
                        && (m.mass - o.mass).abs() < 1e-6
                        && m.name.trim() == o.name.trim()
                });
                if indexed {
                    errors.push(format!(
                        "{} ({:+}) is both indexed and a mass offset: search it one way.",
                        if o.name.trim().is_empty() { "A modification" } else { o.name.trim() },
                        o.mass
                    ));
                }
            }
        }

        let mut fixed_targets = std::collections::HashSet::new();
        for m in &self.mods {
            if let Err(e) = check_target(&m.target) {
                errors.push(format!("Modification target {e}"));
            }
            if let Some(e) = rule_problem(m) {
                errors.push(format!("Modification {e}"));
            }
            if m.kind == ModKind::Fixed && !fixed_targets.insert(m.target.clone()) {
                errors.push(format!(
                    "Two fixed modifications on `{}`: only one fixed modification per target is allowed.",
                    m.target
                ));
            }
        }

        for (name, tol) in [("Precursor", &self.precursor_tol), ("Fragment", &self.fragment_tol)] {
            if tol.lo > tol.hi {
                errors.push(format!("{name} tolerance: the lower bound is above the upper bound."));
            }
        }
        if self.precursor_charge.0 > self.precursor_charge.1 {
            errors.push("Precursor charge: minimum is above maximum.".into());
        }
        if self.isotope_errors.0 > self.isotope_errors.1 {
            errors.push("Isotope errors: minimum is above maximum.".into());
        }

        if self.crosslink.enabled {
            let xl = &self.crosslink;
            if self.prefilter {
                errors.push("Database prefilter cannot be used with crosslink searches.".into());
            }
            if self.wide_window {
                errors.push("Wide-window (DIA) mode is not supported for crosslink searches.".into());
            }
            if self.missed_cleavages == 0 {
                warnings.push(
                    "Crosslink search with 0 missed cleavages finds no linkable peptides; use 2.".into(),
                );
            }
            if xl.multivalent {
                warnings.push(multivalent_warning(xl.valency));
            }
            // A fixed modification on a reactive residue (typically TMT on K) blocks every link there
            let reactive = if xl.two_sites {
                format!("{}{}", xl.site1.residues, xl.site2.residues)
            } else {
                xl.residues.clone()
            };
            let blocked = self
                .mods
                .iter()
                .filter(|m| m.kind == ModKind::Fixed && m.target.len() == 1 && reactive.contains(m.target.as_str()))
                .map(|m| m.target.clone())
                .collect::<Vec<_>>();
            if !blocked.is_empty() {
                warnings.push(format!(
                    "Fixed modification on crosslink-reactive residue {}: those residues can't carry a link. \
                     For TMT applied after crosslinking, make TMT on K a Variable modification \
                     (keep it Fixed on the peptide N-terminus, ^).",
                    blocked.join(", ")
                ));
            }
            if xl.two_sites {
                for (name, site) in [("Site 1", &xl.site1), ("Site 2", &xl.site2)] {
                    if !site.residues.chars().all(|c| VALID_AA.contains(c)) {
                        errors.push(format!("Crosslink {name}: residues must be amino acids from {VALID_AA}."));
                    }
                    if site.residues.is_empty() && !site.protein_nterm && !site.protein_cterm {
                        errors.push(format!("Crosslink {name}: choose residues or a protein terminus."));
                    }
                }
            } else if xl.residues.is_empty() || !xl.residues.chars().all(|c| VALID_AA.contains(c)) {
                errors.push(format!("Crosslink residues must be amino acids from {VALID_AA}."));
            }
            if xl.is_cleavable() {
                let mut arms = vec![("", xl.stub_short, xl.stub_long, &xl.stub_alt)];
                if xl.has_second_arm() {
                    arms[0].0 = "arm 1 ";
                    arms.push(("arm 2 ", xl.stub2_short, xl.stub2_long, &xl.stub2_alt));
                }
                for &(arm, short, long, alt) in &arms {
                    if !(short > 0.0 && long > short) {
                        errors.push(format!("Cleavable linker: {arm}stub masses need 0 < short < long."));
                    }
                    match parse_list(alt) {
                        Err(e) => errors.push(format!("Additional {arm}stub masses: {e}")),
                        Ok(alt) if alt.iter().any(|&m| m <= 0.0) => {
                            errors.push(format!("Additional {arm}stub masses must be positive."))
                        }
                        Ok(_) => {}
                    }
                }
                // a cleavage leaves one arm's short stub and the other arm's long stub
                let last = arms.len() - 1;
                for (i, &(arm, short, ..)) in arms.iter().enumerate() {
                    let (other, _, long, _) = arms[last - i];
                    if (short + long - xl.linker_mass()).abs() > 1e-3 {
                        errors.push(format!(
                            "Cleavable linker: {arm}short + {other}long stub ({:.5} Da) must equal the linker mass ({:.5} Da).",
                            short + long,
                            xl.linker_mass()
                        ));
                    }
                }
            }
        }

        (errors, warnings)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn mod_targets() {
        for ok in ["C", "M", "^", "^Q", "$", "$K", "[", "[M", "]"] {
            assert!(check_target(ok).is_ok(), "{ok}");
        }
        for bad in ["", "B", "CC", "^QQ", "x", "^1"] {
            assert!(check_target(bad).is_err(), "{bad}");
        }
    }

    #[test]
    fn json_round_trip_keeps_settings_and_unknown_keys() {
        let mut s = Settings::default();
        s.mods.push((&PRESETS[4]).into());
        s.mods.last_mut().unwrap().max_count = 2;
        s.max_peaks = Some(250);
        s.fragment_tol = Tolerance { unit: TolUnit::Da, lo: -0.5, hi: 0.5 };
        s.tmt = Tmt::Tmt16;
        s.crosslink.enabled = true;
        s.crosslink.reagent = Reagent::Custom;
        s.crosslink.custom_mass = 158.00376;
        s.crosslink.custom_cleavable = true;
        s.crosslink.stub_alt = "85.98264".into();
        s.crosslink.backend = Backend::Reflect;
        s.clip_n_term_met = true;
        let mut v = s.to_json(&["a.mzML".to_string()]);
        v["bruker_config"] = json!({"ms2": {}, "ion_mobility_scale": "linear"});
        v["database"]["enzyme"]["future_key"] = json!(1);

        let back = Settings::from_json(&v).unwrap();
        let v2 = back.to_json(&back.spectra);
        // everything, including keys the GUI does not edit, must survive
        assert_eq!(v, v2);
        assert_eq!(
            v["database"]["variable_mods"]["S"],
            json!([{"mass": 79.966331, "max_count": 2, "name": "UNIMOD:21", "neutral_losses": [97.976896]}])
        );
        assert_eq!(v["database"]["static_mods"]["C"], json!({"mass": 57.021464, "name": "UNIMOD:4"}));
        assert!(v["database"].get("peptide_min_mass").is_none());
        assert_eq!(v["crosslink"]["cleavable"]["alt"], json!([85.98264]));
        assert_eq!(v["database"]["clip_n_term_met"], json!(true));
        assert!(back.linear_mobility, "the linear 1/K0 scale is read");
        let mut linear = Settings::default();
        linear.linear_alignment = true;
        assert_eq!(linear.to_json(&[])["retention_time_alignment"], json!("linear"));
        assert!(Settings::from_json(&linear.to_json(&[])).unwrap().linear_alignment);
        assert!(Settings::default().to_json(&[]).get("retention_time_alignment").is_none());
        // the long-form LFQ table is on unless switched off
        let mut short = Settings::default();
        assert!(short.to_json(&[])["quant"]["lfq_settings"].get("long_form").is_none());
        short.lfq_long_form = false;
        assert_eq!(short.to_json(&[])["quant"]["lfq_settings"]["long_form"], json!(false));
        assert!(!Settings::from_json(&short.to_json(&[])).unwrap().lfq_long_form);
        assert!(Settings::from_json(&Settings::default().to_json(&[])).unwrap().lfq_long_form);

        // localization: off by default, {"decoy_residues": ...} when on
        let mut localized = Settings::default();
        assert!(localized.to_json(&[]).get("localization").is_none());
        localized.localization = true;
        localized.localization_decoys = "al".into();
        assert_eq!(localized.to_json(&[])["localization"], json!({"decoy_residues": "AL"}));
        let reread = Settings::from_json(&localized.to_json(&[])).unwrap();
        assert!(reread.localization && reread.localization_decoys == "AL");
        let core: sage_core::localization::LocalizationSettings =
            serde_json::from_value(localized.to_json(&[])["localization"].clone()).unwrap();
        assert_eq!(core.decoy_residues, "AL");

        let mut entrapment = Settings::default();
        entrapment.entrapment_tag = "ENTRAP_".into();
        assert_eq!(entrapment.to_json(&[])["entrapment"], json!({"tag": "ENTRAP_"}));
        assert_eq!(Settings::from_json(&entrapment.to_json(&[])).unwrap().entrapment_tag, "ENTRAP_");
        assert!(Settings::default().to_json(&[]).get("entrapment").is_none());
        let mut calibrated = back.clone();
        calibrated.linear_mobility = false;
        assert_eq!(calibrated.to_json(&[])["bruker_config"], json!({"ms2": {}}), "the rest of bruker_config stays");
        assert!(Settings::default().to_json(&[]).get("bruker_config").is_none());
        assert!(back.clip_n_term_met);
        // off is Sage's default: not written
        assert!(Settings::default().to_json(&[]).get("database").and_then(|d| d.get("clip_n_term_met")).is_none());
        // diagnostic ions: off by default, true when on, a custom list kept as written
        assert!(Settings::default().to_json(&[]).get("diagnostic_ions").is_none());
        let mut on = Settings::default();
        on.diagnostic_ions = true;
        assert_eq!(on.to_json(&[])["diagnostic_ions"], json!(true));
        let list = json!([{"name": "TMT126", "mz": 126.127726}]);
        let loaded = Settings::from_json(&json!({"diagnostic_ions": list})).unwrap();
        assert!(loaded.diagnostic_ions);
        assert_eq!(loaded.to_json(&[])["diagnostic_ions"], list);
    }

    #[test]
    fn enzymes_and_two_site_chemistry_round_trip() {
        let mut s = Settings::default();
        s.enzymes = vec!["trypsin".into(), "asp-n".into()];
        s.custom_enzymes.push(EnzymeRule { cleave_at: "DE".into(), restrict: "P".into(), c_terminal: true });
        s.crosslink.enabled = true;
        s.crosslink.reagent = Reagent::Edc;
        (s.crosslink.site1, s.crosslink.site2) = Reagent::Edc.sites().unwrap();
        s.crosslink.two_sites = true;
        let v = s.to_json(&[]);
        assert_eq!(
            v["database"]["enzyme"]["enzymes"],
            json!(["trypsin", "asp-n", {"cleave_at": "DE", "restrict": "P", "c_terminal": true}])
        );
        assert!(v["database"]["enzyme"].get("cleave_at").is_none());
        assert_eq!(v["crosslink"]["reagent"], json!("edc"));
        assert_eq!(v["crosslink"]["site1"], json!("K["));
        assert_eq!(v["crosslink"]["site2"], json!("DE]"));
        assert!(v["crosslink"].get("residues").is_none());
        assert_eq!(Settings::from_json(&v).unwrap().to_json(&[]), v);

        // a preset reagent alone implies its sites
        let edc = Settings::from_json(&json!({"crosslink": {"reagent": "dmtmm"}})).unwrap();
        assert!(edc.crosslink.two_sites);
        assert_eq!(edc.crosslink.site2.residues, "DE");
    }

    /// Isotope-labelled proteins are written as `database.labelled_proteins`,
    /// read back, accepted by Sage, and checked against the FASTA
    #[test]
    fn labelled_proteins_round_trip_and_are_checked() {
        use sage_core::labelled::{Isotope, LabelledProteins};
        let dir = std::env::temp_dir().join(format!("sage-gui-labelled-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let fasta = dir.join("two.fasta");
        std::fs::write(&fasta, ">sp|P1|ONE_HUMAN One protein OS=Homo sapiens GN=ONE\nMKR\n>sp|P2|TWO_HUMAN Two\nMKR\n").unwrap();
        let mut s = Settings { fasta: fasta.display().to_string(), ..Default::default() };
        assert!(s.to_json(&[])["database"].get("labelled_proteins").is_none(), "none: not written");
        assert_eq!(
            fasta_proteins(&s.fasta).unwrap(),
            vec![("sp|P1|ONE_HUMAN".to_string(), "One protein".to_string()), ("sp|P2|TWO_HUMAN".to_string(), "Two".to_string())]
        );
        s.labelled_proteins.insert("sp|P1|ONE_HUMAN".into(), Isotope::C13N15);
        s.labelled_proteins.insert("P2".into(), Isotope::N15);
        s.include_unlabelled = false;
        let v = s.to_json(&[]);
        assert_eq!(v["database"]["labelled_proteins"]["proteins"]["sp|P1|ONE_HUMAN"], json!("13C15N"));
        assert_eq!(v["database"]["labelled_proteins"]["include_unlabelled"], json!(false));
        let back = Settings::from_json(&v).unwrap();
        assert_eq!(back.labelled_proteins, s.labelled_proteins);
        assert!(!back.include_unlabelled);
        let core: LabelledProteins = serde_json::from_value(v["database"]["labelled_proteins"].clone()).unwrap();
        assert!(core.validate().is_ok());
        // in the FASTA (by accession or a part of it): no error, a warning about the labels
        let (errors, warnings) = s.check();
        assert!(errors.iter().all(|e| !e.contains("labelled")), "{errors:?}");
        assert!(warnings.iter().any(|w| w.contains("2 protein(s) are searched isotope-labelled only")));
        // a protein not in the FASTA is an error
        s.labelled_proteins.insert("P99".into(), Isotope::C13);
        assert!(s.check().0.iter().any(|e| e.contains("not in the FASTA: P99")));
        std::fs::remove_dir_all(&dir).unwrap();
    }

    /// Database 1 and 2 replace the FASTA file in the JSON, are read back,
    /// accepted by Sage, and checked
    #[test]
    fn crosslink_databases_replace_the_fasta() {
        let dir = std::env::temp_dir().join(format!("sage-gui-two-db-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let (fasta, one, two) = (dir.join("all.fasta"), dir.join("complex1.fasta"), dir.join("complex2.fasta"));
        std::fs::write(&fasta, ">sp|P9|NINE Nine\nMKR\n").unwrap();
        std::fs::write(&one, ">sp|P1|ONE_HUMAN One OS=Homo sapiens\nMKR\n").unwrap();
        std::fs::write(&two, ">sp|P2|TWO_HUMAN Two OS=Homo sapiens\nMKR\n").unwrap();
        let mut s = Settings { fasta: fasta.display().to_string(), ..Default::default() };
        assert!(!s.two_databases());
        assert_eq!(s.fasta_files(), vec![s.fasta.as_str()]);
        assert!(s.to_json(&[])["database"].get("crosslink_databases").is_none());

        // one database set: the FASTA file is out, and the other is missing
        s.database_1 = one.display().to_string();
        assert!(s.two_databases());
        let v = s.to_json(&[]);
        assert!(v["database"].get("fasta").is_none(), "the FASTA file is excluded");
        assert!(s.check().0.iter().any(|e| e.contains("select Database 2")), "{:?}", s.check().0);

        s.database_2 = two.display().to_string();
        s.crosslink.enabled = true;
        let v = s.to_json(&[]);
        assert!(v["database"].get("fasta").is_none());
        assert_eq!(v["database"]["crosslink_databases"]["database_1"], json!(s.database_1));
        assert_eq!(v["database"]["crosslink_databases"]["database_2"], json!(s.database_2));
        let back = Settings::from_json(&v).unwrap();
        assert_eq!((back.database_1.as_str(), back.database_2.as_str()), (s.database_1.as_str(), s.database_2.as_str()));
        // Sage takes it, without a FASTA file
        let builder: sage_core::database::Builder = serde_json::from_value(v["database"].clone()).unwrap();
        let params = builder.make_parameters().unwrap();
        assert_eq!(params.fasta_files(), vec![s.database_1.as_str(), s.database_2.as_str()]);

        let (errors, warnings) = s.check();
        assert!(errors.iter().all(|e| !e.contains("FASTA") && !e.contains("Database")), "{errors:?}");
        assert!(warnings.iter().any(|w| w.contains("only crosslinks between a Database 1 protein")), "{warnings:?}");
        // the proteins of both, e.g. for isotope labels
        s.labelled_proteins.insert("P2".into(), sage_core::labelled::Isotope::N15);
        assert!(s.check().0.iter().all(|e| !e.contains("not in the FASTA")));
        s.labelled_proteins.insert("P9".into(), sage_core::labelled::Isotope::N15);
        assert!(s.check().0.iter().any(|e| e.contains("not in the FASTA: P9")), "the FASTA file is not searched");
        s.labelled_proteins.clear();

        // crosslinking off: one database, with a warning
        s.crosslink.enabled = false;
        assert!(s.check().1.iter().any(|w| w.contains("crosslinking is off")));
        // cleared: back to the FASTA file
        s.database_1.clear();
        s.database_2.clear();
        assert_eq!(s.to_json(&[])["database"]["fasta"], json!(s.fasta));
        std::fs::remove_dir_all(&dir).unwrap();
    }

    /// SILAC is written as `database.silac`, read back, accepted by Sage, and
    /// checked
    #[test]
    fn silac_round_trips_and_is_checked() {
        let mut s = Settings::default();
        assert!(s.to_json(&[])["database"].get("silac").is_none(), "off: not written");
        s.silac = Silac::Triplex;
        s.silac_medium = (4.025107, 6.020129);
        s.silac_heavy = (8.014199, 10.008269);
        s.crosslink.enabled = true;
        s.silac_mixed = true;
        let v = s.to_json(&[]);
        assert_eq!(v["database"]["silac"]["channels"][0]["name"], json!("medium"));
        assert_eq!(v["database"]["silac"]["channels"][1]["labels"]["R"], json!(10.008269));
        let back = Settings::from_json(&v).unwrap();
        assert_eq!(back.silac, Silac::Triplex);
        assert_eq!((back.silac_medium, back.silac_heavy, back.silac_mixed), (s.silac_medium, s.silac_heavy, true));
        // Sage reads it
        let core: sage_core::silac::SilacSettings = serde_json::from_value(v["database"]["silac"].clone()).unwrap();
        assert!(core.validate().is_ok());
        assert_eq!(core.names(), ["light", "medium", "heavy"]);
        let (errors, warnings) = s.check();
        assert!(errors.iter().all(|e| !e.contains("SILAC")), "{errors:?}");
        assert!(warnings.iter().any(|w| w.contains("3 label states") && w.contains("crosslinks.silac.tsv")));
        // a channel without labels is an error; so are two identical channels
        s.silac_heavy = (0.0, 0.0);
        assert!(s.check().0.iter().any(|e| e.contains("heavy channel labels neither")));
        s.silac_heavy = s.silac_medium;
        assert!(s.check().0.iter().any(|e| e.contains("same labels")));
        // duplex with only heavy lysine
        s.silac = Silac::Duplex;
        s.silac_heavy = (8.014199, 0.0);
        let v = s.to_json(&[]);
        assert_eq!(v["database"]["silac"]["channels"].as_array().unwrap().len(), 1);
        assert!(v["database"]["silac"]["channels"][0]["labels"].get("R").is_none());
        assert_eq!(Settings::from_json(&v).unwrap().silac_heavy, (8.014199, 0.0));
    }

    /// Position, motif and neutral losses are written as Sage reads them, and
    /// named / rule keys of a configuration become rows
    #[test]
    fn mass_offsets_and_site_library_round_trip() {
        let mut s = Settings::default();
        s.mods.push(Modification {
            target: "K".into(),
            mass: 176.01433,
            name: "DSSO_OH".into(),
            kind: ModKind::Variable,
            position: SitePosition::NotCTerm,
            losses: "122.00377".into(),
            search: ModSearch::Offset,
            ..Default::default()
        });
        let mut phospho: Modification = PRESETS.iter().find(|p| p.name == "UNIMOD:21" && p.target == "S").unwrap().into();
        phospho.search = ModSearch::Library;
        s.mods.push(phospho);
        let dir = std::env::temp_dir().join(format!("xnova_sites_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let sites = dir.join("sites.tsv");
        std::fs::write(&sites, "protein	position	residue	modification
P1	3	S	UNIMOD:21
").unwrap();
        s.site_library = sites.display().to_string();
        s.site_library_q = 0.05;
        let v = s.to_json(&[]);
        let db = &v["database"];
        assert_eq!(db["variable_mods"]["K"][0]["search_mode"], json!("mass_offset"));
        assert_eq!(db["variable_mods"]["S"][0]["site_mode"], json!("library"));
        assert_eq!(db["site_library"], json!({"path": s.site_library, "max_q": 0.05}));
        // Sage reads it
        let builder: sage_core::database::Builder = serde_json::from_value(json!({
            "fasta": "x",
            "variable_mods": db["variable_mods"].clone(),
            "site_library": db["site_library"].clone(),
        }))
        .unwrap();
        let params = builder.make_parameters().unwrap();
        assert_eq!(params.modifications.offset_mods.len(), 1);
        assert!(params.modifications.has_library_mods());
        let back = Settings::from_json(&v).unwrap();
        assert_eq!(back.mods.iter().find(|m| m.target == "K").unwrap().search, ModSearch::Offset);
        assert_eq!(back.mods.iter().find(|m| m.target == "S").unwrap().search, ModSearch::Library);
        assert_eq!((back.site_library.as_str(), back.site_library_q), (s.site_library.as_str(), 0.05));
        assert!(back.to_json(&[])["database"].get("site_library").is_some());
        let errors = s.check().0;
        assert!(!errors.iter().any(|e| e.contains("offset") || e.contains("library") || e.contains("Site library")), "{errors:?}");

        // checked
        let mut bad = s.clone();
        bad.site_library.clear();
        bad.crosslink.enabled = true;
        bad.mods.push(Modification { target: "K".into(), mass: 176.01433, name: "DSSO_OH".into(), kind: ModKind::Variable, ..Default::default() });
        bad.mods.push(Modification { target: "^".into(), mass: 42.0106, kind: ModKind::Variable, search: ModSearch::Library, ..Default::default() });
        bad.mods.push(Modification { target: "N".into(), mass: 1.0, kind: ModKind::Variable, motif: "N-x-[ST]".into(), search: ModSearch::Offset, ..Default::default() });
        let errors = bad.check().0;
        for needle in ["need a site library", "linear searches", "both indexed and a mass offset", "library sites need a single residue", "cannot have a motif"] {
            assert!(errors.iter().any(|e| e.contains(needle)), "{needle}: {errors:?}");
        }
        // a fixed modification is always indexed
        let mut fixed = Modification { target: "C".into(), mass: 57.021464, kind: ModKind::Fixed, search: ModSearch::Offset, ..Default::default() };
        assert!(mod_entry(&fixed).is_number());
        fixed.kind = ModKind::Variable;
        assert_eq!(mod_entry(&fixed)["search_mode"], json!("mass_offset"));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn modification_rules_round_trip() {
        let mut s = Settings::default();
        s.mods.push(Modification {
            target: "N".into(),
            mass: 0.984016,
            name: "UNIMOD:7".into(),
            kind: ModKind::Variable,
            motif: "N-{P}-[ST]".into(),
            ..Default::default()
        });
        s.mods.push(Modification {
            target: "K".into(),
            mass: 114.042927,
            name: "UNIMOD:121".into(),
            kind: ModKind::Variable,
            position: SitePosition::NotCTerm,
            ..Default::default()
        });
        s.mods.push(PRESETS.iter().find(|p| p.name == "UNIMOD:21" && p.target == "S").unwrap().into());
        let v = s.to_json(&[]);
        let var = &v["database"]["variable_mods"];
        assert_eq!(var["N"][0]["motif"], json!("N-{P}-[ST]"));
        assert_eq!(var["K"][0]["position"], json!("not_c_term"));
        assert_eq!(var["S"][0]["neutral_losses"], json!([97.976896]));
        // Sage reads it
        let builder: sage_core::database::Builder =
            serde_json::from_value(json!({"fasta": "x", "variable_mods": var.clone()})).unwrap();
        assert!(builder.make_parameters().is_ok());
        let back = Settings::from_json(&v).unwrap();
        let n = back.mods.iter().find(|m| m.target == "N").unwrap();
        assert_eq!(n.motif, "N-{P}-[ST]");
        assert_eq!(back.mods.iter().find(|m| m.target == "K").unwrap().position, SitePosition::NotCTerm);
        assert_eq!(parse_losses(&back.mods.iter().find(|m| m.target == "S").unwrap().losses).unwrap(), vec![97.976896]);
        // keys with rules, from a hand-written configuration
        let mut cfg = v.clone();
        cfg["database"]["variable_mods"] = json!({"internal_residue:K": [114.042927], "motif:[ST]-P": [79.966331]});
        let rows = Settings::from_json(&cfg).unwrap().mods;
        assert!(rows.iter().any(|m| m.target == "K" && m.position == SitePosition::Internal));
        assert_eq!(rows.iter().filter(|m| m.motif == "[ST]-P").count(), 2, "S and T");
        // checked
        let mut bad = Settings::default();
        bad.mods.push(Modification { target: "K".into(), mass: 1.0, motif: "N-x-[ST]".into(), ..Default::default() });
        bad.mods.push(Modification { target: "^".into(), mass: 1.0, position: SitePosition::Internal, ..Default::default() });
        bad.mods.push(Modification { target: "S".into(), mass: 1.0, losses: "abc".into(), ..Default::default() });
        let errors = bad.check().0;
        assert_eq!(errors.iter().filter(|e| e.starts_with("Modification `")).count(), 3, "{errors:?}");
    }

    /// Dimethyl presets: Sage's settings, read back, light label not in the table
    #[test]
    fn dimethyl_round_trips_and_is_checked() {
        let mut s = Settings::default();
        s.silac = Silac::DimethylTriplex;
        let v = s.to_json(&[]);
        let silac = &v["database"]["silac"];
        assert_eq!(silac["light"]["^"], json!(28.0313));
        assert_eq!(silac["channels"][1]["labels"]["K"], json!(36.07567));
        let core: sage_core::silac::SilacSettings = serde_json::from_value(silac.clone()).unwrap();
        assert!(core.validate().is_ok());
        assert_eq!(core, sage_core::silac::SilacSettings::dimethyl(true, true));
        assert_eq!(core.names(), ["light", "medium", "heavy"]);
        assert_eq!(Settings::from_json(&v).unwrap().silac, Silac::DimethylTriplex);
        // a search's own configuration lists the light label as fixed: not shown as a modification
        let mut echoed = v.clone();
        echoed["database"]["static_mods"]["K"] = json!({"mass": 28.0313, "name": "UNIMOD:36"});
        echoed["database"]["static_mods"]["^"] = json!({"mass": 28.0313, "name": "UNIMOD:36"});
        let back = Settings::from_json(&echoed).unwrap();
        assert!(back.mods.iter().all(|m| m.target != "K" && m.target != "^"), "{:?}", back.mods);
        assert_eq!(back.mods.len(), s.mods.len());
        let (errors, warnings) = s.check();
        assert!(errors.iter().all(|e| !e.contains("Dimethyl")), "{errors:?}");
        assert!(warnings.iter().any(|w| w.contains("Dimethyl") && w.contains("3 label states")));
        // another fixed modification on K clashes
        s.mods.push(Modification { target: "K".into(), mass: 229.1629, name: String::new(), kind: ModKind::Fixed, max_count: 0, ..Default::default() });
        assert!(s.check().0.iter().any(|e| e.contains("Dimethyl") && e.contains('K')));
        // duplex: light and medium
        s.silac = Silac::DimethylDuplex;
        let v = s.to_json(&[]);
        assert_eq!(v["database"]["silac"]["channels"].as_array().unwrap().len(), 1);
        assert_eq!(v["database"]["silac"]["channels"][0]["name"], json!("medium"));
        assert_eq!(Settings::from_json(&v).unwrap().silac, Silac::DimethylDuplex);
    }

    #[test]
    fn multivalent_option_round_trips_and_warns() {
        let mut s = Settings::default();
        s.crosslink.enabled = true;
        // off: not written, Sage's default
        assert!(s.to_json(&[])["crosslink"].get("multivalent").is_none());
        s.crosslink.multivalent = true;
        s.crosslink.valency = 3;
        let v = s.to_json(&[]);
        assert_eq!(v["crosslink"]["multivalent"], json!(3));
        let back = Settings::from_json(&v).unwrap();
        assert!(back.crosslink.multivalent);
        assert_eq!(back.crosslink.valency, 3);
        // Sage reads the same setting
        let core: sage_core::crosslink::CrosslinkSettings = serde_json::from_value(v["crosslink"].clone()).unwrap();
        assert_eq!(core.multivalent, Some(3));
        assert!(core.validate().is_ok());

        // the search-space warning, for non-cleavable and MS-cleavable linkers alike
        for reagent in [Reagent::Dss, Reagent::Dsso] {
            s.crosslink.reagent = reagent;
            s.crosslink.reset_stubs();
            let (errors, warnings) = s.check();
            assert!(errors.iter().all(|e| !e.contains("ultivalent")), "{reagent:?}: {errors:?}");
            assert!(warnings.iter().any(|w| w.contains("search space is very large") && w.contains("small protein subset")));
            let core: sage_core::crosslink::CrosslinkSettings =
                serde_json::from_value(s.to_json(&[])["crosslink"].clone()).unwrap();
            assert!(core.validate().is_ok(), "{reagent:?}");
        }
    }

    #[test]
    fn classic_enzyme_blocks_map_to_presets() {
        let enzymes = |e: Value| {
            let s = Settings::from_json(&json!({"database": {"enzyme": e}})).unwrap();
            (s.enzymes, s.custom_enzymes, s.nonspecific)
        };
        assert_eq!(enzymes(json!({"cleave_at": "KR", "restrict": "P"})).0, vec!["trypsin"]);
        assert_eq!(enzymes(json!({"cleave_at": "D", "c_terminal": false})).0, vec!["asp-n"]);
        // Sage's defaults inside an enzyme block: KR without restriction
        assert_eq!(enzymes(json!({"missed_cleavages": 2})).0, vec!["trypsin/p"]);
        assert!(enzymes(json!({"cleave_at": ""})).2);
        let (names, custom, _) = enzymes(json!({"cleave_at": "$"}));
        assert!(names.is_empty() && custom[0].cleave_at == "$");
        assert_eq!(enzymes(json!({"enzymes": ["Lys-C", "GluC"]})).0, vec!["lys-c", "glu-c"]);

        // a single custom rule is written back in the classic form (so `$` works)
        let s = Settings::from_json(&json!({"database": {"enzyme": {"cleave_at": "$"}}})).unwrap();
        assert_eq!(s.to_json(&[])["database"]["enzyme"]["cleave_at"], json!("$"));
        assert!(s.check().0.iter().all(|e| !e.contains("Custom enzyme")));
    }

    #[test]
    fn preset_stub_masses_are_editable() {
        let mut s = Settings::default();
        s.crosslink.enabled = true;
        s.crosslink.reagent = Reagent::Dsso;
        s.crosslink.reset_stubs();
        // unedited preset stubs are implied by the reagent
        let v = s.to_json(&[]);
        assert!(v["crosslink"].get("cleavable").is_none());
        let back = Settings::from_json(&v).unwrap();
        assert_eq!((back.crosslink.stub_short, back.crosslink.stub_long), (54.01056, 103.9932));
        assert_eq!(back.crosslink.stub_alt, "85.98264");
        assert!(!back.crosslink.stubs_edited());

        // dropping the sulfenic acid stub writes an override next to the reagent
        s.crosslink.stub_alt.clear();
        let v = s.to_json(&[]);
        assert_eq!(v["crosslink"]["reagent"], json!("dsso"));
        assert_eq!(v["crosslink"]["cleavable"], json!({"short": 54.01056, "long": 103.9932, "alt": []}));
        let back = Settings::from_json(&v).unwrap();
        assert_eq!(back.crosslink.reagent, Reagent::Dsso);
        assert!(back.crosslink.stub_alt.is_empty());
        assert_eq!(back.to_json(&[]), v);

        // stubs must still add up to the linker mass
        s.crosslink.stub_long = 100.0;
        assert!(s.check().0.iter().any(|e| e.contains("must equal the linker mass")));
    }

    #[test]
    fn asymmetric_linkers_have_two_arms() {
        let mut s = Settings::default();
        s.crosslink.enabled = true;
        s.crosslink.reagent = Reagent::Clixlink;
        s.crosslink.reset_stubs();
        assert!(s.crosslink.has_second_arm());
        assert!(s.check().0.iter().all(|e| !e.contains("Cleavable")), "{:?}", s.check().0);
        // the preset is implied by the reagent
        let v = s.to_json(&[]);
        assert_eq!(v["crosslink"]["reagent"], json!("clixlink"));
        assert!(v["crosslink"].get("cleavable").is_none());

        // an edited second arm writes both arms
        s.crosslink.stub2_alt.clear();
        let v = s.to_json(&[]);
        assert_eq!(v["crosslink"]["cleavable"]["other"], json!({"short": 163.06333, "long": 213.04596, "alt": []}));
        let back = Settings::from_json(&v).unwrap();
        assert_eq!(back.crosslink.reagent, Reagent::Clixlink);
        assert_eq!(back.crosslink.stub_alt, "85.98264");
        assert!(back.crosslink.stub2_alt.is_empty());
        assert_eq!(back.to_json(&[]), v);

        // each arm's short stub + the other arm's long stub must make the linker
        s.crosslink.stub2_long = 200.0;
        assert!(s.check().0.iter().any(|e| e.contains("arm 1 short + arm 2 long")));

        // custom asymmetric linkers
        let mut s = Settings::default();
        s.crosslink.enabled = true;
        s.crosslink.reagent = Reagent::Custom;
        s.crosslink.custom_mass = 267.05653;
        s.crosslink.custom_cleavable = true;
        s.crosslink.asymmetric = true;
        s.crosslink.stub_short = 54.01056;
        s.crosslink.stub_long = 103.9932;
        let v = s.to_json(&[]);
        assert!(v["crosslink"]["cleavable"].get("other").is_some());
        let back = Settings::from_json(&v).unwrap();
        assert!(back.crosslink.custom_cleavable && back.crosslink.asymmetric);
        assert_eq!(back.to_json(&[]), v);
    }

    #[test]
    fn missing_keys_take_sage_defaults() {
        let s = Settings::from_json(&json!({"database": {"fasta": "x.fasta"}})).unwrap();
        assert_eq!(s.missed_cleavages, 0);
        assert!(s.mods.is_empty());
        assert!(!s.crosslink.enabled);
        assert_eq!(s.precursor_charge, (2, 4));
    }

    /// Run the release `sage.exe` on configurations written by the GUI, with
    /// native (e.g. `C:\...`) paths. Skipped when sage has not been built.
    #[test]
    fn generated_configs_run_in_sage() {
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..").join("..").canonicalize().unwrap();
        let sage = root.join("target").join("release").join(if cfg!(windows) { "sage.exe" } else { "sage" });
        if !sage.is_file() {
            eprintln!("skipped: {} not built", sage.display());
            return;
        }
        let data = root.join("tests");
        let mut cases = vec![("linear", Settings::default())];
        let mut dsso = Settings::default();
        dsso.crosslink.enabled = true;
        dsso.crosslink.reagent = Reagent::Dsso;
        cases.push(("dsso", dsso.clone()));
        // edited preset stubs, wider charge range
        dsso.crosslink.stub_alt.clear();
        dsso.precursor_charge = (3, 6);
        cases.push(("dsso_edited_stubs", dsso));
        let mut reflect = Settings::default();
        reflect.crosslink.enabled = true;
        reflect.crosslink.backend = Backend::Reflect;
        cases.push(("reflect", reflect));
        let mut edc = Settings::default();
        edc.enzymes = vec!["trypsin".into(), "asp-n".into(), "glu-c".into()];
        edc.crosslink.enabled = true;
        edc.crosslink.reagent = Reagent::Edc;
        (edc.crosslink.site1, edc.crosslink.site2) = Reagent::Edc.sites().unwrap();
        edc.crosslink.two_sites = true;
        cases.push(("edc_codigest", edc));
        let mut multi = Settings::default();
        multi.enzymes = vec!["trypsin".into(), "lys-c".into(), "chymotrypsin".into(), "proalanase".into()];
        cases.push(("multi_enzyme", multi));
        let mut quant = Settings::default();
        quant.crosslink.enabled = true;
        quant.lfq = true;
        quant.tmt = Tmt::Tmt6;
        quant.tmt_level = 2;
        cases.push(("xl_quant", quant));
        let mut clixlink = Settings::default();
        clixlink.crosslink.enabled = true;
        clixlink.crosslink.reagent = Reagent::Clixlink;
        clixlink.crosslink.reset_stubs();
        clixlink.crosslink.stub2_alt.clear();
        cases.push(("clixlink_edited_stubs", clixlink));
        let mut with_linear = Settings::default();
        with_linear.crosslink.enabled = true;
        with_linear.crosslink.linear_search = true;
        cases.push(("xl_with_linear", with_linear));

        for (name, mut s) in cases {
            let out = std::env::temp_dir().join(format!("sage-gui-test-{name}"));
            let _ = std::fs::remove_dir_all(&out);
            s.fasta = data.join("Q99536.fasta").display().to_string();
            s.spectra = vec![data.join("LQSRPAAPPAPGPGQLTLR.mzML").display().to_string()];
            s.output_directory = out.display().to_string();
            let (errors, _) = s.check();
            assert!(errors.is_empty(), "{name}: {errors:?}");

            std::fs::create_dir_all(&out).unwrap();
            let config = out.join("sage_search.json");
            std::fs::write(&config, serde_json::to_string_pretty(&s.to_json(&s.spectra)).unwrap()).unwrap();
            let run = std::process::Command::new(&sage)
                .arg(&config)
                .args(["--threads", "2"])
                .arg("--disable-telemetry-i-dont-want-to-improve-sage")
                .output()
                .unwrap();
            assert!(run.status.success(), "{name}: {}", String::from_utf8_lossy(&run.stderr));
            let results = if s.crosslink.enabled { "crosslinks.sage.tsv" } else { "results.sage.tsv" };
            let tsv = std::fs::read_to_string(out.join(results)).unwrap();
            if name == "linear" || name == "xl_with_linear" {
                let linear = std::fs::read_to_string(out.join("results.sage.tsv")).unwrap();
                assert!(linear.contains("LQSRPAAPPAPGPGQLTLR"), "{name}: linear search should identify the peptide");
            }
            if name == "xl_quant" {
                // reporter columns on every crosslink match, and an LFQ table
                let header = tsv.lines().next().unwrap();
                assert!(header.ends_with("ion_injection_time\ttmt_1\ttmt_2\ttmt_3\ttmt_4\ttmt_5\ttmt_6"), "{header}");
                let lfq = std::fs::read_to_string(out.join("crosslinks.lfq.tsv")).unwrap();
                assert!(lfq.starts_with("peptide_a\tpeptide_b"));
            }
            let _ = std::fs::remove_dir_all(&out);
        }
    }
}
