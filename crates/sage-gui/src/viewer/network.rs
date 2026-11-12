//! Protein interaction networks of a search's crosslinks: each protein a
//! node, each pair of proteins crosslinked to each other an edge, wider the
//! more crosslinks join them. With what UniProt and STRING know about the
//! proteins ([`super::protein_info`]), the nodes are grouped and coloured by
//! complex (ComplexPortal), pathway (Reactome), GO term, subcellular location
//! or keyword, or by the network's own modules (Louvain) or the user's
//! groups; the edges by whether STRING knows the interaction, or by the
//! complex, pathway or group the two proteins share.
//!
//! Drawn by XNova itself (SVG, PDF, PNG, TIFF), and written for Cytoscape
//! (`.cyjs`, with the layout), as GraphML and as node and edge tables, every
//! annotation looked up with them. Everything about the figure is in
//! [`NetworkStyle`], saved and loaded as JSON.

use super::circos::{self, Background, Palette, TypeStyle};
use super::figure;
use super::mobiligram::Class;
use super::model::{Kind, Results};
use super::protein_info::{Annotations, ProteinInfo, Term};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt::Write;
use std::path::{Path, PathBuf};

/// Folder of a results folder the networks of a search go to
pub const FOLDER: &str = "network";

// ---------------------------------------------------------------- the style

/// What groups (and colours) the proteins
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum GroupBy {
    None,
    Complex,
    Pathway,
    Process,
    Function,
    Component,
    Location,
    Keyword,
    /// The network's own modules (Louvain)
    Module,
    /// A file of the user's: protein, group, colour
    Own,
}

impl GroupBy {
    pub const ALL: [GroupBy; 10] = [
        GroupBy::Complex,
        GroupBy::Pathway,
        GroupBy::Process,
        GroupBy::Function,
        GroupBy::Component,
        GroupBy::Location,
        GroupBy::Keyword,
        GroupBy::Module,
        GroupBy::Own,
        GroupBy::None,
    ];

    pub fn label(self) -> &'static str {
        match self {
            GroupBy::None => "nothing (one colour)",
            GroupBy::Complex => "complex (ComplexPortal)",
            GroupBy::Pathway => "pathway (Reactome)",
            GroupBy::Process => "GO biological process",
            GroupBy::Function => "GO molecular function",
            GroupBy::Component => "GO cellular component",
            GroupBy::Location => "subcellular location (UniProt)",
            GroupBy::Keyword => "UniProt keyword",
            GroupBy::Module => "network module (Louvain)",
            GroupBy::Own => "your own groups (file)",
        }
    }

    /// Its heading in the legend
    pub fn heading(self) -> &'static str {
        match self {
            GroupBy::None => "Proteins",
            GroupBy::Complex => "Complex",
            GroupBy::Pathway => "Pathway",
            GroupBy::Process => "Biological process",
            GroupBy::Function => "Molecular function",
            GroupBy::Component => "Cellular component",
            GroupBy::Location => "Location",
            GroupBy::Keyword => "Keyword",
            GroupBy::Module => "Network module",
            GroupBy::Own => "Group",
        }
    }

    /// Needs the proteins looked up in UniProt
    pub fn looked_up(self) -> bool {
        !matches!(self, GroupBy::None | GroupBy::Module | GroupBy::Own)
    }

    /// A protein's terms of this kind
    pub fn terms(self, p: &ProteinInfo) -> Vec<Term> {
        match self {
            GroupBy::Complex => p.complexes.clone(),
            GroupBy::Pathway => p.pathways.clone(),
            GroupBy::Process => p.go_process.clone(),
            GroupBy::Function => p.go_function.clone(),
            GroupBy::Component => p.go_component.clone(),
            GroupBy::Location => p.locations.iter().map(|l| Term { id: String::new(), name: l.clone() }).collect(),
            GroupBy::Keyword => p.keywords.clone(),
            GroupBy::None | GroupBy::Module | GroupBy::Own => Vec::new(),
        }
    }
}

/// What gives an edge its colour
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EdgeColor {
    /// Known to STRING, or new
    Known,
    /// A complex both proteins are in
    Complex,
    /// A pathway both proteins are in
    Pathway,
    /// Both proteins in one group
    Group,
    Single,
}

impl EdgeColor {
    pub const ALL: [EdgeColor; 5] = [EdgeColor::Known, EdgeColor::Complex, EdgeColor::Pathway, EdgeColor::Group, EdgeColor::Single];

    pub fn label(self) -> &'static str {
        match self {
            EdgeColor::Known => "known (STRING) or new",
            EdgeColor::Complex => "shared complex",
            EdgeColor::Pathway => "shared pathway",
            EdgeColor::Group => "within a group",
            EdgeColor::Single => "one colour",
        }
    }
}

/// Which of STRING's evidence makes an interaction known
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Evidence {
    /// The combined score (all channels, text mining too)
    Combined,
    Experiments,
    /// Curated databases (pathways, complexes)
    Databases,
    ExperimentsOrDatabases,
}

impl Evidence {
    pub const ALL: [Evidence; 4] = [Evidence::Combined, Evidence::Experiments, Evidence::Databases, Evidence::ExperimentsOrDatabases];

    pub fn label(self) -> &'static str {
        match self {
            Evidence::Combined => "any (combined score)",
            Evidence::Experiments => "experiments",
            Evidence::Databases => "curated databases",
            Evidence::ExperimentsOrDatabases => "experiments or databases",
        }
    }

    pub fn score(self, i: &super::protein_info::Interaction) -> f32 {
        match self {
            Evidence::Combined => i.score,
            Evidence::Experiments => i.experiments,
            Evidence::Databases => i.databases,
            Evidence::ExperimentsOrDatabases => i.experiments.max(i.databases),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EdgeWidth {
    Crosslinks,
    Csms,
    Fixed,
}

impl EdgeWidth {
    pub const ALL: [EdgeWidth; 3] = [EdgeWidth::Crosslinks, EdgeWidth::Csms, EdgeWidth::Fixed];

    pub fn label(self) -> &'static str {
        match self {
            EdgeWidth::Crosslinks => "crosslinks (residue pairs)",
            EdgeWidth::Csms => "CSMs",
            EdgeWidth::Fixed => "all the same",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum NodeSize {
    Crosslinks,
    Csms,
    /// Proteins it is crosslinked to
    Partners,
    Fixed,
}

impl NodeSize {
    pub const ALL: [NodeSize; 4] = [NodeSize::Crosslinks, NodeSize::Csms, NodeSize::Partners, NodeSize::Fixed];

    pub fn label(self) -> &'static str {
        match self {
            NodeSize::Crosslinks => "crosslinks",
            NodeSize::Csms => "CSMs",
            NodeSize::Partners => "partners",
            NodeSize::Fixed => "all the same",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum NodeLabel {
    /// Gene name (UniProt's, else the FASTA's), else entry name
    Gene,
    /// `RPB1_YEAST`
    Entry,
    Accession,
    None,
}

impl NodeLabel {
    pub const ALL: [NodeLabel; 4] = [NodeLabel::Gene, NodeLabel::Entry, NodeLabel::Accession, NodeLabel::None];

    pub fn label(self) -> &'static str {
        match self {
            NodeLabel::Gene => "gene name",
            NodeLabel::Entry => "entry name",
            NodeLabel::Accession => "accession",
            NodeLabel::None => "none",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum LabelPlace {
    Center,
    Beside,
    Below,
}

impl LabelPlace {
    pub const ALL: [LabelPlace; 3] = [LabelPlace::Center, LabelPlace::Beside, LabelPlace::Below];

    pub fn label(self) -> &'static str {
        match self {
            LabelPlace::Center => "on the node",
            LabelPlace::Beside => "beside it",
            LabelPlace::Below => "below it",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Layout {
    /// Force-directed (Fruchterman-Reingold)
    Force,
    /// Force-directed, each group a cluster of its own
    Groups,
    /// A circle, group by group
    Circle,
}

impl Layout {
    pub const ALL: [Layout; 3] = [Layout::Force, Layout::Groups, Layout::Circle];

    pub fn label(self) -> &'static str {
        match self {
            Layout::Force => "force-directed",
            Layout::Groups => "clusters by group",
            Layout::Circle => "circle",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum LegendSide {
    Right,
    Left,
    Below,
}

impl LegendSide {
    pub const ALL: [LegendSide; 3] = [LegendSide::Right, LegendSide::Left, LegendSide::Below];

    pub fn label(self) -> &'static str {
        match self {
            LegendSide::Right => "right",
            LegendSide::Left => "left",
            LegendSide::Below => "below",
        }
    }
}

/// Every choice about a network. Missing keys of a JSON style take their default.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct NetworkStyle {
    // which crosslinks (as the crosslink plot)
    /// Highest CSM and unique-crosslink q-value
    pub max_q: f32,
    /// Fewest CSMs of a crosslink (residue pair)
    pub min_csms: u32,
    /// Crosslinks within a protein (homomeric or not) and loop-links, as a
    /// loop on its node
    pub intra: bool,
    pub loops: bool,
    /// Each link of a multivalent species, as an edge
    pub multivalent: bool,
    /// A peptide found in several proteins linked in each of them
    pub all_places: bool,

    // the proteins
    /// Accessions left out (with their edges)
    pub hidden: Vec<String>,
    /// Proteins crosslinked only within themselves, as nodes without edges
    pub lone: bool,
    pub label_by: NodeLabel,
    /// Accession → label of one's own
    pub labels: BTreeMap<String, String>,
    pub label_place: LabelPlace,
    pub size_by: NodeSize,
    /// Smallest and largest node, diameter in points
    pub node_size: f64,
    pub node_size_max: f64,
    pub node_outline: String,

    // the groups
    pub group_by: GroupBy,
    /// The file of [`GroupBy::Own`]: protein, group and an optional colour per line
    pub own_groups: String,
    /// Fewest proteins of a group
    pub min_group: usize,
    /// Most groups coloured (the largest)
    pub max_groups: usize,
    /// Terms not to group by (their keys: `GO:0005634`, `Nucleus`), e.g. too general
    pub excluded: Vec<String>,
    /// Terms to group by first, in this order (their keys), e.g. to name a
    /// group by the pathway one prefers of several shared by the same proteins
    pub preferred: Vec<String>,
    /// Term key → label of one's own
    pub group_labels: BTreeMap<String, String>,
    /// Term key → colour
    pub group_colors: BTreeMap<String, String>,
    pub palette: Palette,
    /// Proteins in no group
    pub ungrouped_color: String,
    /// The groups' names next to them in the figure
    pub group_names: bool,

    // the edges
    pub color_by: EdgeColor,
    pub known_color: String,
    pub new_color: String,
    pub edge_color: String,
    /// Edges sharing no complex, pathway or group
    pub other_color: String,
    pub evidence: Evidence,
    /// Lowest STRING score (0 to 1) of a known interaction
    pub min_score: f32,
    pub width_by: EdgeWidth,
    /// Thinnest and widest edge, points
    pub edge_width: f64,
    pub edge_width_max: f64,
    pub edge_opacity: f64,
    /// The number of crosslinks on each edge
    pub numbers: bool,
    pub self_color: String,

    // the layout
    pub layout: Layout,
    /// Distance between nodes, times the default
    pub spacing: f64,
    /// How strongly the proteins of a group gather (0: not at all)
    pub group_pull: f64,
    /// Another seed, another layout
    pub seed: u32,
    /// Proteins moved by hand: accession → place (0 to 1 across and down the
    /// drawing area)
    pub positions: BTreeMap<String, [f64; 2]>,

    // the canvas
    /// Inches
    pub width: f64,
    pub height: f64,
    pub background: Background,
    pub title: String,
    /// Font families in order of preference, as in CSS
    pub font: String,
    pub title_text: TypeStyle,
    pub label_text: TypeStyle,
    pub group_text: TypeStyle,
    pub legend_text: TypeStyle,
    pub number_text: TypeStyle,
    pub legend: bool,
    pub legend_side: LegendSide,
    /// The numbers under the legend: proteins, pairs, known, crosslinks, CSMs
    pub summary: bool,
}

impl Default for NetworkStyle {
    fn default() -> Self {
        NetworkStyle {
            max_q: 0.01,
            min_csms: 1,
            intra: true,
            loops: true,
            multivalent: true,
            all_places: false,
            hidden: Vec::new(),
            lone: true,
            label_by: NodeLabel::Gene,
            labels: BTreeMap::new(),
            label_place: LabelPlace::Beside,
            size_by: NodeSize::Crosslinks,
            node_size: 9.0,
            node_size_max: 24.0,
            node_outline: "#4d4d4d".into(),
            group_by: GroupBy::Complex,
            own_groups: String::new(),
            min_group: 2,
            max_groups: 12,
            excluded: Vec::new(),
            preferred: Vec::new(),
            group_labels: BTreeMap::new(),
            group_colors: BTreeMap::new(),
            palette: Palette::Tableau20,
            ungrouped_color: "#ffffff".into(),
            group_names: true,
            color_by: EdgeColor::Known,
            known_color: "#e15759".into(),
            new_color: "#8c8c8c".into(),
            edge_color: "#7f7f7f".into(),
            other_color: "#bdbdbd".into(),
            evidence: Evidence::ExperimentsOrDatabases,
            min_score: 0.4,
            width_by: EdgeWidth::Crosslinks,
            edge_width: 0.8,
            edge_width_max: 4.0,
            edge_opacity: 0.85,
            numbers: false,
            self_color: "#7f7f7f".into(),
            layout: Layout::Force,
            spacing: 1.0,
            group_pull: 0.3,
            seed: 1,
            positions: BTreeMap::new(),
            width: 8.0,
            height: 6.0,
            background: Background::White,
            title: String::new(),
            font: "Arial, Helvetica, sans-serif".into(),
            title_text: TypeStyle::new(15.0, true),
            label_text: TypeStyle::new(8.0, false),
            group_text: TypeStyle::new(10.0, true),
            legend_text: TypeStyle::new(8.0, false),
            number_text: TypeStyle::new(6.0, false),
            legend: true,
            legend_side: LegendSide::Right,
            summary: true,
        }
    }
}

impl NetworkStyle {
    /// What the network's data depends on (the rest only changes the drawing)
    pub fn selection(&self) -> (f32, u32, bool, bool, bool, bool) {
        (self.max_q, self.min_csms, self.intra, self.loops, self.multivalent, self.all_places)
    }

    /// What the groups and the layout depend on, beside the data and the lookups
    pub fn plan_key(&self) -> String {
        format!(
            "{:?}",
            (
                (&self.hidden, self.lone),
                self.group_by,
                &self.own_groups,
                self.min_group,
                self.max_groups,
                (&self.excluded, &self.preferred),
                self.layout,
                self.spacing,
                self.group_pull,
                self.seed,
                (self.width / self.height.max(0.1) * 100.0).round(),
                (self.legend, self.legend_side, self.label_place, self.title.is_empty()),
            )
        )
    }

    pub fn to_json(&self) -> String {
        serde_json::to_string_pretty(self).unwrap_or_default()
    }

    pub fn from_json(text: &str) -> Result<NetworkStyle, String> {
        serde_json::from_str(text).map_err(|e| format!("not a network style: {e}"))
    }

    /// Every text size times `factor`
    pub fn scale_text(&mut self, factor: f64) {
        for t in [&mut self.title_text, &mut self.label_text, &mut self.group_text, &mut self.legend_text, &mut self.number_text] {
            t.size = ((t.size * factor) * 10.0).round() / 10.0;
            t.size = t.size.clamp(2.0, 96.0);
        }
    }
}

// ---------------------------------------------------------------- the data

/// A protein of the network
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Node {
    pub accession: String,
    /// The FASTA header's gene (`GN=`)
    pub gene: String,
    /// The FASTA header's protein name
    pub name: String,
    /// `RPB1_YEAST` (UniProt headers), else the accession
    pub entry: String,
    /// Crosslinks (residue pairs) to other proteins, and their CSMs
    pub crosslinks: usize,
    pub csms: u32,
    /// Proteins it is crosslinked to
    pub partners: usize,
    /// Crosslinks within it (and loop-links), and their CSMs
    pub self_links: usize,
    pub self_csms: u32,
}

/// Two proteins crosslinked to each other (`a` < `b`)
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Edge {
    pub a: String,
    pub b: String,
    /// Residue pairs
    pub crosslinks: usize,
    pub csms: u32,
    /// Lowest crosslink q-value, highest score
    pub q: f32,
    pub score: f64,
    /// Residue pairs of multivalent species
    pub multivalent: usize,
    /// The residue pairs: residue in `a`, residue in `b` (1-based)
    pub sites: Vec<(usize, usize)>,
}

/// What a network shows, built from a results folder
#[derive(Clone, Debug, Default, PartialEq)]
pub struct NetworkData {
    /// Most crosslinked first
    pub nodes: Vec<Node>,
    pub edges: Vec<Edge>,
}

impl NetworkData {
    #[cfg(test)]
    pub fn node(&self, accession: &str) -> Option<&Node> {
        self.nodes.iter().find(|n| n.accession == accession)
    }

    /// Every accession, for the lookups
    pub fn accessions(&self) -> Vec<String> {
        self.nodes.iter().map(|n| n.accession.clone()).collect()
    }
}

/// The entry name of a UniProt header's accession (`sp|P04050|RPB1_YEAST`)
fn entry_of(accession: &str) -> String {
    let parts: Vec<&str> = accession.split('|').collect();
    match parts.as_slice() {
        ["sp" | "tr", _, entry, ..] if !entry.is_empty() => entry.to_string(),
        _ => accession.to_string(),
    }
}

/// The network of the crosslinks of a results folder
pub fn build(results: &Results, style: &NetworkStyle) -> NetworkData {
    let selection = circos::CircosStyle {
        max_q: style.max_q,
        min_csms: style.min_csms,
        inter: true,
        intra: style.intra,
        loops: style.loops,
        multivalent: style.multivalent,
        all_places: style.all_places,
        coverage: false,
        ..circos::CircosStyle::default()
    };
    let plot = circos::build(results, &selection);
    from_links(&plot.links, |accession| match results.protein(accession) {
        Some(p) => (p.gene.clone(), p.name.clone()),
        None => (String::new(), String::new()),
    })
}

/// The network of crosslinks between residues; `about` gives a protein's
/// gene and name
pub fn from_links(links: &[circos::Link], about: impl Fn(&str) -> (String, String)) -> NetworkData {
    let mut nodes: BTreeMap<String, Node> = BTreeMap::new();
    let mut edges: BTreeMap<(String, String), Edge> = BTreeMap::new();
    let node = |accession: &str| -> Node {
        let (gene, name) = about(accession);
        Node { accession: accession.to_string(), gene, name, entry: entry_of(accession), ..Default::default() }
    };
    for l in links {
        let (a, b) = (&l.a.0, &l.b.0);
        if a == b {
            let n = nodes.entry(a.clone()).or_insert_with(|| node(a));
            n.self_links += 1;
            n.self_csms += l.csms;
            continue;
        }
        let ((a, ra), (b, rb)) = if a <= b { ((a, l.a.1), (b, l.b.1)) } else { ((b, l.b.1), (a, l.a.1)) };
        for p in [a, b] {
            let n = nodes.entry(p.clone()).or_insert_with(|| node(p));
            n.crosslinks += 1;
            n.csms += l.csms;
        }
        let e = edges.entry((a.clone(), b.clone())).or_insert_with(|| Edge {
            a: a.clone(),
            b: b.clone(),
            q: f32::MAX,
            score: f64::MIN,
            ..Default::default()
        });
        e.crosslinks += 1;
        e.csms += l.csms;
        e.q = e.q.min(l.q);
        e.score = e.score.max(l.score);
        if l.class == Class::Multivalent {
            e.multivalent += 1;
        }
        e.sites.push((ra, rb));
    }
    for e in edges.values() {
        for p in [&e.a, &e.b] {
            if let Some(n) = nodes.get_mut(p) {
                n.partners += 1;
            }
        }
    }
    let mut nodes: Vec<Node> = nodes.into_values().collect();
    nodes.sort_by(|x, y| (y.crosslinks + y.self_links).cmp(&(x.crosslinks + x.self_links)).then(x.accession.cmp(&y.accession)));
    let mut edges: Vec<Edge> = edges.into_values().collect();
    for e in &mut edges {
        e.sites.sort_unstable();
    }
    NetworkData { nodes, edges }
}

// ---------------------------------------------------------------- the user's groups

/// The groups of a file of the user's: each protein's group, and colours
#[derive(Clone, Debug, Default, PartialEq)]
pub struct OwnGroups {
    /// Protein (accession, UniProt accession, gene, entry name or label,
    /// upper case) → group
    pub of: HashMap<String, String>,
    pub colors: BTreeMap<String, String>,
    /// Problems with the file, for the user
    pub notes: Vec<String>,
}

impl OwnGroups {
    /// A tab- or comma-separated file: protein, group and an optional
    /// `#rrggbb` colour per line; `#` comments and a header line
    /// (`protein`) are skipped
    pub fn read(path: &str) -> OwnGroups {
        let mut own = OwnGroups::default();
        if path.trim().is_empty() {
            return own;
        }
        let text = match std::fs::read_to_string(path.trim()) {
            Ok(t) => t,
            Err(e) => {
                own.notes.push(format!("{path}: {e}"));
                return own;
            }
        };
        for (i, line) in text.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = if line.contains('\t') { line.split('\t') } else { line.split(',') }.map(str::trim).collect();
            if fields.len() < 2 || fields[0].is_empty() || fields[1].is_empty() {
                own.notes.push(format!("line {}: needs a protein and a group", i + 1));
                continue;
            }
            if i == 0 && fields[0].eq_ignore_ascii_case("protein") {
                continue;
            }
            own.of.insert(fields[0].to_uppercase(), fields[1].to_string());
            if let Some(color) = fields.get(2).filter(|c| c.starts_with('#') && c.len() == 7) {
                own.colors.insert(fields[1].to_string(), color.to_string());
            }
        }
        own
    }

    /// A protein's group, by any of its names
    fn group(&self, node: &Node, info: Option<&ProteinInfo>, label: &str) -> Option<&String> {
        let mut names = vec![node.accession.clone(), node.entry.clone(), node.gene.clone(), label.to_string()];
        if let Some(p) = info {
            names.extend([p.accession.clone(), p.gene.clone(), p.entry.clone()]);
        }
        names.iter().filter(|n| !n.is_empty()).find_map(|n| self.of.get(&n.to_uppercase()))
    }
}

// ---------------------------------------------------------------- groups and modules

/// A group of proteins, coloured
#[derive(Clone, Debug, PartialEq)]
pub struct Group {
    /// The term's key (`CPX-2662`, `GO:0005634`, `Nucleus`, `module:…`)
    pub key: String,
    /// The term's name (`DNA-directed RNA polymerase II complex`)
    pub name: String,
    /// What the figure calls it (the user's label, else the name)
    pub label: String,
    pub color: String,
    pub members: usize,
}

/// The greedy cover of items by terms: the `preferred` terms first (in
/// their order), then the term most items have takes them, then the one
/// most of the rest have, while it takes at least `min` and at most `max`
/// terms are taken. Returns the terms taken (key, name, items taken), each
/// item's term, and the terms not taken (key, name, items having it)
fn cover(items: &[Vec<Term>], min: usize, max: usize, preferred: &[String]) -> (Vec<(String, String, usize)>, Vec<Option<usize>>, Vec<(String, String, usize)>) {
    let mut names: BTreeMap<&str, &str> = BTreeMap::new();
    let mut having: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
    for (i, terms) in items.iter().enumerate() {
        let mut seen = HashSet::new();
        for t in terms {
            if seen.insert(t.key()) {
                names.entry(t.key()).or_insert(&t.name);
                having.entry(t.key()).or_default().push(i);
            }
        }
    }
    let mut of: Vec<Option<usize>> = vec![None; items.len()];
    let mut taken: Vec<(String, String, usize)> = Vec::new();
    while taken.len() < max {
        let free = |key: &str| having[key].iter().filter(|i| of[**i].is_none()).count();
        let untaken = |key: &str| !taken.iter().any(|(k, _, _)| k == key);
        // a preferred term first; else the most items not yet taken (ties:
        // the more specific, fewer in all, then by name)
        let first = preferred
            .iter()
            .filter_map(|p| having.get_key_value(p.as_str()).map(|(k, list)| (*k, list.len())))
            .find(|(key, _)| untaken(key) && free(key) >= min.max(1))
            .map(|(key, all)| (free(key), all, key));
        let best = first.or_else(|| {
            having
                .iter()
                .filter(|(key, _)| untaken(key))
                .map(|(key, list)| (free(key), list.len(), *key))
                .max_by(|x, y| x.0.cmp(&y.0).then(y.1.cmp(&x.1)).then(names[y.2].cmp(names[x.2])))
        });
        let Some((count, _, key)) = best else { break };
        if count < min.max(1) {
            break;
        }
        for &i in &having[key] {
            if of[i].is_none() {
                of[i] = Some(taken.len());
            }
        }
        taken.push((key.to_string(), names[key].to_string(), count));
    }
    let mut rest: Vec<(String, String, usize)> = having
        .iter()
        .filter(|(key, _)| !taken.iter().any(|(k, _, _)| k == *key))
        .map(|(key, list)| (key.to_string(), names[key].to_string(), list.len()))
        .collect();
    rest.sort_by(|x, y| y.2.cmp(&x.2).then(x.1.cmp(&y.1)));
    (taken, of, rest)
}

/// Communities by Louvain's method (Blondel et al., 2008): each node moved to
/// the neighbouring community that raises the modularity most, until none
/// moves; then the communities taken as nodes, and again, until nothing
/// changes. `edges` are `(node, node, weight)`. Returns each node's
/// community (numbered by size, largest first) and the modularity.
pub fn louvain(n: usize, edges: &[(usize, usize, f64)]) -> (Vec<usize>, f64) {
    if n == 0 {
        return (Vec::new(), 0.0);
    }
    let total: f64 = edges.iter().map(|e| e.2).sum();
    if total <= 0.0 {
        return ((0..n).collect(), 0.0);
    }
    // the graph of the current level: adjacency with weights (self-loops kept)
    let mut member: Vec<usize> = (0..n).collect();
    let mut adjacency: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for &(a, b, w) in edges {
        adjacency[a].push((b, w));
        if a != b {
            adjacency[b].push((a, w));
        }
    }
    let m2 = 2.0 * total;
    loop {
        let size = adjacency.len();
        let degree: Vec<f64> = adjacency
            .iter()
            .enumerate()
            .map(|(i, list)| list.iter().map(|&(j, w)| if j == i { 2.0 * w } else { w }).sum())
            .collect();
        let mut community: Vec<usize> = (0..size).collect();
        let mut tot: Vec<f64> = degree.clone();
        let mut moved_any = false;
        for _pass in 0..100 {
            let mut moved = false;
            for i in 0..size {
                let own = community[i];
                tot[own] -= degree[i];
                let mut links: BTreeMap<usize, f64> = BTreeMap::new();
                for &(j, w) in &adjacency[i] {
                    if j != i {
                        *links.entry(community[j]).or_default() += w;
                    }
                }
                let gain = |c: usize, w_ic: f64| w_ic - tot[c] * degree[i] / m2;
                let mut best = (own, gain(own, links.get(&own).copied().unwrap_or(0.0)));
                for (&c, &w) in &links {
                    let g = gain(c, w);
                    if g > best.1 + 1e-12 {
                        best = (c, g);
                    }
                }
                community[i] = best.0;
                tot[best.0] += degree[i];
                if best.0 != own {
                    moved = true;
                    moved_any = true;
                }
            }
            if !moved {
                break;
            }
        }
        // renumber, and merge each community into a node
        let mut renumber: BTreeMap<usize, usize> = BTreeMap::new();
        for &c in &community {
            let next = renumber.len();
            renumber.entry(c).or_insert(next);
        }
        let count = renumber.len();
        for m in member.iter_mut() {
            *m = renumber[&community[*m]];
        }
        if !moved_any || count == size {
            break;
        }
        let mut merged: Vec<BTreeMap<usize, f64>> = vec![BTreeMap::new(); count];
        for (i, list) in adjacency.iter().enumerate() {
            let ci = renumber[&community[i]];
            for &(j, w) in list {
                let cj = renumber[&community[j]];
                // each edge is listed at both ends (a self-loop once)
                let w = if i == j { w } else { w / 2.0 };
                if ci <= cj {
                    *merged[ci].entry(cj).or_default() += w;
                } else {
                    *merged[cj].entry(ci).or_default() += w;
                }
            }
        }
        let mut next: Vec<Vec<(usize, f64)>> = vec![Vec::new(); count];
        for (a, list) in merged.iter().enumerate() {
            for (&b, &w) in list {
                next[a].push((b, w));
                if a != b {
                    next[b].push((a, w));
                }
            }
        }
        adjacency = next;
    }
    // numbered by size, largest first (then by first member)
    let mut sizes: BTreeMap<usize, (usize, usize)> = BTreeMap::new();
    for (i, &c) in member.iter().enumerate() {
        let e = sizes.entry(c).or_insert((0, i));
        e.0 += 1;
    }
    let mut order: Vec<(usize, (usize, usize))> = sizes.into_iter().collect();
    order.sort_by(|x, y| y.1 .0.cmp(&x.1 .0).then(x.1 .1.cmp(&y.1 .1)));
    let rank: HashMap<usize, usize> = order.iter().enumerate().map(|(r, (c, _))| (*c, r)).collect();
    let member: Vec<usize> = member.iter().map(|c| rank[c]).collect();
    (member.clone(), modularity(&member, edges))
}

/// The modularity of a division of a graph into communities
pub fn modularity(community: &[usize], edges: &[(usize, usize, f64)]) -> f64 {
    let total: f64 = edges.iter().map(|e| e.2).sum();
    if total <= 0.0 {
        return 0.0;
    }
    let count = community.iter().max().map_or(0, |m| m + 1);
    let (mut inside, mut degree) = (vec![0.0; count], vec![0.0; count]);
    for &(a, b, w) in edges {
        if community[a] == community[b] {
            inside[community[a]] += w;
        }
        degree[community[a]] += w;
        degree[community[b]] += w;
    }
    (0..count).map(|c| inside[c] / total - (degree[c] / (2.0 * total)).powi(2)).sum()
}

// ---------------------------------------------------------------- the plan: what is drawn where

/// The nodes and edges drawn, their groups and their places
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Plan {
    /// The nodes drawn: indices into the data's nodes
    pub nodes: Vec<usize>,
    /// The edges drawn (both ends drawn)
    pub edges: Vec<usize>,
    /// The groups coloured, largest first
    pub groups: Vec<Group>,
    /// Each drawn node's group
    pub group_of: Vec<Option<usize>>,
    /// Terms of the grouping not coloured (key, name, proteins having it)
    pub others: Vec<(String, String, usize)>,
    /// The modularity of the network's modules (grouped by module)
    pub modularity: Option<f64>,
    /// Each drawn node's place as laid out, 0 to 1 across and down the box
    /// it was laid out in
    pub layout: Vec<[f64; 2]>,
    /// That box's width over its height (the drawing area's, when laid out)
    pub aspect: f64,
}

/// A node's label as the style asks
pub fn node_label(node: &Node, info: Option<&ProteinInfo>, style: &NetworkStyle) -> String {
    match style.labels.get(&node.accession).filter(|l| !l.trim().is_empty()) {
        Some(label) => label.clone(),
        None => default_label(node, info, style.label_by),
    }
}

/// A node's label of a kind (not one's own)
pub fn default_label(node: &Node, info: Option<&ProteinInfo>, by: NodeLabel) -> String {
    let first = |options: &[&str]| options.iter().find(|s| !s.is_empty()).map(|s| s.to_string()).unwrap_or_default();
    match by {
        NodeLabel::Gene => first(&[info.map_or("", |p| p.gene.as_str()), &node.gene, info.map_or("", |p| p.entry.as_str()), &node.entry]),
        NodeLabel::Entry => first(&[info.map_or("", |p| p.entry.as_str()), &node.entry]),
        NodeLabel::Accession => first(&[info.map_or("", |p| p.accession.as_str()), &node.accession]),
        NodeLabel::None => String::new(),
    }
}

/// Which nodes and edges are drawn, the groups and the layout
pub fn plan(data: &NetworkData, ann: &Annotations, own: &OwnGroups, style: &NetworkStyle) -> Plan {
    plan_with(data, ann, own, style, None)
}

/// [`plan`], with the layout (and its proportions) of an earlier plan for
/// the same nodes (when only the groups' labels or colours changed)
pub fn plan_with(data: &NetworkData, ann: &Annotations, own: &OwnGroups, style: &NetworkStyle, layout: Option<(Vec<[f64; 2]>, f64)>) -> Plan {
    let nodes: Vec<usize> = (0..data.nodes.len())
        .filter(|&i| {
            let n = &data.nodes[i];
            !style.hidden.contains(&n.accession) && (n.partners > 0 || style.lone)
        })
        .collect();
    let index: HashMap<&str, usize> = nodes.iter().enumerate().map(|(k, &i)| (data.nodes[i].accession.as_str(), k)).collect();
    let edges: Vec<usize> =
        (0..data.edges.len()).filter(|&e| index.contains_key(data.edges[e].a.as_str()) && index.contains_key(data.edges[e].b.as_str())).collect();
    let links: Vec<(usize, usize, f64)> =
        edges.iter().map(|&e| (index[data.edges[e].a.as_str()], index[data.edges[e].b.as_str()], data.edges[e].crosslinks as f64)).collect();

    // the groups
    let min = style.min_group.max(1);
    let max = style.max_groups.clamp(1, 60);
    let excluded = |t: &Term| style.excluded.iter().any(|x| x == t.key());
    let (mut taken, of, others, modularity) = match style.group_by {
        GroupBy::None => (Vec::new(), vec![None; nodes.len()], Vec::new(), None),
        GroupBy::Module => {
            let (community, q) = louvain(nodes.len(), &links);
            let items: Vec<Vec<Term>> = community
                .iter()
                .map(|&c| {
                    // a module's key: its first member, which outlives renumbering
                    let first = community.iter().position(|&d| d == c).unwrap_or(0);
                    vec![Term { id: format!("module:{}", data.nodes[nodes[first]].accession), name: format!("Module {}", c + 1) }]
                })
                .collect();
            let items: Vec<Vec<Term>> = items.into_iter().map(|t| t.into_iter().filter(|t| !excluded(t)).collect()).collect();
            let (taken, of, others) = cover(&items, min, max, &style.preferred);
            (taken, of, others, Some(q))
        }
        GroupBy::Own => {
            let items: Vec<Vec<Term>> = nodes
                .iter()
                .map(|&i| {
                    let n = &data.nodes[i];
                    let info = ann.proteins.get(&n.accession);
                    own.group(n, info, &node_label(n, info, style))
                        .map(|g| vec![Term { id: String::new(), name: g.clone() }])
                        .unwrap_or_default()
                        .into_iter()
                        .filter(|t| !excluded(t))
                        .collect()
                })
                .collect();
            let (taken, of, others) = cover(&items, min, max, &style.preferred);
            (taken, of, others, None)
        }
        by => {
            let items: Vec<Vec<Term>> = nodes
                .iter()
                .map(|&i| ann.proteins.get(&data.nodes[i].accession).map(|p| by.terms(p)).unwrap_or_default().into_iter().filter(|t| !excluded(t)).collect())
                .collect();
            let (taken, of, others) = cover(&items, min, max, &style.preferred);
            (taken, of, others, None)
        }
    };
    // modules: renamed in order
    if style.group_by == GroupBy::Module {
        for (i, t) in taken.iter_mut().enumerate() {
            t.1 = format!("Module {}", i + 1);
        }
    }
    let groups: Vec<Group> = taken
        .into_iter()
        .enumerate()
        .map(|(i, (key, name, members))| {
            let label = style.group_labels.get(&key).filter(|l| !l.trim().is_empty()).cloned().unwrap_or_else(|| name.clone());
            let color = style
                .group_colors
                .get(&key)
                .or_else(|| own.colors.get(&key))
                .cloned()
                .unwrap_or_else(|| style.palette.color(i).to_string());
            Group { key, name, label, color, members }
        })
        .collect();

    let mut plan = Plan { nodes, edges, groups, group_of: of, others, modularity, layout: Vec::new(), aspect: 1.0 };
    match layout {
        Some((layout, aspect)) if layout.len() == plan.nodes.len() => {
            plan.layout = layout;
            plan.aspect = aspect;
        }
        _ => {
            // laid out in the proportions of the room the legend and labels leave
            plan.layout = vec![[0.5, 0.5]; plan.nodes.len()];
            let f = frame(data, ann, &plan, style);
            plan.aspect = ((f.area[2] - f.area[0]) / (f.area[3] - f.area[1]).max(1e-6)).clamp(0.2, 5.0);
            plan.layout = lay_out(plan.nodes.len(), &links, &plan.group_of, plan.groups.len(), style, plan.aspect);
        }
    }
    plan
}

// ---------------------------------------------------------------- the layout

/// A small, seeded random number generator (xorshift)
struct Random(u64);

impl Random {
    fn new(seed: u32) -> Random {
        Random(0x9e37_79b9_7f4a_7c15 ^ (seed as u64 + 1).wrapping_mul(0xbf58_476d_1ce4_e5b9))
    }

    /// 0 to 1
    fn next(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// Places for `n` nodes joined by `links`, 0 to 1 across and down a box of
/// proportions `aspect` (width over height): force-directed (Fruchterman &
/// Reingold, 1991), with the nodes of a group drawn together, or on a circle
fn lay_out(n: usize, links: &[(usize, usize, f64)], group_of: &[Option<usize>], groups: usize, style: &NetworkStyle, aspect: f64) -> Vec<[f64; 2]> {
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![[0.5, 0.5]];
    }
    let aspect = aspect.clamp(0.2, 5.0);
    if style.layout == Layout::Circle {
        // group by group (ungrouped last), the most linked first
        let mut degree = vec![0.0; n];
        for &(a, b, w) in links {
            degree[a] += w;
            degree[b] += w;
        }
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&x, &y| {
            (group_of[x].unwrap_or(usize::MAX), -degree[x], x).partial_cmp(&(group_of[y].unwrap_or(usize::MAX), -degree[y], y)).unwrap_or(std::cmp::Ordering::Equal)
        });
        // a circle, not an ellipse, in the box
        let radius = 0.5 * aspect.min(1.0);
        let mut out = vec![[0.5, 0.5]; n];
        for (k, &i) in order.iter().enumerate() {
            let angle = std::f64::consts::TAU * k as f64 / n as f64 - std::f64::consts::FRAC_PI_2;
            out[i] = [0.5 + radius * angle.cos() / aspect, 0.5 + radius * angle.sin()];
        }
        return out;
    }

    let mut random = Random::new(style.seed);
    let centre = [aspect / 2.0, 0.5];
    let pull = match style.layout {
        Layout::Groups => style.group_pull.max(0.0) + 1.0,
        _ => style.group_pull.max(0.0),
    };
    // the start: at random; clusters: each group around a point of a circle
    let mut pos: Vec<[f64; 2]> = (0..n)
        .map(|i| {
            let jitter = [random.next() - 0.5, random.next() - 0.5];
            match (style.layout, group_of[i]) {
                (Layout::Groups, Some(g)) => {
                    let angle = std::f64::consts::TAU * g as f64 / groups.max(1) as f64;
                    [centre[0] + 0.35 * aspect * angle.cos() + jitter[0] * 0.1, centre[1] + 0.35 * angle.sin() + jitter[1] * 0.1]
                }
                _ => [centre[0] + jitter[0] * aspect, centre[1] + jitter[1]],
            }
        })
        .collect();
    let k = 0.75 * style.spacing.clamp(0.2, 5.0) * (aspect / n as f64).sqrt();
    let iterations = if n <= 300 { 500 } else { (150_000 / n).clamp(60, 500) };
    let mut temperature = 0.12 * aspect.max(1.0);
    let cooling = temperature / iterations as f64;
    let gravity = 0.6;
    let mut disp = vec![[0.0f64; 2]; n];
    for _ in 0..iterations {
        for d in disp.iter_mut() {
            *d = [0.0, 0.0];
        }
        // every pair pushes apart (more so across groups, in clusters)
        for i in 0..n {
            for j in (i + 1)..n {
                let (dx, dy) = (pos[i][0] - pos[j][0], pos[i][1] - pos[j][1]);
                let dist = dx.hypot(dy).max(k * 0.01);
                let mut f = k * k / dist;
                if style.layout == Layout::Groups && group_of[i] != group_of[j] {
                    f *= 2.0;
                }
                let (fx, fy) = (dx / dist * f, dy / dist * f);
                disp[i][0] += fx;
                disp[i][1] += fy;
                disp[j][0] -= fx;
                disp[j][1] -= fy;
            }
        }
        // edges pull together, stronger the more crosslinks
        for &(a, b, w) in links {
            let (dx, dy) = (pos[a][0] - pos[b][0], pos[a][1] - pos[b][1]);
            let dist = dx.hypot(dy).max(1e-9);
            let f = dist * dist / k * (1.0 + w.ln().max(0.0) * 0.3);
            let (fx, fy) = (dx / dist * f, dy / dist * f);
            disp[a][0] -= fx;
            disp[a][1] -= fy;
            disp[b][0] += fx;
            disp[b][1] += fy;
        }
        // groups gather around their middle
        if pull > 0.0 && groups > 0 {
            let mut middle = vec![[0.0f64; 3]; groups];
            for i in 0..n {
                if let Some(g) = group_of[i] {
                    middle[g][0] += pos[i][0];
                    middle[g][1] += pos[i][1];
                    middle[g][2] += 1.0;
                }
            }
            for i in 0..n {
                if let Some(g) = group_of[i] {
                    let m = middle[g];
                    let (dx, dy) = (m[0] / m[2] - pos[i][0], m[1] / m[2] - pos[i][1]);
                    let dist = dx.hypot(dy);
                    disp[i][0] += dx * dist / k * pull;
                    disp[i][1] += dy * dist / k * pull;
                }
            }
        }
        // and everything to the centre, so that unconnected parts stay near:
        // more across than down in a tall box, and the other way in a wide
        // one, for the network to take the box's shape
        for i in 0..n {
            disp[i][0] += (centre[0] - pos[i][0]) * gravity / aspect;
            disp[i][1] += (centre[1] - pos[i][1]) * gravity * aspect;
        }
        for i in 0..n {
            let len = disp[i][0].hypot(disp[i][1]);
            if len > 0.0 {
                let step = len.min(temperature);
                pos[i][0] += disp[i][0] / len * step;
                pos[i][1] += disp[i][1] / len * step;
            }
        }
        temperature = (temperature - cooling).max(temperature * 0.02);
    }
    // fitted to the area, its proportions kept
    let (mut x0, mut y0, mut x1, mut y1) = (f64::MAX, f64::MAX, f64::MIN, f64::MIN);
    for p in &pos {
        x0 = x0.min(p[0]);
        x1 = x1.max(p[0]);
        y0 = y0.min(p[1]);
        y1 = y1.max(p[1]);
    }
    let scale = ((x1 - x0) / aspect).max(y1 - y0).max(1e-9);
    let (cx, cy) = ((x0 + x1) / 2.0, (y0 + y1) / 2.0);
    pos.iter().map(|p| [0.5 + (p[0] - cx) / (scale * aspect), 0.5 + (p[1] - cy) / scale]).collect()
}

// ---------------------------------------------------------------- the figure

/// SVG user units per inch (as the other figures)
const UNITS_PER_IN: f64 = 150.0;

fn pt(points: f64) -> f64 {
    points / 72.0 * UNITS_PER_IN
}

fn esc(s: &str) -> String {
    s.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;").replace('"', "&quot;")
}

/// `#rrggbb` mixed with `#rrggbb`, `t` of the way
fn mix(a: &str, b: &str, t: f64) -> String {
    let (a, b) = (circos::rgb(a), circos::rgb(b));
    let t = t.clamp(0.0, 1.0);
    let m = |x: u8, y: u8| (x as f64 + (y as f64 - x as f64) * t).round() as u8;
    format!("#{:02x}{:02x}{:02x}", m(a.0, b.0), m(a.1, b.1), m(a.2, b.2))
}

/// Text's width, roughly (average character widths)
fn text_width(text: &str, t: &TypeStyle) -> f64 {
    text.chars().count() as f64 * pt(t.size.clamp(2.0, 96.0)) * if t.bold { 0.62 } else { 0.56 }
}

/// An edge's colour, and what that colour means (its legend row)
#[derive(Clone, Debug, PartialEq)]
pub struct EdgeLook {
    pub color: String,
    /// Index into the legend's rows of edge colours
    pub row: usize,
    /// STRING's score (the evidence chosen), if asked
    pub string: Option<f32>,
    /// Known to STRING at the style's evidence and score; `None`: not looked up
    pub known: Option<bool>,
    pub width: f64,
}

/// Where everything goes: the canvas, the drawing area, each node and the legend
#[derive(Clone, Debug, PartialEq)]
pub struct Frame {
    pub width: f64,
    pub height: f64,
    /// Where node centres may go: x0, y0, x1, y1
    pub area: [f64; 4],
    /// The box the layout fills (its proportions kept, in the middle of the area)
    pub layout_box: [f64; 4],
    /// The network's room: the canvas but for the title, legend and margins
    pub region: [f64; 4],
    /// Each drawn node: centre and radius
    pub nodes: Vec<(f64, f64, f64)>,
    pub labels: Vec<String>,
    pub colors: Vec<String>,
    pub edges: Vec<EdgeLook>,
    /// The legend's rows of edge colours: colour, text
    pub edge_rows: Vec<(String, String)>,
    /// The legend's blocks (its sections), each at its top left corner:
    /// one below the other at the side, side by side below the network
    pub legend: Vec<(f64, f64, Vec<Row>)>,
    pub title_space: f64,
}

/// A row of the legend
#[derive(Clone, Debug, PartialEq)]
pub enum Row {
    Heading(String),
    Dot { color: String, outline: String, text: String },
    Line { color: String, width: f64, text: String },
    Loop { color: String, text: String },
    Size { diameter: f64, text: String },
    Text(String),
    Gap,
}

/// What makes an interaction known, for the legend
fn known_text(style: &NetworkStyle) -> String {
    format!("known: STRING score ≥ {:.2}, {}", style.min_score, style.evidence.label())
}

/// `text` cut to `width` (with an ellipsis), as far as its width can be told
fn fit(text: &str, width: f64, t: &TypeStyle) -> String {
    if text_width(text, t) <= width {
        return text.to_string();
    }
    let per_char = text_width("x", t).max(1e-9);
    let keep = ((width / per_char) as usize).saturating_sub(1).max(4);
    let mut out: String = text.chars().take(keep).collect();
    out.push('…');
    out
}

/// The colour of each drawn edge, and the legend rows for them
fn edge_looks(data: &NetworkData, ann: &Annotations, plan: &Plan, style: &NetworkStyle) -> (Vec<EdgeLook>, Vec<(String, String)>) {
    let edges: Vec<&Edge> = plan.edges.iter().map(|&e| &data.edges[e]).collect();
    // widths
    let value = |e: &Edge| match style.width_by {
        EdgeWidth::Crosslinks => e.crosslinks as f64,
        EdgeWidth::Csms => e.csms as f64,
        EdgeWidth::Fixed => 1.0,
    };
    let most = edges.iter().map(|e| value(e)).fold(1.0, f64::max);
    let (thin, wide) = (pt(style.edge_width.clamp(0.05, 30.0)), pt(style.edge_width_max.clamp(0.05, 30.0)));
    let width = |v: f64| if most > 1.0 { thin + (wide - thin) * v.max(1.0).ln() / most.ln() } else { thin };

    let mut rows: Vec<(String, String)> = Vec::new();
    let row = |rows: &mut Vec<(String, String)>, color: &str, text: String| -> usize {
        match rows.iter().position(|(c, t)| c == color && *t == text) {
            Some(i) => i,
            None => {
                rows.push((color.to_string(), text));
                rows.len() - 1
            }
        }
    };
    let string_of = |e: &Edge| ann.interaction(&e.a, &e.b).map(|i| i.map_or(0.0, |i| style.evidence.score(i)));
    let mut looks: Vec<EdgeLook> = Vec::new();
    // shared terms (complex, pathway): the group's colour where grouped by
    // the same, else colours of their own, the most shared first
    let shared = |by: GroupBy| -> Vec<Vec<Term>> {
        edges
            .iter()
            .map(|e| match (ann.proteins.get(&e.a), ann.proteins.get(&e.b)) {
                (Some(p), Some(q)) => {
                    let theirs: HashSet<String> = by.terms(q).iter().map(|t| t.key().to_string()).collect();
                    by.terms(p).into_iter().filter(|t| theirs.contains(t.key())).collect()
                }
                _ => Vec::new(),
            })
            .collect()
    };
    let colors: Vec<(String, String)> = match style.color_by {
        EdgeColor::Known => {
            let known_label = "known to STRING".to_string();
            edges
                .iter()
                .map(|e| match string_of(e) {
                    Some(s) if s >= style.min_score => (style.known_color.clone(), known_label.clone()),
                    Some(_) => (style.new_color.clone(), "new: not known to STRING".to_string()),
                    None => (style.new_color.clone(), "not looked up in STRING".to_string()),
                })
                .collect()
        }
        EdgeColor::Complex | EdgeColor::Pathway => {
            let by = if style.color_by == EdgeColor::Complex { GroupBy::Complex } else { GroupBy::Pathway };
            let terms = shared(by);
            let none = if by == GroupBy::Complex { "no complex shared" } else { "no pathway shared" };
            let (taken, of, _) = cover(&terms, 1, style.max_groups.clamp(1, 60), &style.preferred);
            terms
                .iter()
                .zip(&of)
                .map(|(mine, t)| {
                    // grouped by the same: a group's colour
                    if style.group_by == by {
                        if let Some(g) = plan.groups.iter().find(|g| mine.iter().any(|m| m.key() == g.key)) {
                            return (g.color.clone(), format!("within {}", g.label));
                        }
                    }
                    match t {
                        Some(i) => {
                            let (key, name, _) = &taken[*i];
                            let label = style.group_labels.get(key).cloned().unwrap_or_else(|| name.clone());
                            // (after the groups' colours, not to be taken for them)
                            let color = style.group_colors.get(key).cloned().unwrap_or_else(|| style.palette.color(plan.groups.len() + *i).to_string());
                            (color, format!("within {label}"))
                        }
                        None => (style.other_color.clone(), none.to_string()),
                    }
                })
                .collect()
        }
        EdgeColor::Group => {
            let index: HashMap<&str, usize> = plan.nodes.iter().enumerate().map(|(k, &i)| (data.nodes[i].accession.as_str(), k)).collect();
            edges
                .iter()
                .map(|e| match (plan.group_of[index[e.a.as_str()]], plan.group_of[index[e.b.as_str()]]) {
                    (Some(g), Some(h)) if g == h => (plan.groups[g].color.clone(), format!("within {}", plan.groups[g].label)),
                    _ => (style.other_color.clone(), "between groups".to_string()),
                })
                .collect()
        }
        EdgeColor::Single => edges.iter().map(|_| (style.edge_color.clone(), "crosslinked proteins".to_string())).collect(),
    };
    for (e, (color, text)) in edges.iter().zip(colors) {
        let i = row(&mut rows, &color, text);
        let string = string_of(e);
        looks.push(EdgeLook {
            color,
            row: i,
            string,
            known: string.map(|s| s >= style.min_score),
            width: width(value(e)),
        });
    }
    // the legend's rows with their counts, the coloured first
    let counts: Vec<usize> = (0..rows.len()).map(|r| looks.iter().filter(|l| l.row == r).count()).collect();
    let rows = rows.into_iter().zip(counts).map(|((c, t), n)| (c, format!("{t} ({n})"))).collect();
    (looks, rows)
}

/// Where everything goes
pub fn frame(data: &NetworkData, ann: &Annotations, plan: &Plan, style: &NetworkStyle) -> Frame {
    let (w, h) = (style.width.clamp(1.0, 40.0) * UNITS_PER_IN, style.height.clamp(1.0, 40.0) * UNITS_PER_IN);
    let pad = w.min(h) * 0.03;
    let title_space = if style.title.trim().is_empty() { 0.0 } else { pt(style.title_text.size) * 1.8 };
    let nodes: Vec<&Node> = plan.nodes.iter().map(|&i| &data.nodes[i]).collect();
    let info = |n: &Node| ann.proteins.get(&n.accession);

    // node sizes
    let value = |n: &Node| match style.size_by {
        NodeSize::Crosslinks => (n.crosslinks + n.self_links) as f64,
        NodeSize::Csms => (n.csms + n.self_csms) as f64,
        NodeSize::Partners => n.partners as f64,
        NodeSize::Fixed => 1.0,
    };
    let (lo, hi) = nodes.iter().map(|n| value(n)).fold((f64::MAX, f64::MIN), |(lo, hi), v| (lo.min(v), hi.max(v)));
    let (small, large) = (pt(style.node_size.clamp(1.0, 200.0)), pt(style.node_size_max.clamp(1.0, 200.0)).max(pt(style.node_size.clamp(1.0, 200.0))));
    let diameter = |v: f64| if hi > lo { small + (large - small) * ((v - lo) / (hi - lo)).sqrt() } else { small };
    let labels: Vec<String> = nodes.iter().map(|n| node_label(n, info(n), style)).collect();
    let colors: Vec<String> = plan
        .group_of
        .iter()
        .map(|g| g.map(|g| plan.groups[g].color.clone()).unwrap_or_else(|| style.ungrouped_color.clone()))
        .collect();
    let (edges, edge_rows) = edge_looks(data, ann, plan, style);

    // the legend
    let mut legend: Vec<Row> = Vec::new();
    if style.legend {
        if style.group_by != GroupBy::None && !plan.groups.is_empty() {
            legend.push(Row::Heading(style.group_by.heading().to_string()));
            for g in &plan.groups {
                legend.push(Row::Dot { color: g.color.clone(), outline: style.node_outline.clone(), text: format!("{} ({})", g.label, g.members) });
            }
            let ungrouped = plan.group_of.iter().filter(|g| g.is_none()).count();
            if ungrouped > 0 {
                legend.push(Row::Dot { color: style.ungrouped_color.clone(), outline: style.node_outline.clone(), text: format!("other ({ungrouped})") });
            }
            legend.push(Row::Gap);
        }
        if !edges.is_empty() {
            legend.push(Row::Heading("Crosslinked proteins".into()));
            for (color, text) in &edge_rows {
                legend.push(Row::Line { color: color.clone(), width: pt(1.6), text: text.clone() });
            }
            if style.color_by == EdgeColor::Known && edges.iter().any(|e| e.known.is_some()) {
                legend.push(Row::Text(known_text(style)));
            }
            let loops = nodes.iter().filter(|n| n.self_links > 0).count();
            if style.intra && loops > 0 {
                legend.push(Row::Loop { color: style.self_color.clone(), text: format!("crosslinks within a protein ({loops})") });
            }
            if style.width_by != EdgeWidth::Fixed {
                let value = |e: &Edge| if style.width_by == EdgeWidth::Crosslinks { e.crosslinks } else { e.csms as usize };
                let most = plan.edges.iter().map(|&e| value(&data.edges[e])).max().unwrap_or(1);
                let unit = if style.width_by == EdgeWidth::Crosslinks { "crosslink" } else { "CSM" };
                let mut samples = vec![1, ((most as f64).sqrt().round() as usize).max(1), most];
                samples.dedup();
                let (thin, wide) = (pt(style.edge_width.clamp(0.05, 30.0)), pt(style.edge_width_max.clamp(0.05, 30.0)));
                for v in samples {
                    let width = if most > 1 { thin + (wide - thin) * (v as f64).ln() / (most as f64).ln() } else { thin };
                    legend.push(Row::Line { color: style.new_color.clone(), width, text: format!("{v} {unit}{}", if v == 1 { "" } else { "s" }) });
                }
            }
            legend.push(Row::Gap);
        }
        if style.size_by != NodeSize::Fixed && hi > lo {
            legend.push(Row::Heading(format!("Node size: {}", style.size_by.label())));
            legend.push(Row::Size { diameter: diameter(lo), text: format!("{lo}") });
            legend.push(Row::Size { diameter: diameter(hi), text: format!("{hi}") });
            legend.push(Row::Gap);
        }
        if style.summary {
            let known = edges.iter().filter(|e| e.known == Some(true)).count();
            let looked = edges.iter().any(|e| e.known.is_some());
            let crosslinks: usize = plan.edges.iter().map(|&e| data.edges[e].crosslinks).sum();
            let csms: u32 = plan.edges.iter().map(|&e| data.edges[e].csms).sum();
            legend.push(Row::Text(format!("{} proteins, {} crosslinked pairs", nodes.len(), edges.len())));
            if looked {
                legend.push(Row::Text(format!("{known} of them known to STRING")));
            }
            legend.push(Row::Text(format!("{crosslinks} crosslinks, {csms} CSMs, q ≤ {}", style.max_q)));
            if let Some(q) = plan.modularity {
                legend.push(Row::Text(format!("modularity {q:.2}")));
            }
        }
        while legend.last() == Some(&Row::Gap) {
            legend.pop();
        }
    }
    // texts too long for the legend's room shortened
    let line = pt(style.legend_text.size) * 1.45;
    let room = if style.legend_side == LegendSide::Below { w - 2.0 * pad } else { w * 0.45 };
    for row in legend.iter_mut() {
        match row {
            Row::Heading(t) => *t = fit(t, room, &TypeStyle { bold: true, ..style.legend_text }),
            Row::Dot { text, .. } | Row::Line { text, .. } | Row::Loop { text, .. } => *text = fit(text, room - line * 1.8, &style.legend_text),
            Row::Text(t) => *t = fit(t, room, &style.legend_text),
            Row::Size { .. } | Row::Gap => {}
        }
    }
    let row_width = |r: &Row| match r {
        Row::Heading(t) => text_width(t, &TypeStyle { bold: true, ..style.legend_text }),
        Row::Dot { text, .. } | Row::Line { text, .. } | Row::Loop { text, .. } => line * 1.8 + text_width(text, &style.legend_text),
        Row::Size { diameter, text } => diameter.max(line) + line * 0.6 + text_width(text, &style.legend_text),
        Row::Text(t) => text_width(t, &style.legend_text),
        Row::Gap => 0.0,
    };
    let row_height = |r: &Row| match r {
        Row::Size { diameter, .. } => diameter.max(line) + line * 0.2,
        Row::Gap => line * 0.5,
        _ => line,
    };
    // its sections, and their sizes
    let mut blocks: Vec<Vec<Row>> = Vec::new();
    for row in legend {
        match (row, blocks.last_mut()) {
            (Row::Gap, _) => blocks.push(Vec::new()),
            (row, Some(block)) => block.push(row),
            (row, None) => blocks.push(vec![row]),
        }
    }
    blocks.retain(|b| !b.is_empty());
    let sizes: Vec<(f64, f64)> =
        blocks.iter().map(|b| (b.iter().map(row_width).fold(0.0, f64::max).min(room), b.iter().map(row_height).sum())).collect();
    let block_gap = line * 0.5;

    // the drawing area: the canvas but for the title, the legend and margins
    let (mut x0, y0, mut x1, mut y1) = (pad, pad + title_space, w - pad, h - pad);
    let gap = pad;
    let mut placed: Vec<(f64, f64, Vec<Row>)> = Vec::new();
    if !blocks.is_empty() {
        match style.legend_side {
            LegendSide::Right | LegendSide::Left => {
                // one below the other, the whole in the middle of the height
                let width = sizes.iter().map(|s| s.0).fold(0.0, f64::max);
                let height = sizes.iter().map(|s| s.1).sum::<f64>() + block_gap * (sizes.len() - 1) as f64;
                let left = if style.legend_side == LegendSide::Right { x1 - width } else { x0 };
                let mut y = y0 + ((y1 - y0) - height).max(0.0) / 2.0;
                for (block, (_, h)) in blocks.into_iter().zip(&sizes) {
                    placed.push((left, y, block));
                    y += h + block_gap;
                }
                if style.legend_side == LegendSide::Right {
                    x1 -= width + gap;
                } else {
                    x0 += width + gap;
                }
            }
            LegendSide::Below => {
                // side by side, onto a new line when the width is full
                let across = line * 1.5;
                let mut lines: Vec<Vec<usize>> = vec![Vec::new()];
                let mut used = 0.0;
                for (i, (bw, _)) in sizes.iter().enumerate() {
                    let current = lines.last_mut().expect("a line");
                    if !current.is_empty() && used + across + bw > x1 - x0 {
                        lines.push(vec![i]);
                        used = *bw;
                    } else {
                        used += if current.is_empty() { *bw } else { across + bw };
                        current.push(i);
                    }
                }
                let heights: Vec<f64> = lines.iter().map(|l| l.iter().map(|&i| sizes[i].1).fold(0.0, f64::max)).collect();
                let total = heights.iter().sum::<f64>() + block_gap * (heights.len() - 1) as f64;
                let mut y = y1 - total;
                let mut blocks: Vec<Option<Vec<Row>>> = blocks.into_iter().map(Some).collect();
                for (l, height) in lines.iter().zip(&heights) {
                    let mut x = x0;
                    for &i in l {
                        placed.push((x, y, blocks[i].take().unwrap_or_default()));
                        x += sizes[i].0 + across;
                    }
                    y += height + block_gap;
                }
                y1 -= total + gap;
            }
        }
    }
    let legend = placed;
    // node centres keep clear of the area's edges by the largest node and its label
    let label_size = pt(style.label_text.size);
    let widest_label = labels.iter().map(|l| text_width(l, &style.label_text)).fold(0.0, f64::max);
    let r_max = large / 2.0;
    let (mx_left, mx_right, my_top, my_bottom) = match style.label_place {
        LabelPlace::Center => (r_max, r_max, r_max, r_max),
        LabelPlace::Beside => (r_max, r_max + widest_label.min((x1 - x0) * 0.25) + label_size * 0.4, r_max.max(label_size * 0.6), r_max.max(label_size * 0.6)),
        LabelPlace::Below => (r_max.max(widest_label / 2.0), r_max.max(widest_label / 2.0), r_max, r_max + label_size * 1.3),
    };
    let group_room = if style.group_names && !plan.groups.is_empty() { pt(style.group_text.size) * 1.4 } else { 0.0 };
    let area = [x0 + mx_left + group_room * 0.5, y0 + my_top + group_room, x1 - mx_right - group_room * 0.5, y1 - my_bottom - group_room * 0.5];
    let area = [area[0], area[1], area[2].max(area[0] + 1.0), area[3].max(area[1] + 1.0)];
    let (aw, ah) = (area[2] - area[0], area[3] - area[1]);
    let aspect = plan.aspect.clamp(0.2, 5.0);
    let (bw, bh) = if aw / ah > aspect { (ah * aspect, ah) } else { (aw, aw / aspect) };
    let layout_box = [area[0] + (aw - bw) / 2.0, area[1] + (ah - bh) / 2.0, area[0] + (aw + bw) / 2.0, area[1] + (ah + bh) / 2.0];
    // (a node moved by hand may go anywhere in the area)
    let place = |p: [f64; 2]| {
        (
            (layout_box[0] + p[0] * bw).clamp(area[0], area[2]),
            (layout_box[1] + p[1] * bh).clamp(area[1], area[3]),
        )
    };
    let node_places: Vec<(f64, f64, f64)> = nodes
        .iter()
        .enumerate()
        .map(|(k, n)| {
            let p = style.positions.get(&n.accession).copied().or_else(|| plan.layout.get(k).copied()).unwrap_or([0.5, 0.5]);
            let (x, y) = place(p);
            (x, y, diameter(value(n)) / 2.0)
        })
        .collect();
    Frame { width: w, height: h, area, layout_box, region: [x0, y0, x1, y1], nodes: node_places, labels, colors, edges, edge_rows, legend, title_space }
}

/// The network as SVG
pub fn svg(data: &NetworkData, ann: &Annotations, plan: &Plan, style: &NetworkStyle) -> String {
    let f = frame(data, ann, plan, style);
    let (w, h) = (f.width, f.height);
    let ink = if style.background == Background::Black { "#e8e8e8" } else { "#222222" };
    let font = esc(&style.font);
    let mut out = String::new();
    let _ = write!(
        out,
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{:.2}in" height="{:.2}in" viewBox="0 0 {w:.1} {h:.1}">"#,
        w / UNITS_PER_IN,
        h / UNITS_PER_IN
    );
    match style.background {
        Background::White => {
            let _ = write!(out, r##"<rect width="{w:.1}" height="{h:.1}" fill="#ffffff"/>"##);
        }
        Background::Black => {
            let _ = write!(out, r##"<rect width="{w:.1}" height="{h:.1}" fill="#000000"/>"##);
        }
        Background::Transparent => {}
    }
    if f.title_space > 0.0 {
        let _ = write!(
            out,
            r#"<text x="{:.1}" y="{:.1}" text-anchor="middle" {} fill="{ink}">{}</text>"#,
            w / 2.0,
            f.title_space * 0.8,
            style.title_text.attrs(&font),
            esc(&style.title)
        );
    }
    if plan.nodes.is_empty() {
        let _ = write!(
            out,
            r#"<text x="{:.1}" y="{:.1}" text-anchor="middle" {} fill="{ink}">No crosslinked proteins at q ≤ {}</text></svg>"#,
            w / 2.0,
            h / 2.0,
            style.label_text.attrs(&font),
            style.max_q
        );
        return out;
    }
    let index: HashMap<&str, usize> = plan.nodes.iter().enumerate().map(|(k, &i)| (data.nodes[i].accession.as_str(), k)).collect();

    // edges: the coloured on top of the grey, the widest last
    let mut order: Vec<usize> = (0..plan.edges.len()).collect();
    let plain = |c: &str| c == style.new_color || c == style.other_color;
    order.sort_by(|&x, &y| {
        (!plain(&f.edges[x].color), f.edges[x].width).partial_cmp(&(!plain(&f.edges[y].color), f.edges[y].width)).unwrap_or(std::cmp::Ordering::Equal)
    });
    let opacity = style.edge_opacity.clamp(0.05, 1.0);
    let _ = write!(out, r#"<g fill="none" stroke-linecap="round" stroke-opacity="{opacity:.2}">"#);
    for &k in &order {
        let e = &data.edges[plan.edges[k]];
        let (a, b) = (f.nodes[index[e.a.as_str()]], f.nodes[index[e.b.as_str()]]);
        let look = &f.edges[k];
        let _ = write!(
            out,
            r#"<line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}" stroke="{}" stroke-width="{:.2}"/>"#,
            a.0, a.1, b.0, b.1, look.color, look.width
        );
    }
    out.push_str("</g>");
    // crosslinks within a protein: a loop on its node
    if style.intra {
        let _ = write!(out, r#"<g fill="none" stroke="{}" stroke-width="{:.2}">"#, style.self_color, pt(style.edge_width.clamp(0.05, 30.0)).max(pt(0.6)));
        for (k, &i) in plan.nodes.iter().enumerate() {
            if data.nodes[i].self_links > 0 {
                let (x, y, r) = f.nodes[k];
                let lr = r * 0.55;
                let _ = write!(out, r#"<circle cx="{:.2}" cy="{:.2}" r="{lr:.2}"/>"#, x + r * 0.72, y - r * 0.72);
            }
        }
        out.push_str("</g>");
    }
    // the number of crosslinks on each edge
    if style.numbers {
        let size = pt(style.number_text.size);
        let back = if style.background == Background::Black { "#000000" } else { "#ffffff" };
        let _ = write!(out, r#"<g {} text-anchor="middle" fill="{ink}">"#, style.number_text.attrs(&font));
        for &e in &plan.edges {
            let e = &data.edges[e];
            let (a, b) = (f.nodes[index[e.a.as_str()]], f.nodes[index[e.b.as_str()]]);
            let (x, y) = ((a.0 + b.0) / 2.0, (a.1 + b.1) / 2.0);
            let text = e.crosslinks.to_string();
            let half = (text_width(&text, &style.number_text) / 2.0).max(size * 0.35) + size * 0.2;
            let _ = write!(
                out,
                r#"<rect x="{:.2}" y="{:.2}" width="{:.2}" height="{:.2}" rx="{:.2}" fill="{back}" fill-opacity="0.85"/><text x="{x:.2}" y="{:.2}">{text}</text>"#,
                x - half,
                y - size * 0.6,
                half * 2.0,
                size * 1.2,
                size * 0.3,
                y + size * 0.35
            );
        }
        out.push_str("</g>");
    }
    // nodes
    let _ = write!(out, r#"<g stroke="{}" stroke-width="{:.2}">"#, style.node_outline, pt(0.7));
    for (k, (x, y, r)) in f.nodes.iter().enumerate() {
        let _ = write!(out, r#"<circle cx="{x:.2}" cy="{y:.2}" r="{r:.2}" fill="{}"/>"#, f.colors[k]);
    }
    out.push_str("</g>");
    // labels
    let size = pt(style.label_text.size);
    let _ = write!(out, r#"<g {} fill="{ink}">"#, style.label_text.attrs(&font));
    for (k, (x, y, r)) in f.nodes.iter().enumerate() {
        let label = &f.labels[k];
        if label.is_empty() {
            continue;
        }
        let (tx, ty, anchor) = match style.label_place {
            LabelPlace::Center => (*x, y + size * 0.35, "middle"),
            LabelPlace::Beside => (x + r + size * 0.25, y + size * 0.35, "start"),
            LabelPlace::Below => (*x, y + r + size * 1.0, "middle"),
        };
        let _ = write!(out, r#"<text x="{tx:.2}" y="{ty:.2}" text-anchor="{anchor}">{}</text>"#, esc(label));
    }
    out.push_str("</g>");
    // group names: above, below or beside each group (with its labels),
    // where neither nodes nor other names are
    if style.group_names && !plan.groups.is_empty() {
        let gsize = pt(style.group_text.size);
        let lsize = pt(style.label_text.size);
        let [rx0, ry0, rx1, ry1] = f.region;
        let mut taken: Vec<[f64; 4]> = f.nodes.iter().map(|&(x, y, r)| [x - r, y - r, x + r, y + r]).collect();
        let overlaps = |a: &[f64; 4], b: &[f64; 4]| a[0] < b[2] && a[2] > b[0] && a[1] < b[3] && a[3] > b[1];
        let _ = write!(out, r#"<g {} text-anchor="middle">"#, style.group_text.attrs(&font));
        for (g, group) in plan.groups.iter().enumerate() {
            let mut bbox = [f64::MAX, f64::MAX, f64::MIN, f64::MIN];
            for (k, &(x, y, r)) in f.nodes.iter().enumerate() {
                if plan.group_of[k] != Some(g) {
                    continue;
                }
                let tw = text_width(&f.labels[k], &style.label_text);
                let b = match style.label_place {
                    LabelPlace::Center => [x - r.max(tw / 2.0), y - r, x + r.max(tw / 2.0), y + r],
                    LabelPlace::Beside => [x - r, y - r.max(lsize * 0.6), x + r + lsize * 0.25 + tw, y + r.max(lsize * 0.6)],
                    LabelPlace::Below => [x - r.max(tw / 2.0), y - r, x + r.max(tw / 2.0), y + r + lsize * 1.2],
                };
                bbox = [bbox[0].min(b[0]), bbox[1].min(b[1]), bbox[2].max(b[2]), bbox[3].max(b[3])];
            }
            if bbox[0] > bbox[2] {
                continue;
            }
            let tw = text_width(&group.label, &style.group_text);
            let (cx, cy) = ((bbox[0] + bbox[2]) / 2.0, (bbox[1] + bbox[3]) / 2.0);
            let gap = gsize * 0.3;
            // where the name's middle and baseline may go
            let spots = [
                (cx, bbox[1] - gap),
                (cx, bbox[3] + gap + gsize * 0.8),
                (bbox[0] - gap - tw / 2.0, cy + gsize * 0.35),
                (bbox[2] + gap + tw / 2.0, cy + gsize * 0.35),
            ];
            let rect = |(x, y): (f64, f64)| [x - tw / 2.0, y - gsize * 0.8, x + tw / 2.0, y + gsize * 0.2];
            let inside = |r: &[f64; 4]| r[0] >= rx0 - 1.0 && r[2] <= rx1 + 1.0 && r[1] >= ry0 - 1.0 && r[3] <= ry1 + 1.0;
            let spot = spots
                .iter()
                .copied()
                .find(|&s| inside(&rect(s)) && !taken.iter().any(|t| overlaps(t, &rect(s))))
                .or_else(|| spots.iter().copied().find(|&s| inside(&rect(s))))
                .unwrap_or(spots[0]);
            let x = spot.0.clamp(rx0 + tw / 2.0, (rx1 - tw / 2.0).max(rx0 + tw / 2.0));
            let y = spot.1.clamp(ry0 + gsize * 0.8, (ry1 - gsize * 0.2).max(ry0 + gsize * 0.8));
            taken.push(rect((x, y)));
            let color = if style.background == Background::Black { group.color.clone() } else { mix(&group.color, "#000000", 0.25) };
            let _ = write!(out, r#"<text x="{x:.2}" y="{y:.2}" fill="{color}">{}</text>"#, esc(&group.label));
        }
        out.push_str("</g>");
    }
    // the legend, block by block
    if !f.legend.is_empty() {
        let line = pt(style.legend_text.size) * 1.45;
        let _ = write!(out, r#"<g {} fill="{ink}">"#, style.legend_text.attrs(&font));
        let swatch = line * 1.3;
        for (x0, top, rows) in &f.legend {
            let (x0, mut y) = (*x0, *top);
            for row in rows {
                match row {
                    Row::Heading(t) => {
                        let _ = write!(out, r#"<text x="{x0:.1}" y="{:.1}" font-weight="bold">{}</text>"#, y + line * 0.75, esc(t));
                        y += line;
                    }
                    Row::Dot { color, outline, text } => {
                        let _ = write!(
                            out,
                            r#"<circle cx="{:.1}" cy="{:.1}" r="{:.1}" fill="{color}" stroke="{outline}" stroke-width="{:.2}"/>"#,
                            x0 + swatch / 2.0,
                            y + line * 0.5,
                            line * 0.32,
                            pt(0.6)
                        );
                        let _ = write!(out, r#"<text x="{:.1}" y="{:.1}">{}</text>"#, x0 + line * 1.8, y + line * 0.75, esc(text));
                        y += line;
                    }
                    Row::Line { color, width, text } => {
                        let _ = write!(
                            out,
                            r#"<line x1="{x0:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{color}" stroke-width="{width:.2}" stroke-linecap="round"/>"#,
                            y + line * 0.5,
                            x0 + swatch,
                            y + line * 0.5
                        );
                        let _ = write!(out, r#"<text x="{:.1}" y="{:.1}">{}</text>"#, x0 + line * 1.8, y + line * 0.75, esc(text));
                        y += line;
                    }
                    Row::Loop { color, text } => {
                        let _ = write!(
                            out,
                            r#"<circle cx="{:.1}" cy="{:.1}" r="{:.1}" fill="none" stroke="{color}" stroke-width="{:.2}"/>"#,
                            x0 + swatch / 2.0,
                            y + line * 0.5,
                            line * 0.3,
                            pt(0.9)
                        );
                        let _ = write!(out, r#"<text x="{:.1}" y="{:.1}">{}</text>"#, x0 + line * 1.8, y + line * 0.75, esc(text));
                        y += line;
                    }
                    Row::Size { diameter, text } => {
                        let height = diameter.max(line) + line * 0.2;
                        let _ = write!(
                            out,
                            r#"<circle cx="{:.1}" cy="{:.1}" r="{:.1}" fill="{}" stroke="{}" stroke-width="{:.2}"/>"#,
                            x0 + diameter / 2.0,
                            y + height / 2.0,
                            diameter / 2.0,
                            style.ungrouped_color,
                            style.node_outline,
                            pt(0.6)
                        );
                        let _ = write!(out, r#"<text x="{:.1}" y="{:.1}">{}</text>"#, x0 + diameter.max(line) + line * 0.6, y + height / 2.0 + line * 0.25, esc(text));
                        y += height;
                    }
                    Row::Text(t) => {
                        let _ = write!(out, r#"<text x="{x0:.1}" y="{:.1}">{}</text>"#, y + line * 0.75, esc(t));
                        y += line;
                    }
                    Row::Gap => y += line * 0.5,
                }
            }
        }
        out.push_str("</g>");
    }
    out.push_str("</svg>");
    out
}

// ---------------------------------------------------------------- Cytoscape, GraphML and tables

/// Names joined for a table cell
fn names(terms: &[Term]) -> String {
    terms.iter().map(|t| if t.id.is_empty() { t.name.clone() } else { format!("{} [{}]", t.name, t.id) }).collect::<Vec<_>>().join("; ")
}

/// A node's attributes for the exports, in column order
fn node_columns(data: &NetworkData, ann: &Annotations, plan: &Plan, f: &Frame, k: usize) -> Vec<(&'static str, serde_json::Value)> {
    use serde_json::json;
    let n = &data.nodes[plan.nodes[k]];
    let p = ann.proteins.get(&n.accession);
    let text = |s: Option<&str>| json!(s.unwrap_or(""));
    let group = plan.group_of[k].map(|g| plan.groups[g].label.clone()).unwrap_or_default();
    vec![
        ("id", json!(n.accession)),
        ("label", json!(f.labels[k])),
        ("gene", json!(p.map(|p| p.gene.as_str()).filter(|g| !g.is_empty()).unwrap_or(n.gene.as_str()))),
        ("protein", json!(p.map(|p| p.name.as_str()).filter(|g| !g.is_empty()).unwrap_or(n.name.as_str()))),
        ("uniprot", text(p.map(|p| p.accession.as_str()))),
        ("entry", json!(p.map(|p| p.entry.as_str()).filter(|g| !g.is_empty()).unwrap_or(n.entry.as_str()))),
        ("organism", text(p.map(|p| p.organism.as_str()))),
        ("group", json!(group)),
        ("color", json!(f.colors[k])),
        ("crosslinks", json!(n.crosslinks)),
        ("csms", json!(n.csms)),
        ("partners", json!(n.partners)),
        ("self_links", json!(n.self_links)),
        ("location", json!(p.map(|p| p.locations.join("; ")).unwrap_or_default())),
        ("complexes", json!(p.map(|p| names(&p.complexes)).unwrap_or_default())),
        ("pathways", json!(p.map(|p| names(&p.pathways)).unwrap_or_default())),
        ("go_process", json!(p.map(|p| names(&p.go_process)).unwrap_or_default())),
        ("go_function", json!(p.map(|p| names(&p.go_function)).unwrap_or_default())),
        ("go_component", json!(p.map(|p| names(&p.go_component)).unwrap_or_default())),
        ("keywords", json!(p.map(|p| names(&p.keywords)).unwrap_or_default())),
        ("function", text(p.map(|p| p.function.as_str()))),
        ("string_id", text(p.map(|p| p.string_id.as_str()))),
        ("diameter", json!((f.nodes[k].2 * 2.0 * 10.0).round() / 10.0)),
    ]
}

/// A single-precision number as written (`0.999`, not `0.9990000128746033`)
fn clean(x: f32) -> f64 {
    x.to_string().parse().unwrap_or(x as f64)
}

/// An edge's attributes for the exports, in column order
fn edge_columns(data: &NetworkData, ann: &Annotations, plan: &Plan, f: &Frame, k: usize) -> Vec<(&'static str, serde_json::Value)> {
    use serde_json::json;
    let e = &data.edges[plan.edges[k]];
    let interaction = ann.interaction(&e.a, &e.b).flatten();
    let shared = |by: GroupBy| -> String {
        match (ann.proteins.get(&e.a), ann.proteins.get(&e.b)) {
            (Some(p), Some(q)) => {
                let theirs: HashSet<String> = by.terms(q).iter().map(|t| t.key().to_string()).collect();
                names(&by.terms(p).into_iter().filter(|t| theirs.contains(t.key())).collect::<Vec<_>>())
            }
            _ => String::new(),
        }
    };
    let score = |s: Option<f32>| interaction.map_or(json!(null), |_| json!(clean(s.unwrap_or(0.0))));
    let look = &f.edges[k];
    vec![
        ("source", json!(e.a)),
        ("target", json!(e.b)),
        ("crosslinks", json!(e.crosslinks)),
        ("csms", json!(e.csms)),
        ("best_q", json!(clean(e.q))),
        ("best_score", json!((e.score * 100.0).round() / 100.0)),
        ("multivalent", json!(e.multivalent)),
        ("residue_pairs", json!(e.sites.iter().map(|(a, b)| format!("{a}-{b}")).collect::<Vec<_>>().join("; "))),
        ("known", look.known.map_or(json!(null), |k| json!(k))),
        ("string_score", score(interaction.map(|i| i.score))),
        ("string_experiments", score(interaction.map(|i| i.experiments))),
        ("string_databases", score(interaction.map(|i| i.databases))),
        ("string_textmining", score(interaction.map(|i| i.textmining))),
        ("shared_complexes", json!(shared(GroupBy::Complex))),
        ("shared_pathways", json!(shared(GroupBy::Pathway))),
        ("color", json!(look.color)),
        ("width", json!((look.width * 100.0).round() / 100.0)),
    ]
}

/// The network for Cytoscape (`.cyjs`, Cytoscape.js JSON): every attribute,
/// and the layout as drawn
pub fn cytoscape(data: &NetworkData, ann: &Annotations, plan: &Plan, style: &NetworkStyle) -> String {
    use serde_json::{json, Map, Value};
    let f = frame(data, ann, plan, style);
    let nodes: Vec<Value> = (0..plan.nodes.len())
        .map(|k| {
            let mut fields: Map<String, Value> = node_columns(data, ann, plan, &f, k).into_iter().map(|(k, v)| (k.to_string(), v)).collect();
            fields.insert("name".into(), fields["label"].clone());
            json!({"data": fields, "position": {"x": (f.nodes[k].0 * 10.0).round() / 10.0, "y": (f.nodes[k].1 * 10.0).round() / 10.0}})
        })
        .collect();
    let edges: Vec<Value> = (0..plan.edges.len())
        .map(|k| {
            let mut fields: Map<String, Value> = edge_columns(data, ann, plan, &f, k).into_iter().map(|(k, v)| (k.to_string(), v)).collect();
            fields.insert("id".into(), json!(format!("e{}", k + 1)));
            fields.insert("interaction".into(), json!("crosslink"));
            json!({"data": fields})
        })
        .collect();
    let name = if style.title.trim().is_empty() { "XNova crosslink network".to_string() } else { style.title.clone() };
    serde_json::to_string_pretty(&json!({
        "format_version": "1.0",
        "generated_by": concat!("XNova ", env!("CARGO_PKG_VERSION")),
        "target_cytoscapejs_version": "~2.1",
        "data": {"name": name, "uniprot_release": ann.uniprot_release, "string_version": ann.string_version, "looked_up": ann.looked_up},
        "elements": {"nodes": nodes, "edges": edges}
    }))
    .unwrap_or_default()
}

/// The network as GraphML (Cytoscape, Gephi, yEd, igraph, NetworkX)
pub fn graphml(data: &NetworkData, ann: &Annotations, plan: &Plan, style: &NetworkStyle) -> String {
    let f = frame(data, ann, plan, style);
    let kind = |v: &serde_json::Value| match v {
        serde_json::Value::Bool(_) => "boolean",
        serde_json::Value::Number(n) if n.is_u64() || n.is_i64() => "int",
        serde_json::Value::Number(_) => "double",
        _ => "string",
    };
    let plain = |v: &serde_json::Value| match v {
        serde_json::Value::String(s) => esc(s),
        serde_json::Value::Null => String::new(),
        v => v.to_string(),
    };
    let mut out = String::from("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\">\n");
    let node_rows: Vec<_> = (0..plan.nodes.len()).map(|k| node_columns(data, ann, plan, &f, k)).collect();
    let edge_rows: Vec<_> = (0..plan.edges.len()).map(|k| edge_columns(data, ann, plan, &f, k)).collect();
    if let Some(first) = node_rows.first() {
        for (name, v) in first.iter().skip(1) {
            let _ = writeln!(out, r#"  <key id="n_{name}" for="node" attr.name="{name}" attr.type="{}"/>"#, kind(v));
        }
        let _ = writeln!(out, r#"  <key id="n_x" for="node" attr.name="x" attr.type="double"/>"#);
        let _ = writeln!(out, r#"  <key id="n_y" for="node" attr.name="y" attr.type="double"/>"#);
    }
    if let Some(first) = edge_rows.first() {
        for (name, v) in first.iter().skip(2) {
            // known and scores may be empty for some edges: typed by their kind when present
            let t = match *name {
                "known" => "boolean",
                n if n.starts_with("string_") => "double",
                _ => kind(v),
            };
            let _ = writeln!(out, r#"  <key id="e_{name}" for="edge" attr.name="{name}" attr.type="{t}"/>"#);
        }
    }
    let _ = writeln!(out, r#"  <graph id="crosslinks" edgedefault="undirected">"#);
    for (k, row) in node_rows.iter().enumerate() {
        let _ = writeln!(out, r#"    <node id="{}">"#, esc(row[0].1.as_str().unwrap_or("")));
        for (name, v) in row.iter().skip(1) {
            let _ = writeln!(out, r#"      <data key="n_{name}">{}</data>"#, plain(v));
        }
        let _ = writeln!(out, r#"      <data key="n_x">{:.1}</data>"#, f.nodes[k].0);
        let _ = writeln!(out, r#"      <data key="n_y">{:.1}</data>"#, f.nodes[k].1);
        out.push_str("    </node>\n");
    }
    for (k, row) in edge_rows.iter().enumerate() {
        let _ = writeln!(
            out,
            r#"    <edge id="e{}" source="{}" target="{}">"#,
            k + 1,
            esc(row[0].1.as_str().unwrap_or("")),
            esc(row[1].1.as_str().unwrap_or(""))
        );
        for (name, v) in row.iter().skip(2) {
            if !v.is_null() {
                let _ = writeln!(out, r#"      <data key="e_{name}">{}</data>"#, plain(v));
            }
        }
        out.push_str("    </edge>\n");
    }
    out.push_str("  </graph>\n</graphml>\n");
    out
}

fn tsv(rows: Vec<Vec<(&'static str, serde_json::Value)>>) -> String {
    let mut out = String::new();
    let cell = |v: &serde_json::Value| match v {
        serde_json::Value::String(s) => s.replace(['\t', '\n', '\r'], " "),
        serde_json::Value::Null => String::new(),
        v => v.to_string(),
    };
    if let Some(first) = rows.first() {
        out.push_str(&first.iter().map(|(k, _)| *k).collect::<Vec<_>>().join("\t"));
        out.push('\n');
    }
    for row in &rows {
        out.push_str(&row.iter().map(|(_, v)| cell(v)).collect::<Vec<_>>().join("\t"));
        out.push('\n');
    }
    out
}

/// The proteins and what is known about them, a row each
pub fn node_table(data: &NetworkData, ann: &Annotations, plan: &Plan, style: &NetworkStyle) -> String {
    let f = frame(data, ann, plan, style);
    tsv((0..plan.nodes.len()).map(|k| node_columns(data, ann, plan, &f, k)).collect())
}

/// The crosslinked pairs of proteins, a row each
pub fn edge_table(data: &NetworkData, ann: &Annotations, plan: &Plan, style: &NetworkStyle) -> String {
    let f = frame(data, ann, plan, style);
    tsv((0..plan.edges.len()).map(|k| edge_columns(data, ann, plan, &f, k)).collect())
}

/// Write a network into `dir`: `crosslink_network.svg`, `.pdf` and `.png`
/// (300 dpi), the style (`network_style.json`), `crosslink_network.cyjs`
/// (Cytoscape), `.graphml` and the node and edge tables. Returns the files.
pub fn write_all(data: &NetworkData, ann: &Annotations, plan: &Plan, style: &NetworkStyle, dir: &Path) -> Result<Vec<PathBuf>, String> {
    std::fs::create_dir_all(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
    let svg = svg(data, ann, plan, style);
    let mut written = Vec::new();
    for format in [figure::Format::Svg, figure::Format::Pdf, figure::Format::Png] {
        let path = dir.join(format!("crosslink_network.{}", format.extension()));
        figure::export(&svg, format, 300, style.background != Background::Transparent, &path)?;
        written.push(path);
    }
    for (name, text) in [
        ("network_style.json", style.to_json()),
        ("crosslink_network.cyjs", cytoscape(data, ann, plan, style)),
        ("crosslink_network.graphml", graphml(data, ann, plan, style)),
        ("network_proteins.tsv", node_table(data, ann, plan, style)),
        ("network_pairs.tsv", edge_table(data, ann, plan, style)),
    ] {
        let path = dir.join(name);
        std::fs::write(&path, text).map_err(|e| format!("{}: {e}", path.display()))?;
        written.push(path);
    }
    Ok(written)
}

/// Write a search's network with `style` into the results folder's
/// [`FOLDER`], with what was looked up about its proteins (kept there; and
/// looked up first if `look_up`). A search without crosslinks gets none.
/// Returns the files written, and the lookups' problems.
pub fn write_figures(folder: &Path, style: &NetworkStyle, look_up: bool) -> Result<(Vec<PathBuf>, Vec<String>), String> {
    let results = Results::load(folder)?;
    if !results.matches.iter().any(|m| m.kind != Kind::Linear) {
        return Ok((Vec::new(), Vec::new()));
    }
    let data = build(&results, style);
    if data.nodes.is_empty() {
        return Ok((Vec::new(), Vec::new()));
    }
    let mut ann = Annotations::load(folder).unwrap_or_default();
    let mut problems = Vec::new();
    let mut written = Vec::new();
    if look_up {
        let (found, issues) = super::protein_info::look_up(&data.accessions(), &ann, false, &|step: String| eprintln!("{step}"));
        ann = found;
        problems = issues;
        written.push(ann.save(folder)?);
    }
    let own = OwnGroups::read(&style.own_groups);
    problems.extend(own.notes.iter().cloned());
    let plan = plan(&data, &ann, &own, style);
    written.extend(write_all(&data, &ann, &plan, style, &folder.join(FOLDER))?);
    Ok((written, problems))
}

#[cfg(test)]
mod test {
    use super::super::protein_info::{Interaction, ProteinInfo};
    use super::*;

    fn link(a: (&str, usize), b: (&str, usize), class: Class, csms: u32) -> circos::Link {
        circos::Link { a: (a.0.into(), a.1), b: (b.0.into(), b.1), class, csms, q: 0.001, score: 20.0 }
    }

    /// Two complexes (A, B, C and D, E) joined by one crosslink, with loops
    /// on A and a protein crosslinked only within itself (F)
    fn data() -> NetworkData {
        let links = vec![
            link(("A", 10), ("B", 20), Class::InterProtein, 3),
            link(("A", 12), ("B", 25), Class::InterProtein, 1),
            link(("B", 5), ("A", 40), Class::InterProtein, 2),
            link(("A", 1), ("C", 2), Class::InterProtein, 1),
            link(("B", 7), ("C", 9), Class::Multivalent, 1),
            link(("C", 3), ("D", 4), Class::InterProtein, 1),
            link(("D", 5), ("E", 6), Class::InterProtein, 4),
            link(("A", 3), ("A", 90), Class::IntraProtein, 2),
            link(("F", 3), ("F", 9), Class::Loop, 1),
        ];
        from_links(&links, |a| (format!("G{a}"), format!("protein {a}")))
    }

    fn info(gene: &str, complexes: &[(&str, &str)], location: &str) -> ProteinInfo {
        ProteinInfo {
            gene: gene.into(),
            complexes: complexes.iter().map(|(id, name)| Term { id: id.to_string(), name: name.to_string() }).collect(),
            locations: vec![location.into()],
            ..Default::default()
        }
    }

    fn annotations() -> Annotations {
        let mut ann = Annotations::default();
        let pol = ("CPX-1", "Polymerase");
        let spl = ("CPX-2", "Spliceosome");
        for (a, gene, complexes) in [("A", "RPA", vec![pol]), ("B", "RPB", vec![pol]), ("C", "RPC", vec![pol, spl]), ("D", "SPD", vec![spl]), ("E", "SPE", vec![spl])] {
            ann.proteins.insert(a.into(), info(gene, &complexes, "Nucleus"));
        }
        ann.string_queried.extend(["A", "B", "C", "D", "E"].map(String::from));
        ann.interactions.push(Interaction { a: "A".into(), b: "B".into(), score: 0.99, experiments: 0.9, databases: 0.9, ..Default::default() });
        ann.interactions.push(Interaction { a: "D".into(), b: "E".into(), score: 0.5, textmining: 0.5, ..Default::default() });
        ann
    }

    #[test]
    fn nodes_and_edges_from_crosslinks() {
        let d = data();
        assert_eq!(d.nodes.len(), 6);
        assert_eq!(d.edges.len(), 5, "A-B, A-C, B-C, C-D, D-E: {:?}", d.edges.iter().map(|e| (&e.a, &e.b)).collect::<Vec<_>>());
        let ab = d.edges.iter().find(|e| e.a == "A" && e.b == "B").unwrap();
        assert_eq!((ab.crosslinks, ab.csms), (3, 6));
        // the residue in A first, whichever way round the link was
        assert_eq!(ab.sites, [(10, 20), (12, 25), (40, 5)]);
        let a = d.node("A").unwrap();
        assert_eq!((a.crosslinks, a.partners, a.self_links, a.self_csms), (4, 2, 1, 2));
        assert_eq!(d.edges.iter().find(|e| e.a == "B" && e.b == "C").unwrap().multivalent, 1);
        let f = d.node("F").unwrap();
        assert_eq!((f.partners, f.self_links), (0, 1));
        assert_eq!(entry_of("sp|P04050|RPB1_YEAST"), "RPB1_YEAST");
        assert_eq!(entry_of("P04050"), "P04050");
    }

    #[test]
    fn groups_cover_the_most_proteins_first() {
        let terms = |names: &[&str]| names.iter().map(|n| Term { id: String::new(), name: n.to_string() }).collect::<Vec<_>>();
        let items = vec![terms(&["nucleus", "pol"]), terms(&["nucleus", "pol"]), terms(&["nucleus", "pol", "spl"]), terms(&["nucleus", "spl"]), terms(&["spl"]), terms(&[])];
        let (taken, of, rest) = cover(&items, 2, 10, &[]);
        // the general term takes most; the one left with spl is no group of 2
        assert_eq!(taken, [("nucleus".into(), "nucleus".into(), 4)]);
        assert_eq!(of, [Some(0), Some(0), Some(0), Some(0), None, None]);
        assert_eq!(rest, [("pol".into(), "pol".into(), 3), ("spl".into(), "spl".into(), 3)]);
        // without the general term: two groups
        let items: Vec<Vec<Term>> = items.into_iter().map(|t| t.into_iter().filter(|t| t.name != "nucleus").collect()).collect();
        let (taken, of, rest) = cover(&items, 2, 10, &[]);
        assert_eq!(taken.iter().map(|t| t.0.as_str()).collect::<Vec<_>>(), ["pol", "spl"]);
        assert_eq!(of, [Some(0), Some(0), Some(0), Some(1), Some(1), None]);
        assert!(rest.is_empty());
        // at most one group
        let (taken, ..) = cover(&items, 2, 1, &[]);
        assert_eq!(taken.len(), 1);
        // a preferred term first, though another has more
        let (taken, of, _) = cover(&items, 2, 10, &["spl".to_string(), "missing".to_string()]);
        assert_eq!(taken.iter().map(|t| (t.0.as_str(), t.2)).collect::<Vec<_>>(), [("spl", 3), ("pol", 2)]);
        assert_eq!(of, [Some(1), Some(1), Some(0), Some(0), Some(0), None]);
    }

    #[test]
    fn modules() {
        // two triangles joined by one edge
        let edges = [(0, 1, 1.0), (1, 2, 1.0), (0, 2, 1.0), (3, 4, 1.0), (4, 5, 1.0), (3, 5, 1.0), (2, 3, 1.0)];
        let (community, q) = louvain(6, &edges);
        assert_eq!(community[0], community[1]);
        assert_eq!(community[1], community[2]);
        assert_eq!(community[3], community[5]);
        assert_ne!(community[0], community[3]);
        assert!((q - 5.0 / 14.0).abs() < 1e-9, "{q}");
        assert_eq!(modularity(&[0; 6], &edges), 0.0);
        // nothing to divide
        assert_eq!(louvain(3, &[]).0, [0, 1, 2]);
    }

    #[test]
    fn plans_and_layouts() {
        let d = data();
        let ann = annotations();
        let style = NetworkStyle::default();
        let p = plan(&d, &ann, &OwnGroups::default(), &style);
        assert_eq!(p.nodes.len(), 6);
        assert_eq!(p.groups.iter().map(|g| (g.label.as_str(), g.members)).collect::<Vec<_>>(), [("Polymerase", 3), ("Spliceosome", 2)]);
        let group = |a: &str| p.group_of[p.nodes.iter().position(|&i| d.nodes[i].accession == a).unwrap()];
        assert_eq!((group("A"), group("C"), group("D"), group("F")), (Some(0), Some(0), Some(1), None));
        // laid out within the area, the same every time
        assert!(p.layout.iter().all(|q| (-1e-9..=1.0 + 1e-9).contains(&q[0]) && (-1e-9..=1.0 + 1e-9).contains(&q[1])));
        assert_eq!(plan(&d, &ann, &OwnGroups::default(), &style).layout, p.layout);
        assert_ne!(plan(&d, &ann, &OwnGroups::default(), &NetworkStyle { seed: 7, ..style.clone() }).layout, p.layout);
        // no lone proteins, one hidden, a term excluded
        let fewer = NetworkStyle { lone: false, hidden: vec!["E".into()], excluded: vec!["CPX-2".into()], ..style.clone() };
        let p2 = plan(&d, &ann, &OwnGroups::default(), &fewer);
        assert_eq!(p2.nodes.len(), 4);
        assert_eq!(p2.edges.len(), 4);
        assert_eq!(p2.groups.len(), 1);
        // modules
        let by_module = plan(&d, &ann, &OwnGroups::default(), &NetworkStyle { group_by: GroupBy::Module, ..style.clone() });
        assert!(by_module.modularity.is_some());
        assert_eq!(by_module.groups[0].label, "Module 1");
        // circles and clusters
        for layout in [Layout::Circle, Layout::Groups] {
            let p3 = plan(&d, &ann, &OwnGroups::default(), &NetworkStyle { layout, ..style.clone() });
            assert_eq!(p3.layout.len(), 6);
        }
    }

    #[test]
    fn own_groups_from_a_file() {
        let path = std::env::temp_dir().join(format!("xnova_groups_{}.tsv", std::process::id()));
        std::fs::write(&path, "protein\tgroup\tcolour\nGA\tcore\t#ff0000\nB\tcore\n# a comment\nGD\tsplicing\nnothing\n").unwrap();
        let own = OwnGroups::read(path.to_str().unwrap());
        let _ = std::fs::remove_file(&path);
        assert_eq!(own.of.len(), 3);
        assert_eq!(own.colors.get("core").map(String::as_str), Some("#ff0000"));
        assert_eq!(own.notes.len(), 1, "{:?}", own.notes);
        let d = data();
        let style = NetworkStyle { group_by: GroupBy::Own, min_group: 1, ..NetworkStyle::default() };
        let p = plan(&d, &Annotations::default(), &own, &style);
        assert_eq!(p.groups.iter().map(|g| (g.label.as_str(), g.members, g.color.as_str())).collect::<Vec<_>>()[0], ("core", 2, "#ff0000"));
        assert_eq!(p.groups.len(), 2);
    }

    #[test]
    fn edge_colours() {
        let d = data();
        let ann = annotations();
        let style = NetworkStyle::default();
        let p = plan(&d, &ann, &OwnGroups::default(), &style);
        let f = frame(&d, &ann, &p, &style);
        let look = |a: &str, b: &str| &f.edges[p.edges.iter().position(|&e| d.edges[e].a == a && d.edges[e].b == b).unwrap()];
        // known by experiments or databases; D-E by text mining only: new
        assert_eq!(look("A", "B").known, Some(true));
        assert_eq!(look("A", "B").color, style.known_color);
        assert_eq!(look("D", "E").known, Some(false));
        assert_eq!(look("A", "C").known, Some(false));
        // the combined score counts text mining
        let combined = NetworkStyle { evidence: Evidence::Combined, ..style.clone() };
        assert_eq!(frame(&d, &ann, &p, &combined).edges[p.edges.iter().position(|&e| d.edges[e].a == "D").unwrap()].known, Some(true));
        // nothing looked up
        let bare = frame(&d, &Annotations::default(), &p, &style);
        assert!(bare.edges.iter().all(|e| e.known.is_none()));
        // shared complexes: the groups' colours (grouped by complex too)
        let by_complex = NetworkStyle { color_by: EdgeColor::Complex, ..style.clone() };
        let f = frame(&d, &ann, &p, &by_complex);
        let row = |a: &str, b: &str| f.edge_rows[f.edges[p.edges.iter().position(|&e| d.edges[e].a == a && d.edges[e].b == b).unwrap()].row].1.clone();
        assert_eq!(row("A", "B"), "within Polymerase (3)");
        assert_eq!(row("C", "D"), "within Spliceosome (2)");
        assert!(f.edge_rows.iter().all(|(_, t)| !t.contains("no complex")), "{:?}", f.edge_rows);
        // wider the more crosslinks
        assert!(look("A", "B").width > look("A", "C").width);
    }

    #[test]
    fn figure_and_exports() {
        let d = data();
        let ann = annotations();
        let style = NetworkStyle { title: "Two complexes".into(), numbers: true, ..NetworkStyle::default() };
        let p = plan(&d, &ann, &OwnGroups::default(), &style);
        let svg = svg(&d, &ann, &p, &style);
        assert!(svg.contains(">Two complexes<"));
        assert!(svg.contains(">RPA<"), "UniProt's gene names as labels");
        assert!(svg.contains(">GF<"), "the FASTA's where UniProt has none");
        assert!(svg.contains(">Polymerase (3)<") && svg.contains(">Spliceosome (2)<") && svg.contains(">other (1)<"));
        assert!(svg.contains(">known to STRING (1)<"));
        assert!(svg.contains("known: STRING score ≥ 0.40, experiments or databases"));
        // a long name, shortened to the legend's room
        let long = fit(&"x".repeat(500), 100.0, &style.legend_text);
        assert!(long.ends_with('…') && text_width(&long, &style.legend_text) <= 100.0 + text_width("x", &style.legend_text));
        assert!(svg.contains("6 proteins, 5 crosslinked pairs"));
        assert!(figure::rasterize(&svg, 30).is_ok());
        // a node moved by hand
        let mut moved = style.clone();
        moved.positions.insert("A".into(), [0.0, 0.0]);
        let f = frame(&d, &ann, &p, &moved);
        let a = p.nodes.iter().position(|&i| d.nodes[i].accession == "A").unwrap();
        assert_eq!((f.nodes[a].0, f.nodes[a].1), (f.layout_box[0], f.layout_box[1]));
        assert!(f.layout_box[0] >= f.area[0] - 1e-9 && f.layout_box[2] <= f.area[2] + 1e-9);
        // Cytoscape: the attributes and the places
        let cy: serde_json::Value = serde_json::from_str(&cytoscape(&d, &ann, &p, &style)).unwrap();
        let nodes = cy["elements"]["nodes"].as_array().unwrap();
        assert_eq!(nodes.len(), 6);
        let node_a = nodes.iter().find(|n| n["data"]["id"] == "A").unwrap();
        assert_eq!(node_a["data"]["complexes"], "Polymerase [CPX-1]");
        assert_eq!(node_a["data"]["group"], "Polymerase");
        assert!(node_a["position"]["x"].is_number());
        let edges = cy["elements"]["edges"].as_array().unwrap();
        let ab = edges.iter().find(|e| e["data"]["source"] == "A" && e["data"]["target"] == "B").unwrap();
        assert_eq!(ab["data"]["known"], true);
        assert_eq!(ab["data"]["string_score"], 0.99, "as STRING gave it");
        assert_eq!(ab["data"]["best_q"], 0.001);
        assert_eq!(ab["data"]["residue_pairs"], "10-20; 12-25; 40-5");
        assert_eq!(ab["data"]["shared_complexes"], "Polymerase [CPX-1]");
        // GraphML that reads back
        let xml = graphml(&d, &ann, &p, &style);
        let doc = roxmltree::Document::parse(&xml).unwrap();
        assert_eq!(doc.descendants().filter(|n| n.has_tag_name("node")).count(), 6);
        assert_eq!(doc.descendants().filter(|n| n.has_tag_name("edge")).count(), 5);
        // tables
        let nodes = node_table(&d, &ann, &p, &style);
        assert_eq!(nodes.lines().count(), 7);
        assert!(nodes.starts_with("id\tlabel\tgene\tprotein\tuniprot"));
        assert_eq!(edge_table(&d, &ann, &p, &style).lines().count(), 6);
        // styles as JSON, missing keys their default
        let back = NetworkStyle::from_json(&moved.to_json()).unwrap();
        assert_eq!(back, moved);
        assert_eq!(NetworkStyle::from_json(r#"{"width": 5}"#).unwrap().height, 6.0);
        // no proteins at all
        let empty = svg_of_nothing();
        assert!(empty.contains("No crosslinked proteins"));
    }

    fn svg_of_nothing() -> String {
        let d = NetworkData::default();
        let p = plan(&d, &Annotations::default(), &OwnGroups::default(), &NetworkStyle::default());
        svg(&d, &Annotations::default(), &p, &NetworkStyle::default())
    }

    /// A larger network lays out in reasonable time, its nodes spread out
    #[test]
    fn large_layout() {
        let n = 400;
        let mut random = Random::new(3);
        let mut links = Vec::new();
        for i in 0..n {
            for _ in 0..2 {
                let j = (random.next() * n as f64) as usize % n;
                if i != j {
                    links.push(link((&format!("P{i}"), 1), (&format!("P{j}"), 2), Class::InterProtein, 1));
                }
            }
        }
        let d = from_links(&links, |_| (String::new(), String::new()));
        let started = std::time::Instant::now();
        let p = plan(&d, &Annotations::default(), &OwnGroups::default(), &NetworkStyle { group_by: GroupBy::Module, ..NetworkStyle::default() });
        let took = started.elapsed();
        assert!(took.as_secs_f64() < 20.0, "{took:?}");
        assert!(p.groups.len() > 2, "modules found: {}", p.groups.len());
        // spread: few nodes on top of each other
        let close = (0..p.layout.len())
            .flat_map(|i| (i + 1..p.layout.len()).map(move |j| (i, j)))
            .filter(|&(i, j)| (p.layout[i][0] - p.layout[j][0]).hypot(p.layout[i][1] - p.layout[j][1]) < 0.003)
            .count();
        assert!(close < n / 20, "{close} pairs on top of each other");
    }
}
