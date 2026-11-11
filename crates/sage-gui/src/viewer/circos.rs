//! Circos-style circular plots of a search's crosslinks, drawn by XNova
//! itself (no Circos, Perl or internet needed): each protein an arc of the
//! circle, scaled to its length, with residue ticks; each crosslink a curve
//! between its two linked residues; tracks inside the ring for the
//! crosslinked sites, the linker's reactive residues, the linear search's
//! sequence coverage and the user's own annotations (domains, subunits).
//! Everything about the figure is in [`CircosStyle`], which is saved and
//! loaded as JSON, and kept under a name as a [`Template`] in sage-gui's
//! library, to give other results the same look.
//!
//! The same plot can be written as Circos input files (karyotype, links,
//! tracks and a `circos.conf`), for Circos itself (circos.ca, Martin
//! Krzywinski), and as a protein-by-protein table for Circos' table viewer.

use super::figure;
use super::mobiligram::Class;
use super::model::{parse_species, Kind, Match, Results};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::fmt::Write;
use std::path::{Path, PathBuf};

/// Folder of a results folder the plots of a search go to
pub const FOLDER: &str = "circos";
/// SVG user units per inch (as the other figures)
const UNITS_PER_IN: f64 = 150.0;

// ---------------------------------------------------------------- the style

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Palette {
    /// Tableau's 20 colours (pairs of a strong and a light tone)
    Tableau20,
    /// Okabe & Ito's 8 colours, distinguishable with colour-vision deficiency
    OkabeIto,
    /// Every protein grey
    Grey,
}

impl Palette {
    pub const ALL: [Palette; 3] = [Palette::Tableau20, Palette::OkabeIto, Palette::Grey];

    pub fn label(self) -> &'static str {
        match self {
            Palette::Tableau20 => "Tableau 20",
            Palette::OkabeIto => "Okabe-Ito (colour-blind safe)",
            Palette::Grey => "grey",
        }
    }

    fn colors(self) -> &'static [&'static str] {
        match self {
            Palette::Tableau20 => &[
                "#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f", "#edc948", "#b07aa1", "#ff9da7", "#9c755f",
                "#bab0ac", "#a0cbe8", "#ffbe7d", "#ff9d9a", "#86bcb6", "#8cd17d", "#f1ce63", "#d4a6c8", "#fabfd2",
                "#d7b5a6", "#79706e",
            ],
            Palette::OkabeIto => &["#e69f00", "#56b4e9", "#009e73", "#f0e442", "#0072b2", "#d55e00", "#cc79a7", "#999999"],
            Palette::Grey => &["#b8b8b8"],
        }
    }

    pub fn color(self, i: usize) -> &'static str {
        let colors = self.colors();
        colors[i % colors.len()]
    }
}

/// What gives a link its colour
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ColorBy {
    /// Between proteins, within a protein, loop-link, multivalent
    Class,
    /// The colour of the protein of its first end
    Protein,
    /// Its best hyperscore, on a gradient
    Score,
    /// Its crosslink q-value, on a gradient (lowest q: high colour)
    QValue,
    /// Its number of crosslink-spectrum matches, on a gradient
    Csms,
    /// One colour for all
    Single,
}

impl ColorBy {
    pub const ALL: [ColorBy; 6] = [ColorBy::Class, ColorBy::Protein, ColorBy::Score, ColorBy::QValue, ColorBy::Csms, ColorBy::Single];

    pub fn label(self) -> &'static str {
        match self {
            ColorBy::Class => "kind of crosslink",
            ColorBy::Protein => "protein",
            ColorBy::Score => "score (hyperscore)",
            ColorBy::QValue => "crosslink q-value",
            ColorBy::Csms => "number of CSMs",
            ColorBy::Single => "one colour",
        }
    }

}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Background {
    White,
    Transparent,
    Black,
}

impl Background {
    pub const ALL: [Background; 3] = [Background::White, Background::Transparent, Background::Black];

    pub fn label(self) -> &'static str {
        match self {
            Background::White => "white",
            Background::Transparent => "transparent",
            Background::Black => "black",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum LabelMode {
    /// Along the arc
    Parallel,
    /// Pointing out from the centre
    Radial,
    None,
}

impl LabelMode {
    pub const ALL: [LabelMode; 3] = [LabelMode::Parallel, LabelMode::Radial, LabelMode::None];

    pub fn label(self) -> &'static str {
        match self {
            LabelMode::Parallel => "along the arc",
            LabelMode::Radial => "radial",
            LabelMode::None => "none",
        }
    }
}

/// A kind of crosslink, for a highlight
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum LinkKind {
    Between,
    Within,
    Loop,
    #[default]
    Multivalent,
}

impl LinkKind {
    pub const ALL: [LinkKind; 4] = [LinkKind::Between, LinkKind::Within, LinkKind::Loop, LinkKind::Multivalent];

    pub fn label(self) -> &'static str {
        match self {
            LinkKind::Between => "between proteins",
            LinkKind::Within => "within a protein",
            LinkKind::Loop => "loop-links",
            LinkKind::Multivalent => "multivalent",
        }
    }

    fn class(self) -> Class {
        match self {
            LinkKind::Between => Class::InterProtein,
            LinkKind::Within => Class::IntraProtein,
            LinkKind::Loop => Class::Loop,
            LinkKind::Multivalent => Class::Multivalent,
        }
    }
}

/// The crosslinks a highlight picks; proteins by accession, residues
/// 1-based
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "pick", rename_all = "snake_case")]
pub enum Pick {
    /// One crosslink, by its two residues (either way round)
    Link { a: (String, usize), b: (String, usize) },
    /// Every crosslink with an end in the protein
    Protein { protein: String },
    /// The crosslinks between two proteins (the same one twice: within it)
    Between { a: String, b: String },
    /// Every crosslink of a kind
    Kind { kind: LinkKind },
    /// The crosslinks with an end in a stretch of a protein (inclusive)
    Region { protein: String, start: usize, end: usize },
}

impl Pick {
    pub fn picks(&self, l: &Link) -> bool {
        match self {
            Pick::Link { a, b } => (&l.a == a && &l.b == b) || (&l.a == b && &l.b == a),
            Pick::Protein { protein } => &l.a.0 == protein || &l.b.0 == protein,
            Pick::Between { a, b } => (&l.a.0 == a && &l.b.0 == b) || (&l.a.0 == b && &l.b.0 == a),
            Pick::Kind { kind } => l.class == kind.class(),
            Pick::Region { protein, start, end } => {
                [&l.a, &l.b].iter().any(|(accession, residue)| accession == protein && (*start..=*end).contains(residue))
            }
        }
    }

    /// Picks by protein, so it belongs to a plot's proteins (not a kind of
    /// crosslink, which suits any plot)
    pub fn of_proteins(&self) -> bool {
        !matches!(self, Pick::Kind { .. })
    }

    /// What it picks, in words, the proteins by their labels in the plot
    pub fn describe(&self, data: &CircosData, style: &CircosStyle) -> String {
        let name = |accession: &str| protein_label(data, style, accession);
        let residue = |(accession, at): &(String, usize)| {
            let letter = data
                .proteins
                .iter()
                .find(|p| &p.accession == accession)
                .and_then(|p| p.sequence.chars().nth(at.checked_sub(1)?))
                .map(String::from)
                .unwrap_or_default();
            format!("{} {letter}{at}", name(accession))
        };
        match self {
            Pick::Link { a, b } => format!("{} – {}", residue(a), residue(b)),
            Pick::Protein { protein } => format!("crosslinks of {}", name(protein)),
            Pick::Between { a, b } if a == b => format!("within {}", name(a)),
            Pick::Between { a, b } => format!("{} – {}", name(a), name(b)),
            Pick::Kind { kind } => kind.label().to_string(),
            Pick::Region { protein, start, end } => format!("{} {start}–{end}", name(protein)),
        }
    }
}

/// Crosslinks picked out in a colour of their own
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Highlight {
    #[serde(flatten)]
    pub pick: Pick,
    pub color: String,
    /// The legend's text; empty: what it picks
    #[serde(default)]
    pub label: String,
    /// Drawn (unticked: kept, but not drawn)
    #[serde(default = "drawn_by_default")]
    pub on: bool,
}

fn drawn_by_default() -> bool {
    true
}

impl Highlight {
    pub fn new(pick: Pick, color: &str) -> Self {
        Highlight { pick, color: color.to_string(), label: String::new(), on: true }
    }

    /// Its text in the legend
    pub fn text(&self, data: &CircosData, style: &CircosStyle) -> String {
        if self.label.trim().is_empty() {
            self.pick.describe(data, style)
        } else {
            self.label.trim().to_string()
        }
    }
}

/// Colours for new highlights, in turn: strong, and apart from the kinds'
pub const HIGHLIGHT_COLORS: [&str; 8] = ["#0072b2", "#009e73", "#d55e00", "#cc79a7", "#e69f00", "#56b4e9", "#7b3294", "#000000"];

/// The highlight that colours a crosslink: the first drawn one that picks it
pub fn highlight_of(style: &CircosStyle, l: &Link) -> Option<usize> {
    style.highlights.iter().position(|h| h.on && h.pick.picks(l))
}

/// A protein's label in the plot: the style's, else its gene name
fn protein_label(data: &CircosData, style: &CircosStyle, accession: &str) -> String {
    style
        .labels
        .get(accession)
        .filter(|l| !l.trim().is_empty())
        .cloned()
        .or_else(|| data.proteins.iter().find(|p| p.accession == accession).map(|p| p.label.clone()))
        .unwrap_or_else(|| accession.split('|').nth(1).unwrap_or(accession).to_string())
}

/// The corner of the legend
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum LegendCorner {
    BottomLeft,
    BottomRight,
    TopLeft,
    TopRight,
}

impl LegendCorner {
    pub const ALL: [LegendCorner; 4] = [LegendCorner::BottomLeft, LegendCorner::BottomRight, LegendCorner::TopLeft, LegendCorner::TopRight];

    pub fn label(self) -> &'static str {
        match self {
            LegendCorner::BottomLeft => "bottom left",
            LegendCorner::BottomRight => "bottom right",
            LegendCorner::TopLeft => "top left",
            LegendCorner::TopRight => "top right",
        }
    }

    fn right(self) -> bool {
        matches!(self, LegendCorner::BottomRight | LegendCorner::TopRight)
    }

    fn top(self) -> bool {
        matches!(self, LegendCorner::TopLeft | LegendCorner::TopRight)
    }
}

/// Size and style of one kind of text in the plot
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct TypeStyle {
    /// Points
    pub size: f64,
    pub bold: bool,
    pub italic: bool,
}

impl Default for TypeStyle {
    fn default() -> Self {
        TypeStyle::new(11.0, false)
    }
}

impl TypeStyle {
    pub const fn new(size: f64, bold: bool) -> Self {
        TypeStyle { size, bold, italic: false }
    }

    /// Its size in SVG units
    fn units(&self) -> f64 {
        pt(self.size.clamp(2.0, 96.0))
    }

    /// Its SVG attributes, in `font` (escaped)
    fn attrs(&self, font: &str) -> String {
        format!(
            r#"font-family="{font}" font-size="{:.1}"{}{}"#,
            self.units(),
            if self.bold { r#" font-weight="bold""# } else { "" },
            if self.italic { r#" font-style="italic""# } else { "" }
        )
    }

    /// Circos' font of this style (etc/fonts.conf)
    fn circos_font(&self) -> &'static str {
        match (self.bold, self.italic) {
            (true, true) => "bolditalic",
            (true, false) => "bold",
            (false, true) => "italic",
            (false, false) => "default",
        }
    }
}

/// The generic font families, for whatever font is installed
pub const GENERIC_FONTS: [&str; 3] = ["sans-serif", "serif", "monospace"];

/// A font list for the SVG: `name`, then the generic family to fall back on
/// where it is not installed (the SVG opened on another computer)
pub fn font_list(name: &str) -> String {
    let name = name.trim();
    if GENERIC_FONTS.contains(&name) {
        return name.to_string();
    }
    let lower = name.to_lowercase();
    let has = |keys: &[&str]| keys.iter().any(|k| lower.contains(k));
    let fallback = if has(&["mono", "courier", "consol", "code"]) {
        "monospace"
    } else if !lower.contains("sans")
        && has(&["serif", "times", "georgia", "garamond", "cambria", "palatino", "book", "roman", "minion", "baskerville", "bodoni", "century"])
    {
        "serif"
    } else {
        "sans-serif"
    };
    // plain names as they are, others quoted (CSS)
    let quoted = if name.chars().all(|c| c.is_ascii_alphabetic() || c == ' ' || c == '-') {
        name.to_string()
    } else {
        format!("'{}'", name.replace('\'', ""))
    };
    format!("{quoted}, {fallback}")
}

/// The first font of a font list (the one used where installed)
pub fn first_font(list: &str) -> String {
    list.split(',').next().unwrap_or("").trim().trim_matches(['"', '\'']).to_string()
}

/// Every choice about a plot: which crosslinks, the proteins, the ring,
/// tracks, links and canvas. Missing keys of a JSON style take their default.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct CircosStyle {
    // which crosslinks
    /// Highest CSM and unique-crosslink q-value shown
    pub max_q: f32,
    /// Fewest crosslink-spectrum matches of a link shown
    pub min_csms: u32,
    /// Crosslinks between proteins
    pub inter: bool,
    /// Crosslinks between two peptides of one protein
    pub intra: bool,
    pub loops: bool,
    /// Each link of a multivalent species
    pub multivalent: bool,
    /// A peptide found in several proteins (or places) linked in each of
    /// them; otherwise only in its first
    pub all_places: bool,

    // the proteins
    /// Accessions in drawing order (clockwise); proteins not listed follow,
    /// most crosslinked first
    pub order: Vec<String>,
    /// Accessions left out (with their crosslinks)
    pub hidden: Vec<String>,
    /// Accession → label (default: gene, entry name or accession)
    pub labels: BTreeMap<String, String>,
    /// Accession → colour (`#rrggbb`; default: the palette's)
    pub colors: BTreeMap<String, String>,
    pub palette: Palette,

    // the ring
    /// Where the first protein starts, degrees clockwise from the top
    pub start_angle: f64,
    /// Space between proteins, degrees
    pub gap: f64,
    /// Outer radius of the protein ring, as a fraction of half the figure
    pub radius: f64,
    /// Thickness of the protein ring, as a fraction of half the figure
    pub thickness: f64,
    pub label_mode: LabelMode,
    pub ticks: bool,
    pub tick_labels: bool,
    /// Residues between labelled ticks (0: chosen from the longest protein)
    pub tick_spacing: u32,

    // the tracks, inside the ring from the outside in
    /// Residues covered by the linear search's peptides at `max_q`
    pub coverage: bool,
    /// A tab-separated file of regions (protein, start, end, label, colour)
    pub annotations: String,
    /// The residues the linker reacts with
    pub reactive: bool,
    /// Residues for `reactive` (empty: the search's crosslinker's)
    pub reactive_residues: String,
    /// A bar at each crosslinked residue, its height the number of CSMs
    pub sites: bool,
    /// Width of each track, as a fraction of half the figure
    pub track_size: f64,

    // the links
    pub color_by: ColorBy,
    pub inter_color: String,
    pub intra_color: String,
    pub loop_color: String,
    pub multi_color: String,
    pub single_color: String,
    /// Gradient ends (score, q-value, CSMs)
    pub low_color: String,
    pub high_color: String,
    /// Line width (points)
    pub link_width: f64,
    /// Wider lines for links seen in more CSMs
    pub width_by_csms: bool,
    pub link_opacity: f64,
    /// How far the curves bend towards the centre (Circos' `bezier_radius`):
    /// 0 through the centre, 1 straight
    pub bezier_radius: f64,
    /// Crosslinks picked out in colours of their own; a crosslink takes the
    /// colour of the first highlight that picks it
    pub highlights: Vec<Highlight>,
    /// Width of the highlighted crosslinks, times their own
    pub highlight_width: f64,
    /// The crosslinks no highlight picks in `dim_color`, fainter (while a
    /// highlight is drawn)
    pub dim_others: bool,
    pub dim_color: String,

    // the canvas
    /// Width and height, inches
    pub size: f64,
    pub background: Background,
    pub title: String,
    /// Font families in order of preference, as in CSS: the first one
    /// installed is used
    pub font: String,
    /// Size (points), bold and italic of the title, the protein labels, the
    /// tick numbers, the legend and the annotation labels
    pub title_text: TypeStyle,
    pub label_text: TypeStyle,
    pub tick_text: TypeStyle,
    pub legend_text: TypeStyle,
    pub region_text: TypeStyle,
    pub legend: bool,
    /// Where the legend goes; the circle keeps clear of it
    pub legend_corner: LegendCorner,
}

impl Default for CircosStyle {
    fn default() -> Self {
        CircosStyle {
            max_q: 0.01,
            min_csms: 1,
            inter: true,
            intra: true,
            loops: true,
            multivalent: true,
            all_places: false,
            order: Vec::new(),
            hidden: Vec::new(),
            labels: BTreeMap::new(),
            colors: BTreeMap::new(),
            palette: Palette::Tableau20,
            start_angle: 0.0,
            gap: 2.0,
            radius: 0.78,
            thickness: 0.045,
            label_mode: LabelMode::Parallel,
            ticks: true,
            tick_labels: true,
            tick_spacing: 0,
            coverage: true,
            annotations: String::new(),
            reactive: false,
            reactive_residues: String::new(),
            sites: true,
            track_size: 0.04,
            color_by: ColorBy::Class,
            inter_color: hex(Class::InterProtein.color()),
            intra_color: hex(Class::IntraProtein.color()),
            loop_color: hex(Class::Loop.color()),
            multi_color: hex(Class::Multivalent.color()),
            single_color: "#555555".into(),
            low_color: "#c6dbef".into(),
            high_color: "#08306b".into(),
            link_width: 1.2,
            width_by_csms: true,
            link_opacity: 0.8,
            bezier_radius: 0.15,
            size: 8.0,
            background: Background::White,
            title: String::new(),
            font: "Arial, Helvetica, sans-serif".into(),
            title_text: TypeStyle::new(15.0, true),
            label_text: TypeStyle::new(11.0, true),
            tick_text: TypeStyle::new(8.0, false),
            legend_text: TypeStyle::new(8.0, false),
            region_text: TypeStyle::new(7.0, false),
            legend: true,
            legend_corner: LegendCorner::BottomLeft,
            highlights: Vec::new(),
            highlight_width: 1.5,
            dim_others: false,
            dim_color: "#c8c8c8".into(),
        }
    }
}

/// Ready-made styles
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Preset {
    Default,
    /// Black and grey: links in dark grey, proteins in light grey
    Monochrome,
    /// White on black, for slides
    Dark,
    /// Links coloured by score, for reviewing a search
    Review,
}

impl Preset {
    pub const ALL: [Preset; 4] = [Preset::Default, Preset::Monochrome, Preset::Dark, Preset::Review];

    pub fn label(self) -> &'static str {
        match self {
            Preset::Default => "default",
            Preset::Monochrome => "monochrome",
            Preset::Dark => "dark (slides)",
            Preset::Review => "review (by score)",
        }
    }

    /// The preset's look, keeping `base`'s selection (q-value, proteins, labels, tracks' files)
    pub fn apply(self, base: &CircosStyle) -> CircosStyle {
        let d = CircosStyle::default();
        let keep = |mut s: CircosStyle| {
            s.max_q = base.max_q;
            s.min_csms = base.min_csms;
            s.inter = base.inter;
            s.intra = base.intra;
            s.loops = base.loops;
            s.multivalent = base.multivalent;
            s.all_places = base.all_places;
            s.order = base.order.clone();
            s.hidden = base.hidden.clone();
            s.labels = base.labels.clone();
            s.annotations = base.annotations.clone();
            s.title = base.title.clone();
            s.font = base.font.clone();
            s.title_text = base.title_text;
            s.label_text = base.label_text;
            s.tick_text = base.tick_text;
            s.legend_text = base.legend_text;
            s.region_text = base.region_text;
            s.legend = base.legend;
            s.legend_corner = base.legend_corner;
            s.highlights = base.highlights.clone();
            s.highlight_width = base.highlight_width;
            s.dim_others = base.dim_others;
            s.dim_color = base.dim_color.clone();
            s
        };
        keep(match self {
            Preset::Default => d,
            Preset::Monochrome => CircosStyle {
                palette: Palette::Grey,
                color_by: ColorBy::Single,
                single_color: "#333333".into(),
                link_opacity: 0.7,
                ..d
            },
            Preset::Dark => CircosStyle { background: Background::Black, link_opacity: 0.9, ..d },
            Preset::Review => CircosStyle { color_by: ColorBy::Score, width_by_csms: true, ..d },
        })
    }
}

impl CircosStyle {
    /// The parts that decide which crosslinks and proteins are in the plot
    /// (the data is rebuilt only when they change)
    pub fn selection(&self) -> (f32, u32, [bool; 5], Vec<String>, String, bool, String) {
        (
            self.max_q,
            self.min_csms,
            [self.inter, self.intra, self.loops, self.multivalent, self.all_places],
            self.hidden.clone(),
            self.annotations.clone(),
            self.coverage,
            self.reactive_residues.clone(),
        )
    }

    pub fn to_json(&self) -> String {
        serde_json::to_string_pretty(self).unwrap_or_default()
    }

    pub fn from_json(text: &str) -> Result<CircosStyle, String> {
        serde_json::from_str(text).map_err(|e| e.to_string())
    }

    /// Every text larger or smaller by `factor`
    pub fn scale_text(&mut self, factor: f64) {
        for t in [&mut self.title_text, &mut self.label_text, &mut self.tick_text, &mut self.legend_text, &mut self.region_text] {
            t.size = ((t.size * factor).clamp(2.0, 96.0) * 10.0).round() / 10.0;
        }
    }

    /// Take the settings that belong to one set of proteins from `from`:
    /// which are shown, their order, labels and colours, the annotation file
    /// and the title
    pub fn set_proteins(&mut self, from: &CircosStyle) {
        self.order = from.order.clone();
        self.hidden = from.hidden.clone();
        self.labels = from.labels.clone();
        self.colors = from.colors.clone();
        self.annotations = from.annotations.clone();
        self.title = from.title.clone();
        // highlights of proteins too; those of a kind of crosslink stay
        let kinds: Vec<Highlight> = self.highlights.iter().filter(|h| !h.pick.of_proteins()).cloned().collect();
        self.highlights = from.highlights.iter().filter(|h| h.pick.of_proteins()).cloned().chain(kinds).collect();
    }
}

/// A style saved under a name in sage-gui's library, to give other results
/// the same look: the crosslinks drawn, ring, tracks, links, colours, sizes
/// and fonts, and, if asked, the proteins
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Template {
    pub name: String,
    /// It also sets the proteins (matched by accession): which are shown,
    /// their order, labels and colours, the annotation file and the title.
    /// For other results of the same proteins.
    #[serde(default)]
    pub proteins: bool,
    #[serde(default)]
    pub style: CircosStyle,
}

impl Template {
    /// `style` saved as `name`. With `proteins`, the proteins of `data` keep
    /// the place and colour they have now, palette colours included, so they
    /// look the same in other results; without, the protein settings are
    /// left out.
    pub fn new(name: &str, style: &CircosStyle, proteins: bool, data: &CircosData) -> Template {
        let mut style = style.clone();
        if proteins {
            let shown: Vec<(String, String)> =
                shown_proteins(data, &style).into_iter().map(|(p, _, color)| (p.accession.clone(), color)).collect();
            let mut order: Vec<String> = shown.iter().map(|(a, _)| a.clone()).collect();
            let palette = style.palette;
            for (i, p) in data.proteins.iter().enumerate() {
                if !order.contains(&p.accession) {
                    order.push(p.accession.clone());
                }
                style.colors.entry(p.accession.clone()).or_insert_with(|| palette.color(i).to_string());
            }
            for (accession, color) in shown {
                style.colors.entry(accession).or_insert(color);
            }
            // proteins of earlier results, not in these, keep their place after them
            for a in &style.order {
                if !order.contains(a) {
                    order.push(a.clone());
                }
            }
            style.order = order;
        } else {
            style.set_proteins(&CircosStyle::default());
        }
        Template { name: name.trim().to_string(), proteins, style }
    }

    /// `current` in this template's style. A template without proteins
    /// leaves the current proteins, labels, colours, annotations and title.
    pub fn apply(&self, current: &CircosStyle) -> CircosStyle {
        let mut style = self.style.clone();
        if !self.proteins {
            style.set_proteins(current);
        }
        style
    }

    /// What it sets, in a line
    pub fn summary(&self) -> String {
        let s = &self.style;
        let font = first_font(&s.font);
        let mut tracks: Vec<&str> = Vec::new();
        for (on, name) in [(s.sites, "sites"), (s.coverage, "coverage"), (s.reactive, "reactive residues")] {
            if on {
                tracks.push(name);
            }
        }
        if !s.annotations.trim().is_empty() {
            tracks.push("annotations");
        }
        let mut parts = vec![
            format!("{} in, {font} {} pt, {} background", s.size, s.label_text.size, s.background.label()),
            format!("links by {}", s.color_by.label()),
            format!("{} palette", s.palette.label()),
            format!("q ≤ {}, ≥ {} CSM{}", s.max_q, s.min_csms, if s.min_csms == 1 { "" } else { "s" }),
            if tracks.is_empty() { "no tracks".to_string() } else { format!("tracks: {}", tracks.join(", ")) },
        ];
        if self.proteins {
            parts.push(format!("{} proteins with their order, labels and colours", s.order.len()));
        }
        parts.join(" · ")
    }
}

// ---------------------------------------------------------------- the data

/// A protein of the plot
#[derive(Clone, Debug, PartialEq)]
pub struct Segment {
    pub accession: String,
    pub label: String,
    pub sequence: String,
    /// Crosslinks touching it
    pub links: usize,
}

/// A crosslink between two residues (1-based positions in their proteins),
/// with every CSM of it at the style's q-value
#[derive(Clone, Debug, PartialEq)]
pub struct Link {
    pub a: (String, usize),
    pub b: (String, usize),
    pub class: Class,
    pub csms: u32,
    /// Lowest crosslink q-value of its CSMs
    pub q: f32,
    /// Highest hyperscore of its CSMs
    pub score: f64,
}

/// A region of an annotation file
#[derive(Clone, Debug, PartialEq)]
pub struct Region {
    pub accession: String,
    pub start: usize,
    pub end: usize,
    pub label: String,
    pub color: String,
}

/// What a plot shows, built from a results folder
#[derive(Clone, Debug, Default, PartialEq)]
pub struct CircosData {
    /// Every protein with a crosslink (hidden ones too), most crosslinked first
    pub proteins: Vec<Segment>,
    pub links: Vec<Link>,
    /// Accession → covered residue ranges (1-based, inclusive) of the linear search
    pub coverage: HashMap<String, Vec<(usize, usize)>>,
    pub regions: Vec<Region>,
    /// Residues the linker reacts with (`K`, `KSTY`)
    pub reactive: String,
    /// Annotation file problems, for the user
    pub notes: Vec<String>,
}

/// A peptide in Sage's notation without its modifications
fn bare(peptide: &str) -> String {
    let mut out = String::with_capacity(peptide.len());
    let mut depth = 0;
    for c in peptide.chars() {
        match c {
            '[' => depth += 1,
            ']' => depth -= 1,
            c if depth == 0 && c.is_ascii_alphabetic() => out.push(c),
            _ => {}
        }
    }
    out
}

/// Where a peptide sits: `(accession, 0-based start)` in the proteins Sage
/// listed for it (else any FASTA protein holding it); only the first unless `all`
fn places(results: &Results, sequence: &str, listed: Option<&String>, all: bool) -> Vec<(String, usize)> {
    let mut accessions: Vec<String> = listed
        .map(|l| l.split(';').map(str::trim).filter(|a| !a.is_empty()).map(String::from).collect())
        .unwrap_or_default();
    if accessions.is_empty() {
        accessions = results.proteins().filter(|(_, p)| p.sequence.contains(sequence)).map(|(a, _)| a.clone()).collect();
        accessions.sort();
    }
    let mut out = Vec::new();
    for accession in accessions {
        if accession.starts_with(results.decoy_tag()) {
            continue;
        }
        let Some(protein) = results.protein(&accession) else { continue };
        for (start, _) in protein.sequence.match_indices(sequence) {
            out.push((accession.clone(), start));
            if !all {
                return out;
            }
        }
    }
    out
}

/// A protein's default label: gene, else entry name (`RPB1_YEAST`), else accession
fn default_label(results: &Results, accession: &str) -> String {
    match results.protein(accession) {
        Some(p) if !p.gene.is_empty() => p.gene.clone(),
        _ => accession.rsplit('|').next().unwrap_or(accession).to_string(),
    }
}

/// The crosslinks of a match as pairs of `(peptide, 0-based residue)` and
/// the class of each pair's kind
fn ends_of(m: &Match) -> Vec<((String, Option<String>, usize), (String, Option<String>, usize))> {
    let list = |i: usize| m.protein_lists.get(i).cloned().filter(|l| !l.is_empty());
    let site = |s: usize| s.saturating_sub(1);
    match m.kind {
        Kind::Linear => Vec::new(),
        Kind::Crosslink => vec![((bare(&m.peptide), list(0), site(m.sites.0)), (bare(&m.partner), list(1), site(m.sites.1)))],
        Kind::Loop => {
            let p = bare(&m.peptide);
            vec![((p.clone(), list(0), site(m.sites.0)), (p, list(0), site(m.sites.1)))]
        }
        Kind::Multi => {
            let (peptides, links) = parse_species(&m.species);
            let peptides: Vec<String> = peptides.iter().map(|p| bare(p)).collect();
            // the protein lists Sage wrote belong to the species' first two peptides
            let listed = |k: usize| -> Option<String> {
                let sequence = peptides.get(k)?;
                [&m.peptide, &m.partner].iter().position(|p| &bare(p) == sequence).and_then(list)
            };
            links
                .into_iter()
                .filter_map(|(pa, ra, pb, rb)| {
                    Some(((peptides.get(pa)?.clone(), listed(pa), ra), (peptides.get(pb)?.clone(), listed(pb), rb)))
                })
                .collect()
        }
    }
}

/// The plot's data from a loaded results folder
pub fn build(results: &Results, style: &CircosStyle) -> CircosData {
    let mut merged: BTreeMap<((String, usize), (String, usize), u8), Link> = BTreeMap::new();
    let class_rank = |c: Class| c as u8;
    for m in &results.matches {
        if !m.target || m.rank > 1 || m.spectrum_q > style.max_q || m.q > style.max_q {
            continue;
        }
        let wanted = match m.kind {
            Kind::Linear => continue,
            Kind::Loop => style.loops,
            Kind::Multi => style.multivalent,
            Kind::Crosslink => style.inter || style.intra,
        };
        if !wanted {
            continue;
        }
        for ((pa, la, ra), (pb, lb, rb)) in ends_of(m) {
            let here = places(results, &pa, la.as_ref(), style.all_places);
            let there = places(results, &pb, lb.as_ref(), style.all_places);
            for (acc_a, start_a) in &here {
                for (acc_b, start_b) in &there {
                    let a = (acc_a.clone(), start_a + ra + 1);
                    let b = (acc_b.clone(), start_b + rb + 1);
                    let class = match m.kind {
                        Kind::Multi => Class::Multivalent,
                        Kind::Loop => Class::Loop,
                        _ if acc_a == acc_b => Class::IntraProtein,
                        _ => Class::InterProtein,
                    };
                    let shown = match class {
                        Class::InterProtein => style.inter,
                        Class::IntraProtein => style.intra,
                        _ => true,
                    };
                    if !shown || (a == b && class != Class::Multivalent) {
                        continue;
                    }
                    let (a, b) = if a <= b { (a, b) } else { (b, a) };
                    let link = merged.entry((a.clone(), b.clone(), class_rank(class))).or_insert(Link {
                        a,
                        b,
                        class,
                        csms: 0,
                        q: f32::MAX,
                        score: f64::MIN,
                    });
                    link.csms += 1;
                    link.q = link.q.min(m.q);
                    link.score = link.score.max(m.score);
                }
            }
        }
    }
    let links: Vec<Link> = merged
        .into_values()
        .filter(|l| l.csms >= style.min_csms)
        .filter(|l| !style.hidden.contains(&l.a.0) && !style.hidden.contains(&l.b.0))
        .collect();

    // the proteins of the links (and of hidden ones, for the protein list)
    let mut counts: HashMap<String, usize> = HashMap::new();
    for l in &links {
        *counts.entry(l.a.0.clone()).or_default() += 1;
        if l.b.0 != l.a.0 {
            *counts.entry(l.b.0.clone()).or_default() += 1;
        }
    }
    for h in &style.hidden {
        counts.entry(h.clone()).or_default();
    }
    let mut proteins: Vec<Segment> = counts
        .into_iter()
        .filter_map(|(accession, links)| {
            let protein = results.protein(&accession)?;
            Some(Segment { label: default_label(results, &accession), sequence: protein.sequence.clone(), accession, links })
        })
        .collect();
    proteins.sort_by(|a, b| b.links.cmp(&a.links).then(a.label.cmp(&b.label)));

    // linear sequence coverage
    let mut coverage: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
    if style.coverage {
        let shown: std::collections::HashSet<&str> = proteins.iter().map(|p| p.accession.as_str()).collect();
        for m in results.matches.iter().filter(|m| m.kind == Kind::Linear && m.target && m.q <= style.max_q && m.spectrum_q <= style.max_q) {
            let sequence = bare(&m.peptide);
            for (accession, start) in places(results, &sequence, m.protein_lists.first(), true) {
                if shown.contains(accession.as_str()) {
                    coverage.entry(accession).or_default().push((start + 1, start + sequence.len()));
                }
            }
        }
        for ranges in coverage.values_mut() {
            ranges.sort_unstable();
            let mut merged: Vec<(usize, usize)> = Vec::new();
            for &(s, e) in ranges.iter() {
                match merged.last_mut() {
                    Some(last) if s <= last.1 + 1 => last.1 = last.1.max(e),
                    _ => merged.push((s, e)),
                }
            }
            *ranges = merged;
        }
    }

    let reactive = if style.reactive_residues.trim().is_empty() {
        results
            .crosslink_settings()
            .map(|s| s.site1.clone().unwrap_or_else(|| s.residues.clone()) + s.site2.as_deref().unwrap_or(""))
            .unwrap_or_else(|| "K".into())
    } else {
        style.reactive_residues.clone()
    };
    let reactive: String = reactive.chars().filter(|c| c.is_ascii_uppercase()).collect();

    let (regions, notes) = match style.annotations.trim() {
        "" => (Vec::new(), Vec::new()),
        path => read_regions(Path::new(path), &proteins, results),
    };
    CircosData { proteins, links, coverage, regions, reactive, notes }
}

/// Regions from a tab-separated file: protein (accession, gene, entry name
/// or label), start, end, label and an optional `#rrggbb` colour; lines
/// starting with `#` and a header line (not a number in column 2) are skipped
fn read_regions(path: &Path, proteins: &[Segment], results: &Results) -> (Vec<Region>, Vec<String>) {
    let text = match std::fs::read_to_string(path) {
        Ok(t) => t,
        Err(e) => return (Vec::new(), vec![format!("annotation file {}: {e}", path.display())]),
    };
    let find = |name: &str| -> Option<String> {
        let name = name.trim();
        proteins
            .iter()
            .find(|p| {
                p.accession == name
                    || p.label.eq_ignore_ascii_case(name)
                    || p.accession.split('|').any(|part| part.eq_ignore_ascii_case(name))
                    || results.protein(&p.accession).is_some_and(|q| q.gene.eq_ignore_ascii_case(name))
            })
            .map(|p| p.accession.clone())
    };
    let mut regions = Vec::new();
    let mut unknown: Vec<String> = Vec::new();
    for line in text.lines().map(str::trim).filter(|l| !l.is_empty() && !l.starts_with('#')) {
        let cols: Vec<&str> = line.split('\t').collect();
        let (Some(protein), Some(start), Some(end)) = (cols.first(), cols.get(1), cols.get(2)) else { continue };
        let (Ok(start), Ok(end)) = (start.trim().parse::<usize>(), end.trim().parse::<usize>()) else { continue };
        let Some(accession) = find(protein) else {
            if !unknown.contains(&protein.to_string()) {
                unknown.push(protein.to_string());
            }
            continue;
        };
        let color = cols.get(4).map(|c| c.trim()).filter(|c| c.starts_with('#') && c.len() == 7).unwrap_or("").to_string();
        regions.push(Region {
            accession,
            start: start.min(end),
            end: start.max(end),
            label: cols.get(3).map(|l| l.trim().to_string()).unwrap_or_default(),
            color,
        });
    }
    let notes = match unknown.len() {
        0 => Vec::new(),
        _ => vec![format!("annotation file: not among the plotted proteins: {}", unknown.join(", "))],
    };
    (regions, notes)
}

// ---------------------------------------------------------------- drawing

fn hex(c: eframe::egui::Color32) -> String {
    format!("#{:02x}{:02x}{:02x}", c.r(), c.g(), c.b())
}

/// `#rrggbb` → RGB (grey for anything else)
pub fn rgb(color: &str) -> (u8, u8, u8) {
    let c = color.trim().trim_start_matches('#');
    let part = |i: usize| c.get(i..i + 2).and_then(|h| u8::from_str_radix(h, 16).ok());
    match (c.len(), part(0), part(2), part(4)) {
        (6, Some(r), Some(g), Some(b)) => (r, g, b),
        _ => (128, 128, 128),
    }
}

fn mix(low: &str, high: &str, t: f64) -> String {
    let (a, b) = (rgb(low), rgb(high));
    let t = t.clamp(0.0, 1.0);
    let m = |x: u8, y: u8| (x as f64 + (y as f64 - x as f64) * t).round() as u8;
    format!("#{:02x}{:02x}{:02x}", m(a.0, b.0), m(a.1, b.1), m(a.2, b.2))
}

fn esc(s: &str) -> String {
    s.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;").replace('"', "&quot;")
}

/// Points → SVG units
fn pt(points: f64) -> f64 {
    points / 72.0 * UNITS_PER_IN
}

/// The proteins drawn, in order (the style's order first), with their colours
pub fn shown_proteins<'a>(data: &'a CircosData, style: &CircosStyle) -> Vec<(&'a Segment, String, String)> {
    let mut shown: Vec<&Segment> = data.proteins.iter().filter(|p| !style.hidden.contains(&p.accession)).collect();
    shown.sort_by_key(|p| style.order.iter().position(|a| a == &p.accession).unwrap_or(usize::MAX));
    // palette colours follow the order of every protein, so hiding one keeps the others' colours
    shown
        .into_iter()
        .map(|p| {
            let i = data.proteins.iter().position(|q| q.accession == p.accession).unwrap_or(0);
            let color = style.colors.get(&p.accession).cloned().unwrap_or_else(|| style.palette.color(i).to_string());
            let label = style.labels.get(&p.accession).cloned().filter(|l| !l.trim().is_empty()).unwrap_or_else(|| p.label.clone());
            (p, label, color)
        })
        .collect()
}

/// Residues between ticks: labelled (major) and plain (minor)
fn tick_steps(longest: usize, spacing: u32) -> (usize, usize) {
    let major = if spacing > 0 { spacing as usize } else { figure::tick_step(longest as f64, 5.0).max(1.0) as usize };
    let minor = if major % 5 == 0 { major / 5 } else if major % 2 == 0 { major / 2 } else { major };
    (major.max(1), minor.max(1))
}

/// The plot as SVG
pub fn svg(data: &CircosData, style: &CircosStyle) -> String {
    let size = style.size.clamp(2.0, 40.0);
    let s = size * UNITS_PER_IN;
    let half = s / 2.0;
    let title_space = if style.title.trim().is_empty() { 0.0 } else { style.title_text.units() * 1.8 };
    let (cx0, cy0) = (half, half + title_space * 0.5);
    let ink = if style.background == Background::Black { "#e8e8e8" } else { "#222222" };
    let faint = if style.background == Background::Black { "#666666" } else { "#bbbbbb" };
    let font = esc(&style.font);
    let tick_size = style.tick_text.units();
    let legend_size = style.legend_text.units();
    let region_size = style.region_text.units();

    let track = style.track_size.clamp(0.005, 0.2) * half;
    let pad = half * 0.008;

    let proteins = shown_proteins(data, style);
    let mut out = String::new();
    let _ = write!(
        out,
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{size}in" height="{size}in" viewBox="0 0 {s:.1} {s:.1}">"#
    );
    match style.background {
        Background::White => {
            let _ = write!(out, r##"<rect width="{s:.1}" height="{s:.1}" fill="#ffffff"/>"##);
        }
        Background::Black => {
            let _ = write!(out, r##"<rect width="{s:.1}" height="{s:.1}" fill="#000000"/>"##);
        }
        Background::Transparent => {}
    }
    if title_space > 0.0 {
        let _ = write!(
            out,
            r#"<text x="{half:.1}" y="{:.1}" text-anchor="middle" {} fill="{ink}">{}</text>"#,
            title_space * 0.75,
            style.title_text.attrs(&font),
            esc(&style.title)
        );
    }
    if proteins.is_empty() {
        let _ = write!(
            out,
            r#"<text x="{half:.1}" y="{cy0:.1}" text-anchor="middle" {} fill="{ink}">No crosslinks at q ≤ {}</text></svg>"#,
            style.label_text.attrs(&font),
            style.max_q
        );
        return out;
    }

    // the crosslinks drawn (both proteins shown) and their range
    let shown: std::collections::HashSet<&str> = proteins.iter().map(|(p, _, _)| p.accession.as_str()).collect();
    let shown_links: Vec<&Link> = data.links.iter().filter(|l| shown.contains(l.a.0.as_str()) && shown.contains(l.b.0.as_str())).collect();
    let (score_lo, score_hi) = shown_links.iter().fold((f64::MAX, f64::MIN), |(lo, hi), l| (lo.min(l.score), hi.max(l.score)));
    let most_csms = shown_links.iter().map(|l| l.csms).max().unwrap_or(1).max(1);
    let (legend_rows, legend_summary) = legend_rows(data, style, &shown_links, (score_lo, score_hi), most_csms);

    // how far ticks and labels reach beyond the ring
    let longest = proteins.iter().map(|(p, _, _)| p.sequence.len()).max().unwrap_or(1);
    let (major, minor) = tick_steps(longest, style.tick_spacing);
    let tick_len = half * 0.018;
    let tick_reach = match (style.ticks, style.tick_labels) {
        (false, _) => pad * 2.0,
        (true, false) => tick_len + pad,
        (true, true) => tick_len + pad + tick_size * (major.max(longest).to_string().len() as f64) * 0.62 + pad,
    };
    let label_size = style.label_text.units();
    let label_reach = match style.label_mode {
        LabelMode::None => 0.0,
        LabelMode::Parallel => pad * 2.0 + label_size * 1.1,
        LabelMode::Radial => {
            let longest_label = proteins.iter().map(|(_, label, _)| label.chars().count()).max().unwrap_or(0);
            pad * 2.0 + label_size * char_width(&style.label_text) * longest_label as f64
        }
    };
    // the legend's box, in its corner: the circle keeps clear of it
    let legend_line = legend_size * 1.45;
    let legend_box = style.legend.then(|| {
        let em = legend_size * char_width(&style.legend_text);
        let widest = legend_rows
            .iter()
            .map(|(_, text)| legend_size * 2.2 + em * text.chars().count() as f64)
            .fold(em * legend_summary.chars().count() as f64, f64::max);
        let height = legend_size * 1.1 + legend_line * legend_rows.len() as f64;
        let margin = s * 0.03;
        let x0 = if style.legend_corner.right() { (s - margin - widest).max(0.0) } else { margin };
        let y0 = if style.legend_corner.top() { margin.max(title_space + pad) } else { (s - margin - height).max(0.0) };
        [x0, y0, x0 + widest, y0 + height]
    });
    let r1_wanted = (style.radius.clamp(0.2, 0.95) * half - title_space * 0.5).max(half * 0.2);
    let gap = pad * 3.0 + legend_size * 0.3;
    let (cx, cy, r1) = fit_circle(s, title_space, (cx0, cy0), r1_wanted, half * 0.2, tick_reach + label_reach, legend_box, gap);
    let r0 = (r1 - style.thickness.clamp(0.005, 0.3) * half).max(r1 * 0.5);

    // angles: degrees clockwise from the top
    let total: usize = proteins.iter().map(|(p, _, _)| p.sequence.len().max(1)).sum();
    let n = proteins.len() as f64;
    let gap = style.gap.clamp(0.0, 30.0);
    let available = (360.0 - gap * n).max(10.0);
    let per = available / total as f64;
    let mut starts: HashMap<&str, (f64, f64)> = HashMap::new();
    let mut angle = style.start_angle + gap / 2.0;
    for (p, _, _) in &proteins {
        let span = p.sequence.len().max(1) as f64 * per;
        starts.insert(p.accession.as_str(), (angle, span));
        angle += span + gap;
    }
    let at = |accession: &str, residue: usize| -> Option<f64> {
        starts.get(accession).map(|&(a, _)| a + (residue as f64 - 0.5) * per)
    };
    let xy = |deg: f64, r: f64| {
        let t = deg.to_radians();
        (cx + r * t.sin(), cy - r * t.cos())
    };
    let band = |out: &mut String, a0: f64, a1: f64, ra: f64, rb: f64, attrs: &str| {
        let (x0, y0) = xy(a0, rb);
        let (x1, y1) = xy(a1, rb);
        let (x2, y2) = xy(a1, ra);
        let (x3, y3) = xy(a0, ra);
        let large = if a1 - a0 > 180.0 { 1 } else { 0 };
        let _ = write!(
            out,
            r#"<path d="M{x0:.2},{y0:.2} A{rb:.2},{rb:.2} 0 {large} 1 {x1:.2},{y1:.2} L{x2:.2},{y2:.2} A{ra:.2},{ra:.2} 0 {large} 0 {x3:.2},{y3:.2} Z" {attrs}/>"#
        );
    };
    // text along the arc (rotated to read from outside) or radial
    let text_at = |out: &mut String, deg: f64, r: f64, text: &str, look: &TypeStyle, radial: bool| {
        let (x, y) = xy(deg, r);
        let a = deg.rem_euclid(360.0);
        let attrs = look.attrs(&font);
        let (rotate, anchor, baseline) = if radial {
            if a < 180.0 { (a - 90.0, "start", "middle") } else { (a + 90.0, "end", "middle") }
        } else if a > 90.0 && a < 270.0 {
            (a + 180.0, "middle", "hanging")
        } else {
            (a, "middle", "auto")
        };
        let _ = write!(
            out,
            r#"<text x="{x:.2}" y="{y:.2}" transform="rotate({rotate:.2} {x:.2} {y:.2})" text-anchor="{anchor}" dominant-baseline="{baseline}" {attrs} fill="{ink}">{}</text>"#,
            esc(text)
        );
    };

    // the protein ring
    out.push_str("<g>");
    for (p, _, color) in &proteins {
        let (a0, span) = starts[p.accession.as_str()];
        band(&mut out, a0, a0 + span, r0, r1, &format!(r#"fill="{color}" stroke="{ink}" stroke-width="0.6""#));
    }
    out.push_str("</g>");

    // ticks and their numbers
    let mut label_r = r1 + pad * 2.0;
    if style.ticks {
        let _ = write!(out, r#"<g stroke="{ink}" stroke-width="0.7">"#);
        for (p, _, _) in &proteins {
            let len = p.sequence.len();
            let mut pos = minor;
            while pos <= len {
                let a = at(&p.accession, pos).unwrap_or_default();
                let l = if pos % major == 0 { tick_len } else { tick_len * 0.5 };
                let (x0, y0) = xy(a, r1);
                let (x1, y1) = xy(a, r1 + l);
                let _ = write!(out, r#"<line x1="{x0:.2}" y1="{y0:.2}" x2="{x1:.2}" y2="{y1:.2}"/>"#);
                pos += minor;
            }
        }
        out.push_str("</g>");
        label_r = r1 + tick_len + pad;
        if style.tick_labels {
            for (p, _, _) in &proteins {
                let mut pos = major;
                while pos <= p.sequence.len() {
                    let a = at(&p.accession, pos).unwrap_or_default();
                    text_at(&mut out, a, r1 + tick_len + pad, &pos.to_string(), &style.tick_text, true);
                    pos += major;
                }
            }
            // the numbers' length, radially
            label_r += tick_size * (major.max(longest).to_string().len() as f64) * 0.62 + pad;
        }
    }
    // protein labels
    if style.label_mode != LabelMode::None {
        for (p, label, _) in &proteins {
            let (a0, span) = starts[p.accession.as_str()];
            let radial = style.label_mode == LabelMode::Radial;
            text_at(&mut out, a0 + span / 2.0, label_r + pad * 2.0, label, &style.label_text, radial);
        }
    }

    // tracks inside the ring, from the outside in
    let mut r = r0 - pad * 2.0;
    if style.coverage && !data.coverage.is_empty() {
        let (ra, rb) = (r - track * 0.6, r);
        let _ = write!(out, r#"<g fill="{faint}">"#);
        for (p, _, _) in &proteins {
            let (a0, span) = starts[p.accession.as_str()];
            band(&mut out, a0, a0 + span, ra, rb, r#"fill-opacity="0.25""#);
            for &(s0, s1) in data.coverage.get(&p.accession).map(Vec::as_slice).unwrap_or(&[]) {
                let a = at(&p.accession, s0).unwrap_or(a0) - per / 2.0;
                let b = at(&p.accession, s1).unwrap_or(a0) + per / 2.0;
                band(&mut out, a, b, ra, rb, r#"fill-opacity="0.9""#);
            }
        }
        out.push_str("</g>");
        r = ra - pad * 2.0;
    }
    let regions: Vec<&Region> = data.regions.iter().filter(|g| starts.contains_key(g.accession.as_str())).collect();
    if !regions.is_empty() {
        let (ra, rb) = (r - track, r);
        for (i, g) in regions.iter().enumerate() {
            let color = if g.color.is_empty() { Palette::OkabeIto.color(i).to_string() } else { g.color.clone() };
            let a = at(&g.accession, g.start).unwrap_or_default() - per / 2.0;
            let b = at(&g.accession, g.end).unwrap_or_default() + per / 2.0;
            band(&mut out, a, b, ra, rb, &format!(r#"fill="{color}" fill-opacity="0.85""#));
            // its label, where it fits
            let arc = (b - a).to_radians() * (ra + rb) / 2.0;
            if !g.label.is_empty() && arc > region_size * 0.66 * g.label.chars().count() as f64 {
                text_at(&mut out, (a + b) / 2.0, (ra + rb) / 2.0, &g.label, &style.region_text, false);
            }
        }
        r = ra - pad * 2.0;
    }
    if style.reactive && !data.reactive.is_empty() {
        let (ra, rb) = (r - track * 0.5, r);
        let _ = write!(out, r#"<g stroke="{ink}" stroke-width="0.5" stroke-opacity="0.6">"#);
        for (p, _, _) in &proteins {
            for (i, residue) in p.sequence.chars().enumerate() {
                if data.reactive.contains(residue) {
                    let a = at(&p.accession, i + 1).unwrap_or_default();
                    let (x0, y0) = xy(a, ra);
                    let (x1, y1) = xy(a, rb);
                    let _ = write!(out, r#"<line x1="{x0:.2}" y1="{y0:.2}" x2="{x1:.2}" y2="{y1:.2}"/>"#);
                }
            }
        }
        out.push_str("</g>");
        r = ra - pad * 2.0;
    }
    if style.sites && !shown_links.is_empty() {
        let mut per_site: HashMap<(&str, usize), u32> = HashMap::new();
        for l in &shown_links {
            for (acc, pos) in [&l.a, &l.b] {
                *per_site.entry((acc.as_str(), *pos)).or_default() += l.csms;
            }
        }
        let most = per_site.values().copied().max().unwrap_or(1).max(1) as f64;
        let colors: HashMap<&str, &str> = proteins.iter().map(|(p, _, c)| (p.accession.as_str(), c.as_str())).collect();
        let (ra, rb) = (r - track, r);
        let width = (per * 0.9).max(0.25);
        let mut sites: Vec<_> = per_site.into_iter().collect();
        sites.sort();
        for ((acc, pos), csms) in sites {
            let a = at(acc, pos).unwrap_or_default();
            let h = (ra + (rb - ra) * (csms as f64 / most).sqrt()).max(ra + pad * 0.5);
            band(&mut out, a - width / 2.0, a + width / 2.0, ra, h, &format!(r#"fill="{}" stroke="{ink}" stroke-width="0.3""#, colors.get(acc).unwrap_or(&"#888888")));
        }
        r = ra - pad * 2.0;
    }

    // the links
    let rl = r.max(r1 * 0.3);
    let protein_color: HashMap<&str, &str> = proteins.iter().map(|(p, _, c)| (p.accession.as_str(), c.as_str())).collect();
    let color_of = |l: &Link| -> String {
        match style.color_by {
            ColorBy::Class => match l.class {
                Class::InterProtein => style.inter_color.clone(),
                Class::IntraProtein => style.intra_color.clone(),
                Class::Loop => style.loop_color.clone(),
                _ => style.multi_color.clone(),
            },
            ColorBy::Protein => protein_color.get(l.a.0.as_str()).unwrap_or(&"#888888").to_string(),
            ColorBy::Single => style.single_color.clone(),
            ColorBy::Score => {
                let t = if score_hi > score_lo { (l.score - score_lo) / (score_hi - score_lo) } else { 1.0 };
                mix(&style.low_color, &style.high_color, t)
            }
            ColorBy::QValue => {
                let t = if style.max_q > 0.0 { 1.0 - (l.q / style.max_q) as f64 } else { 1.0 };
                mix(&style.low_color, &style.high_color, t)
            }
            ColorBy::Csms => mix(&style.low_color, &style.high_color, (l.csms as f64).ln() / (most_csms as f64).ln().max(1e-9)),
        }
    };
    let mut drawn: Vec<(&Link, Option<usize>)> = shown_links.iter().map(|l| (*l, highlight_of(style, l))).collect();
    let dim = style.dim_others && drawn.iter().any(|(_, h)| h.is_some());
    // the highlighted on top, and the strongest: weakest first
    drawn.sort_by(|(a, ha), (b, hb)| (ha.is_some(), a.csms, a.score).partial_cmp(&(hb.is_some(), b.csms, b.score)).unwrap_or(std::cmp::Ordering::Equal));
    let base_width = pt(style.link_width.max(0.05));
    let opacity = style.link_opacity.clamp(0.05, 1.0);
    let _ = write!(out, r#"<g fill="none" stroke-linecap="round" stroke-opacity="{opacity:.2}">"#);
    for &(l, highlight) in &drawn {
        let (Some(a), Some(b)) = (at(&l.a.0, l.a.1), at(&l.b.0, l.b.1)) else { continue };
        let (x0, y0) = xy(a, rl);
        let (x1, y1) = xy(b, rl);
        // the control point on the bisector of the shorter arc; close ends bend less
        let mut d = (b - a).rem_euclid(360.0);
        let mut mid = a + d / 2.0;
        if d > 180.0 {
            d = 360.0 - d;
            mid += 180.0;
        }
        let bend = style.bezier_radius.clamp(0.0, 1.0);
        let rc = rl * (bend + (1.0 - bend) * (1.0 - d / 180.0).powi(2));
        let (qx, qy) = xy(mid, rc);
        let mut width = if style.width_by_csms { base_width * (1.0 + (l.csms as f64).log2() * 0.6) } else { base_width };
        let (color, faint) = match highlight {
            Some(i) => {
                width *= style.highlight_width.clamp(0.2, 10.0);
                (style.highlights[i].color.clone(), format!(r#" stroke-opacity="{:.2}""#, opacity.max(0.9)))
            }
            None if dim => (style.dim_color.clone(), format!(r#" stroke-opacity="{:.2}""#, opacity * 0.6)),
            None => (color_of(l), String::new()),
        };
        let _ = write!(
            out,
            r#"<path d="M{x0:.2},{y0:.2} Q{qx:.2},{qy:.2} {x1:.2},{y1:.2}" stroke="{color}" stroke-width="{width:.2}"{faint}/>"#
        );
    }
    out.push_str("</g>");

    // the legend, in its box: a row per colour, then the summary
    if let Some([x0, y0, _, _]) = legend_box {
        let _ = write!(out, r#"<g {} fill="{ink}">"#, style.legend_text.attrs(&font));
        let mut y = y0 + legend_size * 0.85;
        for (color, text) in &legend_rows {
            let _ = write!(
                out,
                r#"<line x1="{x0:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{color}" stroke-width="{:.1}" stroke-linecap="round"/>"#,
                y - legend_size * 0.35,
                x0 + legend_size * 1.6,
                y - legend_size * 0.35,
                pt(1.6)
            );
            let _ = write!(out, r#"<text x="{:.1}" y="{y:.1}">{}</text>"#, x0 + legend_size * 2.2, esc(text));
            y += legend_line;
        }
        let _ = write!(out, r#"<text x="{x0:.1}" y="{y:.1}">{}</text>"#, esc(&legend_summary));
        out.push_str("</g>");
    }
    out.push_str("</svg>");
    out
}

/// The legend: a colour and a text per row (the highlights, then the other
/// crosslinks), and the summary line under them
fn legend_rows(
    data: &CircosData,
    style: &CircosStyle,
    all: &[&Link],
    (score_lo, score_hi): (f64, f64),
    most_csms: u32,
) -> (Vec<(String, String)>, String) {
    let mut rows: Vec<(String, String)> = Vec::new();
    let marks: Vec<Option<usize>> = all.iter().map(|l| highlight_of(style, l)).collect();
    for (i, h) in style.highlights.iter().enumerate() {
        let count = marks.iter().filter(|m| **m == Some(i)).count();
        if count > 0 {
            rows.push((h.color.clone(), format!("{} ({count})", h.text(data, style))));
        }
    }
    let links: Vec<&Link> = all.iter().zip(&marks).filter(|(_, m)| m.is_none()).map(|(l, _)| *l).collect();
    if style.dim_others && marks.iter().any(Option::is_some) {
        if !links.is_empty() {
            rows.push((style.dim_color.clone(), format!("other crosslinks ({})", links.len())));
        }
        let csms: u32 = all.iter().map(|l| l.csms).sum();
        return (rows, format!("{} crosslinks, {csms} CSMs, q ≤ {}", all.len(), style.max_q));
    }
    match style.color_by {
        _ if links.is_empty() => {}
        ColorBy::Class => {
            for (class, color, name) in [
                (Class::InterProtein, &style.inter_color, "between proteins"),
                (Class::IntraProtein, &style.intra_color, "within a protein"),
                (Class::Loop, &style.loop_color, "loop-links"),
                (Class::Multivalent, &style.multi_color, "multivalent"),
            ] {
                let count = links.iter().filter(|l| l.class == class).count();
                if count > 0 {
                    rows.push((color.clone(), format!("{name} ({count})")));
                }
            }
        }
        ColorBy::Single | ColorBy::Protein => {}
        gradient => {
            let (lo, hi) = match gradient {
                ColorBy::Score => (format!("{score_lo:.1}"), format!("{score_hi:.1}")),
                ColorBy::QValue => (format!("{}", style.max_q), "0".into()),
                _ => ("1".into(), most_csms.to_string()),
            };
            rows.push((style.low_color.clone(), format!("{}: {lo}", gradient.label())));
            rows.push((style.high_color.clone(), format!("{}: {hi}", gradient.label())));
        }
    }
    let csms: u32 = all.iter().map(|l| l.csms).sum();
    (rows, format!("{} crosslinks, {csms} CSMs, q ≤ {}", all.len(), style.max_q))
}

/// A character's average width in a style, in ems (to lay text out)
fn char_width(t: &TypeStyle) -> f64 {
    if t.bold {
        0.62
    } else {
        0.56
    }
}

/// The circle's centre and ring radius: the radius as large as it may be
/// (`r1` at most, `r_min` at least) and the centre moved as little as it
/// must, so that the circle and all around it (`outer` beyond the ring:
/// ticks, numbers and labels) keep `gap` clear of the legend's box
/// `[x0, y0, x1, y1]`, within the square canvas (`side`, below `top`)
#[allow(clippy::too_many_arguments)]
fn fit_circle(side: f64, top: f64, centre: (f64, f64), r1: f64, r_min: f64, outer: f64, legend: Option<[f64; 4]>, gap: f64) -> (f64, f64, f64) {
    let Some([x0, y0, x1, y1]) = legend else { return (centre.0, centre.1, r1) };
    let clear = |(cx, cy): (f64, f64), r: f64| {
        let dx = (x0 - cx).max(cx - x1).max(0.0);
        let dy = (y0 - cy).max(cy - y1).max(0.0);
        dx.hypot(dy) >= r + outer + gap
    };
    // away from the legend, by a share `t` of the room on each side
    let right = (x0 + x1) / 2.0 < centre.0;
    let down = (y0 + y1) / 2.0 < centre.1;
    let place = |r: f64, t: f64| {
        let reach = r + outer;
        let room_x = if right { side - (centre.0 + reach) } else { centre.0 - reach }.max(0.0);
        let room_y = if down { side - (centre.1 + reach) } else { centre.1 - reach - top }.max(0.0);
        let (sx, sy) = (if right { 1.0 } else { -1.0 }, if down { 1.0 } else { -1.0 });
        (centre.0 + sx * room_x * t, centre.1 + sy * room_y * t)
    };
    let fits = |r: f64| (0..=40).map(|i| place(r, i as f64 / 40.0)).find(|&c| clear(c, r));
    if let Some((cx, cy)) = fits(r1) {
        return (cx, cy, r1);
    }
    if fits(r_min).is_none() {
        let (cx, cy) = place(r_min, 1.0);
        return (cx, cy, r_min);
    }
    // the largest radius that fits
    let (mut lo, mut hi) = (r_min, r1);
    for _ in 0..30 {
        let mid = (lo + hi) / 2.0;
        if fits(mid).is_some() {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let (cx, cy) = fits(lo).unwrap_or(centre);
    (cx, cy, lo)
}

// ---------------------------------------------------------------- Circos files

/// A Circos id: letters, digits and `_` only
fn circos_id(i: usize) -> String {
    format!("p{}", i + 1)
}

fn circos_label(label: &str) -> String {
    let s: String = label.chars().map(|c| if c.is_ascii_alphanumeric() || c == '_' || c == '-' || c == '.' { c } else { '_' }).collect();
    if s.is_empty() { "protein".into() } else { s }
}

fn circos_rgb(color: &str) -> String {
    let (r, g, b) = rgb(color);
    format!("{r},{g},{b}")
}

/// The plot as Circos input files in `dir`: `karyotype.txt` (the proteins),
/// `links.txt`, the tracks, and a `circos.conf` with the style's colours and
/// layout (run `circos -conf circos.conf` there). Returns the files written.
pub fn write_circos_files(data: &CircosData, style: &CircosStyle, dir: &Path) -> Result<Vec<PathBuf>, String> {
    std::fs::create_dir_all(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
    let proteins = shown_proteins(data, style);
    let ids: HashMap<&str, String> = proteins.iter().enumerate().map(|(i, (p, _, _))| (p.accession.as_str(), circos_id(i))).collect();
    let mut written = Vec::new();
    let mut put = |name: &str, text: String| -> Result<(), String> {
        let path = dir.join(name);
        std::fs::write(&path, text).map_err(|e| format!("{}: {e}", path.display()))?;
        written.push(path);
        Ok(())
    };

    let mut karyotype = String::from("# XNova crosslink plot: one Circos chromosome per protein (residues)\n");
    for (i, (p, label, _)) in proteins.iter().enumerate() {
        let _ = writeln!(karyotype, "chr - {} {} 0 {} {}", circos_id(i), circos_label(label), p.sequence.len(), circos_id(i));
    }
    put("karyotype.txt", karyotype)?;

    let mut links = String::from("# protein residue protein residue options (one crosslink per line)\n");
    let plotted = |l: &&Link| ids.contains_key(l.a.0.as_str()) && ids.contains_key(l.b.0.as_str());
    let dim = style.dim_others && data.links.iter().filter(plotted).any(|l| highlight_of(style, l).is_some());
    for l in data.links.iter().filter(plotted) {
        let (Some(a), Some(b)) = (ids.get(l.a.0.as_str()), ids.get(l.b.0.as_str())) else { continue };
        let highlight = highlight_of(style, l);
        let color = match (highlight, l.class) {
            (Some(i), _) => format!("hl{}", i + 1),
            (None, _) if dim => "xl_dim".into(),
            (None, Class::InterProtein) => "xl_inter".into(),
            (None, Class::IntraProtein) => "xl_intra".into(),
            (None, Class::Loop) => "xl_loop".into(),
            (None, _) => "xl_multi".into(),
        };
        let mut thickness = if style.width_by_csms { 1 + ((l.csms as f64).log2() * 1.5).round() as u32 } else { 2 };
        // highlighted: wider, and drawn on top
        let mut z = "";
        if highlight.is_some() {
            thickness = ((thickness as f64) * style.highlight_width.clamp(0.2, 10.0)).round().max(1.0) as u32;
            z = ",z=10";
        }
        let _ = writeln!(
            links,
            "{a} {} {} {b} {} {} color={color},thickness={thickness}{z},csms={},q={},score={:.2}",
            l.a.1, l.a.1, l.b.1, l.b.1, l.csms, l.q, l.score
        );
    }
    put("links.txt", links)?;

    let mut sites: BTreeMap<(String, usize), u32> = BTreeMap::new();
    for l in &data.links {
        for (acc, pos) in [&l.a, &l.b] {
            if let Some(id) = ids.get(acc.as_str()) {
                *sites.entry((id.clone(), *pos)).or_default() += l.csms;
            }
        }
    }
    let mut text = String::from("# crosslinked residues: protein start end CSMs\n");
    for ((id, pos), csms) in &sites {
        let _ = writeln!(text, "{id} {} {pos} {csms}", pos - 1);
    }
    put("sites.txt", text)?;

    let mut text = String::from("# residues the linker reacts with\n");
    for (p, _, _) in &proteins {
        for (i, residue) in p.sequence.chars().enumerate() {
            if data.reactive.contains(residue) {
                let _ = writeln!(text, "{} {} {}", ids[p.accession.as_str()], i, i + 1);
            }
        }
    }
    put("reactive.txt", text)?;

    let mut text = String::from("# sequence covered by the linear search\n");
    for (p, _, _) in &proteins {
        for &(s0, s1) in data.coverage.get(&p.accession).map(Vec::as_slice).unwrap_or(&[]) {
            let _ = writeln!(text, "{} {} {}", ids[p.accession.as_str()], s0 - 1, s1);
        }
    }
    put("coverage.txt", text)?;

    let mut regions = String::from("# annotated regions\n");
    let mut region_labels = String::from("# annotated regions' labels\n");
    let mut region_colors = String::new();
    for (i, g) in data.regions.iter().enumerate() {
        let Some(id) = ids.get(g.accession.as_str()) else { continue };
        let color = if g.color.is_empty() { Palette::OkabeIto.color(i).to_string() } else { g.color.clone() };
        let _ = writeln!(region_colors, "region{} = {}", i + 1, circos_rgb(&color));
        let _ = writeln!(regions, "{id} {} {} fill_color=region{}", g.start - 1, g.end, i + 1);
        if !g.label.is_empty() {
            let _ = writeln!(region_labels, "{id} {} {} {}", g.start - 1, g.end, circos_label(&g.label));
        }
    }
    put("regions.txt", regions)?;
    put("region_labels.txt", region_labels)?;

    // colours: each protein's, the links', the regions'
    let mut colors = String::new();
    for (i, (_, _, color)) in proteins.iter().enumerate() {
        let _ = writeln!(colors, "{} = {}", circos_id(i), circos_rgb(color));
    }
    for (name, color) in [
        ("xl_inter", &style.inter_color),
        ("xl_intra", &style.intra_color),
        ("xl_loop", &style.loop_color),
        ("xl_multi", &style.multi_color),
    ] {
        let color = match style.color_by {
            ColorBy::Single => &style.single_color,
            _ => color,
        };
        let _ = writeln!(colors, "{name} = {}", circos_rgb(color));
    }
    let _ = writeln!(colors, "xl_dim = {}", circos_rgb(&style.dim_color));
    for (i, h) in style.highlights.iter().enumerate() {
        let _ = writeln!(colors, "hl{} = {}", i + 1, circos_rgb(&h.color));
    }
    colors.push_str(&region_colors);
    let longest = proteins.iter().map(|(p, _, _)| p.sequence.len()).max().unwrap_or(1);
    let (major, minor) = tick_steps(longest, style.tick_spacing);
    let bezier = style.bezier_radius.clamp(0.0, 1.0);
    let radius_px = (style.size.clamp(2.0, 40.0) * 150.0 * style.radius.clamp(0.2, 0.95)).round();
    let yes = |b: bool| if b { "yes" } else { "no" };
    // text sizes in the image's pixels: the image is twice `radius_px` wide
    let px_per_pt = 2.0 * radius_px / style.size.clamp(2.0, 40.0) / 72.0;
    let circos_px = |t: &TypeStyle| (t.size.clamp(2.0, 96.0) * px_per_pt).round().max(4.0);
    let mut conf = format!(
        r#"# Circos configuration written by XNova (see README.txt). Run:  circos -conf circos.conf
karyotype = karyotype.txt
chromosomes_units = 1
chromosomes_display_default = yes

<colors>
<<include etc/colors.conf>>
<<include etc/brewer.conf>>
{colors}</colors>
# text in Circos' own fonts (etc/fonts.conf), bold and italic as in XNova
<fonts>
<<include etc/fonts.conf>>
</fonts>
<patterns>
<<include etc/patterns.conf>>
</patterns>

<ideogram>
<spacing>
default = {gap}u
</spacing>
radius = 0.82r
thickness = 25p
fill = yes
stroke_color = black
stroke_thickness = 1p
show_label = {labels}
label_font = {label_font}
label_radius = dims(ideogram,radius_outer) + {label_offset}p
label_size = {label_px}p
label_parallel = {parallel}
</ideogram>

show_ticks = {ticks}
show_tick_labels = {tick_labels}
<ticks>
radius = dims(ideogram,radius_outer)
color = black
thickness = 2p
multiplier = 1
format = %d
<tick>
spacing = {minor}u
size = 6p
</tick>
<tick>
spacing = {major}u
size = 12p
show_label = yes
label_size = {tick_px}p
label_font = {tick_font}
label_offset = 4p
</tick>
</ticks>

<plots>
<plot>
type = highlight
file = coverage.txt
r0 = 0.94r
r1 = 0.97r
fill_color = grey
show = {coverage}
</plot>
<plot>
type = highlight
file = regions.txt
r0 = 0.89r
r1 = 0.93r
</plot>
<plot>
type = text
file = region_labels.txt
r0 = 0.89r
r1 = 0.99r
label_size = {region_px}p
label_font = {region_font}
</plot>
<plot>
type = highlight
file = reactive.txt
r0 = 0.86r
r1 = 0.88r
fill_color = dgrey
show = {reactive}
</plot>
<plot>
type = histogram
file = sites.txt
r0 = 0.80r
r1 = 0.85r
fill_color = vdgrey
color = black
thickness = 1p
show = {sites}
</plot>
</plots>

<links>
<link>
file = links.txt
radius = 0.79r
bezier_radius = {bezier}r
crest = 0.4
ribbon = no
</link>
</links>

<image>
<<include etc/image.conf>>
radius* = {radius_px}p
background* = {background}
</image>

<<include etc/housekeeping.conf>>
"#,
        gap = (style.gap * longest as f64 / 60.0).round().max(1.0),
        label_font = style.label_text.circos_font(),
        label_px = circos_px(&style.label_text),
        tick_font = style.tick_text.circos_font(),
        tick_px = circos_px(&style.tick_text),
        region_font = style.region_text.circos_font(),
        region_px = circos_px(&style.region_text),
        labels = yes(style.label_mode != LabelMode::None),
        parallel = yes(style.label_mode == LabelMode::Parallel),
        label_offset = if style.tick_labels { 70 } else { 20 },
        ticks = yes(style.ticks),
        tick_labels = yes(style.tick_labels),
        coverage = yes(style.coverage),
        reactive = yes(style.reactive),
        sites = yes(style.sites),
        background = match style.background {
            Background::Black => "black",
            Background::Transparent => "transparent",
            Background::White => "white",
        },
    );
    if data.regions.is_empty() {
        conf = conf.replace("file = regions.txt\n", "file = regions.txt\nshow = no\n").replace(
            "file = region_labels.txt\n",
            "file = region_labels.txt\nshow = no\n",
        );
    }
    put("circos.conf", conf)?;
    put(
        "README.txt",
        "Circos input files written by XNova from its crosslink results.\n\n\
         Circos (http://circos.ca, Martin Krzywinski; GPL) draws them: install it (it needs Perl and the GD library),\n\
         then run, in this folder:\n\n    circos -conf circos.conf\n\n\
         It writes circos.png and circos.svg. Every file is plain text, so the configuration can be edited freely:\n\
         karyotype.txt (one 'chromosome' per protein, in residues), links.txt (one crosslink per line, with its CSMs,\n\
         q-value and score), sites.txt (CSMs per crosslinked residue), reactive.txt (the linker's reactive residues),\n\
         coverage.txt (the linear search's sequence coverage), regions.txt and region_labels.txt (annotations).\n"
            .to_string(),
    )?;
    Ok(written)
}

/// The number of crosslinks between each pair of proteins, as a table for
/// Circos' table viewer (rows and columns are proteins; each pair once)
pub fn crosslink_table(data: &CircosData, style: &CircosStyle) -> String {
    let proteins = shown_proteins(data, style);
    let index: HashMap<&str, usize> = proteins.iter().enumerate().map(|(i, (p, _, _))| (p.accession.as_str(), i)).collect();
    let n = proteins.len();
    let mut counts = vec![vec![0u32; n]; n];
    for l in &data.links {
        let (Some(&i), Some(&j)) = (index.get(l.a.0.as_str()), index.get(l.b.0.as_str())) else { continue };
        let (i, j) = (i.min(j), i.max(j));
        counts[i][j] += 1;
    }
    let labels: Vec<String> = proteins.iter().map(|(_, l, _)| circos_label(l)).collect();
    let mut out = format!("data\t{}\n", labels.join("\t"));
    for (i, row) in counts.iter().enumerate() {
        let cells: Vec<String> = row.iter().map(|&c| if c == 0 { "-".into() } else { c.to_string() }).collect();
        let _ = writeln!(out, "{}\t{}", labels[i], cells.join("\t"));
    }
    out
}

/// Write a search's plot with `style` into the results folder's
/// [`FOLDER`]: `crosslinks_circos.svg`, `.pdf` and `.png` (300 dpi), the
/// style (`circos_style.json`, to load and change in sage-gui), the Circos
/// input files (`circos_input/`) and the table-viewer table. A search
/// without crosslinks gets none. Returns the files written.
pub fn write_figures(folder: &Path, style: &CircosStyle) -> Result<Vec<PathBuf>, String> {
    let results = Results::load(folder)?;
    if !results.matches.iter().any(|m| m.kind != Kind::Linear) {
        return Ok(Vec::new());
    }
    let data = build(&results, style);
    if data.links.is_empty() {
        return Ok(Vec::new());
    }
    let dir = folder.join(FOLDER);
    std::fs::create_dir_all(&dir).map_err(|e| format!("{}: {e}", dir.display()))?;
    let svg = svg(&data, style);
    let mut written = Vec::new();
    for format in [figure::Format::Svg, figure::Format::Pdf, figure::Format::Png] {
        let path = dir.join(format!("crosslinks_circos.{}", format.extension()));
        figure::export(&svg, format, 300, style.background != Background::Transparent, &path)?;
        written.push(path);
    }
    let path = dir.join("circos_style.json");
    std::fs::write(&path, style.to_json()).map_err(|e| format!("{}: {e}", path.display()))?;
    written.push(path);
    let path = dir.join("crosslink_table.txt");
    std::fs::write(&path, crosslink_table(&data, style)).map_err(|e| format!("{}: {e}", path.display()))?;
    written.push(path);
    written.extend(write_circos_files(&data, style, &dir.join("circos_input"))?);
    Ok(written)
}

#[cfg(test)]
mod test {
    use super::*;

    fn segment(accession: &str, label: &str, length: usize) -> Segment {
        Segment { accession: accession.into(), label: label.into(), sequence: "K".repeat(length), links: 1 }
    }

    fn link(a: (&str, usize), b: (&str, usize), class: Class, csms: u32) -> Link {
        Link { a: (a.0.into(), a.1), b: (b.0.into(), b.1), class, csms, q: 0.001, score: 20.0 + csms as f64 }
    }

    fn data() -> CircosData {
        CircosData {
            proteins: vec![segment("sp|P1|ONE", "ONE", 400), segment("sp|P2|TWO", "TWO", 200), segment("sp|P3|THREE", "THREE", 100)],
            links: vec![
                link(("sp|P1|ONE", 10), ("sp|P2|TWO", 50), Class::InterProtein, 3),
                link(("sp|P1|ONE", 100), ("sp|P1|ONE", 300), Class::IntraProtein, 1),
                link(("sp|P2|TWO", 20), ("sp|P2|TWO", 25), Class::Loop, 2),
            ],
            coverage: HashMap::from([("sp|P1|ONE".to_string(), vec![(5, 60)])]),
            regions: vec![Region { accession: "sp|P2|TWO".into(), start: 1, end: 80, label: "domain".into(), color: String::new() }],
            reactive: "K".into(),
            notes: Vec::new(),
        }
    }

    /// The SVG holds every protein, link and track the style asks for
    #[test]
    fn draws_proteins_links_and_tracks() {
        let style = CircosStyle { reactive: true, title: "Pol II".into(), ..CircosStyle::default() };
        let drawn = svg(&data(), &style);
        assert!(drawn.starts_with("<svg") && drawn.ends_with("</svg>"));
        assert!(drawn.contains(r#"width="8in""#));
        for label in ["ONE", "TWO", "THREE", "Pol II", "domain"] {
            assert!(drawn.contains(&format!(">{label}<")), "{label}");
        }
        assert_eq!(drawn.matches(" Q").count(), 3, "three links");
        assert!(drawn.contains(&style.inter_color) && drawn.contains(&style.loop_color));
        assert!(drawn.contains("3 crosslinks, 6 CSMs"), "legend");
        // the SVG renders
        let image = figure::rasterize(&drawn, 30).unwrap();
        assert_eq!((image.width(), image.height()), (240, 240));

        // a hidden protein takes its links along; others keep their colours
        let hidden = CircosStyle { hidden: vec!["sp|P3|THREE".into()], ..CircosStyle::default() };
        let plotted = data();
        let shown = shown_proteins(&plotted, &hidden);
        assert_eq!(shown.len(), 2);
        assert_eq!(shown[1].2, Palette::Tableau20.color(1));
        assert!(!svg(&data(), &hidden).contains(">THREE<"));
        // order and labels from the style
        let ordered = CircosStyle {
            order: vec!["sp|P3|THREE".into()],
            labels: BTreeMap::from([("sp|P2|TWO".to_string(), "Rpb2".to_string())]),
            ..CircosStyle::default()
        };
        let shown = shown_proteins(&plotted, &ordered);
        assert_eq!(shown[0].0.accession, "sp|P3|THREE");
        assert_eq!(shown[2].1, "Rpb2");
        // nothing to draw
        let empty = svg(&CircosData::default(), &CircosStyle::default());
        assert!(empty.contains("No crosslinks"));
    }

    /// Link colours: by class, single, gradients
    #[test]
    fn link_colours() {
        let single = CircosStyle { color_by: ColorBy::Single, single_color: "#123456".into(), ..CircosStyle::default() };
        assert_eq!(svg(&data(), &single).matches("#123456").count(), 3);
        let by_score = CircosStyle { color_by: ColorBy::Score, ..CircosStyle::default() };
        let drawn = svg(&data(), &by_score);
        assert!(drawn.contains(&by_score.high_color), "the best link at the high end");
        assert!(drawn.contains(&by_score.low_color), "the worst at the low end");
        assert_eq!(mix("#000000", "#ffffff", 0.5), "#808080");
        assert_eq!(rgb("#ff8000"), (255, 128, 0));
        assert_eq!(rgb("nonsense"), (128, 128, 128));
    }

    /// Styles round-trip through JSON; missing keys take their defaults
    #[test]
    fn style_files() {
        let style = CircosStyle {
            max_q: 0.05,
            labels: BTreeMap::from([("sp|P1|ONE".to_string(), "Rpb1".to_string())]),
            color_by: ColorBy::QValue,
            background: Background::Black,
            ..CircosStyle::default()
        };
        assert_eq!(CircosStyle::from_json(&style.to_json()).unwrap(), style);
        let partial = CircosStyle::from_json(r#"{"gap": 5, "palette": "okabe_ito"}"#).unwrap();
        assert_eq!((partial.gap, partial.palette, partial.max_q), (5.0, Palette::OkabeIto, 0.01));
        assert!(CircosStyle::from_json(r#"{"palette": "neon"}"#).is_err());
        // presets keep the selection
        let dark = Preset::Dark.apply(&style);
        assert_eq!((dark.background, dark.max_q), (Background::Black, 0.05));
        assert_eq!(dark.labels, style.labels);
    }

    /// Templates: the style alone, or with the proteins as they look now
    #[test]
    fn templates() {
        let plotted = data();
        let styled = CircosStyle {
            size: 6.5,
            font: "Helvetica".into(),
            color_by: ColorBy::Score,
            palette: Palette::OkabeIto,
            title: "Pol II".into(),
            order: vec!["sp|P2|TWO".into()],
            hidden: vec!["sp|P3|THREE".into()],
            labels: BTreeMap::from([("sp|P1|ONE".to_string(), "Rpb1".to_string())]),
            colors: BTreeMap::from([("sp|P2|TWO".to_string(), "#123456".to_string())]),
            ..CircosStyle::default()
        };
        // the style alone: other results keep their own proteins
        let look = Template::new(" Paper ", &styled, false, &plotted);
        assert_eq!(look.name, "Paper");
        assert!(look.style.order.is_empty() && look.style.colors.is_empty() && look.style.title.is_empty());
        let other = CircosStyle { title: "Mediator".into(), labels: BTreeMap::from([("sp|M|MED".to_string(), "Med1".to_string())]), ..CircosStyle::default() };
        let applied = look.apply(&other);
        assert_eq!((applied.size, applied.color_by, applied.palette), (6.5, ColorBy::Score, Palette::OkabeIto));
        assert_eq!((applied.title.as_str(), &applied.labels), ("Mediator", &other.labels));
        assert!(look.summary().contains("6.5 in, Helvetica 11 pt"), "{}", look.summary());

        // with the proteins: order and every colour as now, the hidden one last
        let same = Template::new("Pol II", &styled, true, &plotted);
        assert_eq!(same.style.order, ["sp|P2|TWO", "sp|P1|ONE", "sp|P3|THREE"]);
        assert_eq!(same.style.colors["sp|P2|TWO"], "#123456");
        assert_eq!(same.style.colors["sp|P1|ONE"], Palette::OkabeIto.color(0));
        assert_eq!(same.style.colors["sp|P3|THREE"], Palette::OkabeIto.color(2));
        let applied = same.apply(&other);
        assert_eq!((applied.title.as_str(), &applied.labels, &applied.hidden), ("Pol II", &styled.labels, &styled.hidden));
        // drawn the same from other results that list the proteins differently
        let mut reordered = data();
        reordered.proteins.reverse();
        let a: Vec<_> = shown_proteins(&plotted, &applied).into_iter().map(|(p, l, c)| (p.accession.clone(), l, c)).collect();
        let b: Vec<_> = shown_proteins(&reordered, &applied).into_iter().map(|(p, l, c)| (p.accession.clone(), l, c)).collect();
        assert_eq!(a, b);
        assert!(same.summary().contains("3 proteins"));
        // through JSON, as the library keeps them
        let text = serde_json::to_string(&same).unwrap();
        assert_eq!(serde_json::from_str::<Template>(&text).unwrap(), same);
        let old: Template = serde_json::from_str(r#"{"name": "x", "style": {"gap": 4}}"#).unwrap();
        assert!(!old.proteins && old.style.gap == 4.0);
    }

    /// Each kind of text in its own size, bold and italic, in the font chosen
    #[test]
    fn text_styles() {
        let mut style = CircosStyle { title: "Pol II".into(), font: font_list("Times New Roman"), ..CircosStyle::default() };
        style.label_text = TypeStyle { size: 14.0, bold: false, italic: true };
        style.tick_text.bold = true;
        style.region_text.italic = true;
        let drawn = svg(&data(), &style);
        let text_of = |label: &str| drawn.split("<text").find(|t| t.contains(&format!(">{label}<"))).unwrap_or_else(|| panic!("{label}")).to_string();
        let one = text_of("ONE");
        assert!(one.contains(&format!(r#"font-size="{:.1}""#, pt(14.0))), "{one}");
        assert!(one.contains(r#"font-style="italic""#) && !one.contains("bold"), "{one}");
        assert!(one.contains(r#"font-family="Times New Roman, serif""#), "{one}");
        assert!(text_of("Pol II").contains(r#"font-weight="bold""#));
        assert!(text_of("100").contains(r#"font-weight="bold""#), "tick numbers");
        assert!(text_of("domain").contains("italic"));
        assert!(!text_of("domain").contains("bold"));
        assert!(figure::rasterize(&drawn, 30).is_ok());

        style.scale_text(1.1);
        assert_eq!((style.label_text.size, style.tick_text.size), (15.4, 8.8));
        // presets keep the text
        assert_eq!(Preset::Dark.apply(&style).label_text, style.label_text);
        assert_eq!(Preset::Dark.apply(&style).font, style.font);

        assert_eq!(font_list("Courier New"), "Courier New, monospace");
        assert_eq!(font_list("Open Sans"), "Open Sans, sans-serif");
        assert_eq!(font_list("Georgia"), "Georgia, serif");
        assert_eq!(font_list("serif"), "serif");
        assert_eq!(font_list("Noto Sans CJK 2"), "'Noto Sans CJK 2', sans-serif");
        assert_eq!(first_font("'Noto Sans CJK 2', sans-serif"), "Noto Sans CJK 2");

        // the Circos files: the same styles in Circos' fonts
        let dir = std::env::temp_dir().join(format!("xnova_circos_text_{}", std::process::id()));
        write_circos_files(&data(), &style, &dir).unwrap();
        let conf = std::fs::read_to_string(dir.join("circos.conf")).unwrap();
        for needed in ["label_font = italic", "label_font = bold\nlabel_offset", "label_font = italic\n</plot>"] {
            assert!(conf.contains(needed), "{needed}");
        }
        let _ = std::fs::remove_dir_all(&dir);
        // older styles without the text keys take the defaults
        let old = CircosStyle::from_json(r#"{"font": "Arial", "font_size": 9}"#).unwrap();
        assert_eq!(old.label_text, CircosStyle::default().label_text);
    }

    /// The circle keeps clear of the legend: moved, then smaller if need be
    #[test]
    fn legend_keeps_clear() {
        let side = 1200.0;
        let clear = |(cx, cy, r): (f64, f64, f64), b: [f64; 4]| {
            let dx = (b[0] - cx).max(cx - b[2]).max(0.0);
            let dy = (b[1] - cy).max(cy - b[3]).max(0.0);
            dx.hypot(dy) >= r + 60.0 + 5.0 - 1e-6
        };
        let fit = |legend: Option<[f64; 4]>| fit_circle(side, 0.0, (600.0, 600.0), 450.0, 120.0, 60.0, legend, 5.0);
        // no legend, or one out of the way: as it was
        assert_eq!(fit(None), (600.0, 600.0, 450.0));
        assert_eq!(fit(Some([36.0, 1130.0, 100.0, 1164.0])), (600.0, 600.0, 450.0));
        // in the way: the circle moves up and to the right, as large as before
        let legend = [36.0, 1000.0, 300.0, 1164.0];
        let moved = fit(Some(legend));
        assert!(clear(moved, legend), "{moved:?}");
        assert!(moved.0 > 600.0 && moved.1 < 600.0 && moved.2 == 450.0, "{moved:?}");
        // a top-right legend: down and to the left
        let legend = [900.0, 36.0, 1164.0, 200.0];
        let moved = fit(Some(legend));
        assert!(clear(moved, legend) && moved.0 < 600.0 && moved.1 > 600.0, "{moved:?}");
        // too large for that: smaller as well
        let legend = [36.0, 700.0, 700.0, 1164.0];
        let smaller = fit(Some(legend));
        assert!(clear(smaller, legend) && smaller.2 < 450.0 && smaller.2 > 120.0, "{smaller:?}");

        // in the SVG: a large legend at the top right, everything drawn
        let style = CircosStyle { legend_text: TypeStyle::new(14.0, false), legend_corner: LegendCorner::TopRight, ..CircosStyle::default() };
        let drawn = svg(&data(), &style);
        let summary = drawn.split("<text").find(|t| t.contains("3 crosslinks, 6 CSMs")).unwrap();
        let x: f64 = summary.split("x=\"").nth(1).unwrap().split('"').next().unwrap().parse().unwrap();
        let y: f64 = summary.split("y=\"").nth(1).unwrap().split('"').next().unwrap().parse().unwrap();
        assert!(x > 600.0 && y < 300.0, "({x}, {y})");
        assert!(figure::rasterize(&drawn, 30).is_ok());
        assert_eq!(CircosStyle::from_json(&style.to_json()).unwrap().legend_corner, LegendCorner::TopRight);
    }

    /// Highlights: crosslinks picked out by residue, protein, pair, kind or
    /// stretch, in their colour, wider and on top; the others dimmed
    #[test]
    fn highlights() {
        let plotted = data();
        let (inter, intra, looped) = (&plotted.links[0], &plotted.links[1], &plotted.links[2]);
        let one = |p: &Pick| plotted.links.iter().filter(|l| p.picks(l)).count();
        let link = Pick::Link { a: ("sp|P2|TWO".into(), 50), b: ("sp|P1|ONE".into(), 10) };
        assert!(link.picks(inter) && !link.picks(intra));
        assert_eq!(one(&Pick::Protein { protein: "sp|P1|ONE".into() }), 2);
        assert_eq!(one(&Pick::Between { a: "sp|P2|TWO".into(), b: "sp|P1|ONE".into() }), 1);
        assert_eq!(one(&Pick::Between { a: "sp|P2|TWO".into(), b: "sp|P2|TWO".into() }), 1);
        assert!(Pick::Kind { kind: LinkKind::Loop }.picks(looped));
        assert_eq!(one(&Pick::Region { protein: "sp|P1|ONE".into(), start: 90, end: 110 }), 1);
        assert_eq!(link.describe(&plotted, &CircosStyle::default()), "TWO K50 – ONE K10");
        assert_eq!(Pick::Between { a: "sp|P2|TWO".into(), b: "sp|P2|TWO".into() }.describe(&plotted, &CircosStyle::default()), "within TWO");

        let mut style = CircosStyle::default();
        style.highlights = vec![
            Highlight::new(Pick::Protein { protein: "sp|P1|ONE".into() }, "#0072b2"),
            Highlight { label: "loops".into(), ..Highlight::new(Pick::Kind { kind: LinkKind::Loop }, "#009e73") },
            // picks nothing new: the first highlight has both of ONE's links
            Highlight::new(Pick::Region { protein: "sp|P1|ONE".into(), start: 1, end: 400 }, "#d55e00"),
        ];
        assert_eq!(highlight_of(&style, inter), Some(0));
        assert_eq!(highlight_of(&style, looped), Some(1));
        let drawn = svg(&plotted, &style);
        let paths: Vec<&str> = drawn.split("<path d=\"M").filter(|p| p.contains(" Q")).collect();
        assert_eq!(paths.len(), 3);
        assert!(paths.iter().all(|p| p.contains("#0072b2") || p.contains("#009e73")), "highlight colours");
        assert!(paths.last().unwrap().contains("#0072b2"), "highlighted with the most CSMs on top");
        assert!(drawn.contains("crosslinks of ONE (2)") && drawn.contains("loops (1)"), "legend");
        assert!(!drawn.contains("#d55e00") && !drawn.contains("between proteins"), "no empty rows");

        // the others dimmed
        style.highlights.truncate(1);
        style.dim_others = true;
        let drawn = svg(&plotted, &style);
        assert!(drawn.contains(&format!(r#"stroke="{}""#, style.dim_color)) && drawn.contains("other crosslinks (1)"));
        // unticked: drawn as before, and nothing dimmed
        style.highlights[0].on = false;
        let drawn = svg(&plotted, &style);
        assert!(!drawn.contains(&style.dim_color) && drawn.contains("between proteins (1)"));
        style.highlights[0].on = true;

        // through JSON
        assert_eq!(CircosStyle::from_json(&style.to_json()).unwrap(), style);
        assert!(style.to_json().contains(r#""pick": "protein""#));
        // in the Circos files: their colour, wider, on top
        let dir = std::env::temp_dir().join(format!("xnova_circos_highlights_{}", std::process::id()));
        write_circos_files(&plotted, &style, &dir).unwrap();
        let links = std::fs::read_to_string(dir.join("links.txt")).unwrap();
        let conf = std::fs::read_to_string(dir.join("circos.conf")).unwrap();
        assert!(links.contains("color=hl1,thickness=5,z=10") && links.contains("color=xl_dim"), "{links}");
        assert!(conf.contains("hl1 = 0,114,178") && conf.contains("xl_dim = 200,200,200"));
        let _ = std::fs::remove_dir_all(&dir);

        // templates: highlights of a kind go with the style, those of proteins with the proteins
        let kinds = CircosStyle { highlights: vec![Highlight::new(Pick::Kind { kind: LinkKind::Multivalent }, "#000000")], ..CircosStyle::default() };
        let look = Template::new("look", &kinds, false, &plotted);
        let applied = look.apply(&style);
        let picks: Vec<_> = applied.highlights.iter().map(|h| h.pick.clone()).collect();
        assert_eq!(picks, [Pick::Protein { protein: "sp|P1|ONE".into() }, Pick::Kind { kind: LinkKind::Multivalent }]);
        assert!(Template::new("x", &style, false, &plotted).style.highlights.is_empty());
        assert_eq!(Template::new("x", &style, true, &plotted).style.highlights, style.highlights);
    }

    /// Circos input files and the table-viewer table
    #[test]
    fn circos_files_and_table() {
        let dir = std::env::temp_dir().join(format!("xnova_circos_{}", std::process::id()));
        let written = write_circos_files(&data(), &CircosStyle::default(), &dir).unwrap();
        let read = |name: &str| std::fs::read_to_string(dir.join(name)).unwrap();
        assert!(written.iter().any(|p| p.ends_with("circos.conf")));
        assert_eq!(read("karyotype.txt").lines().filter(|l| l.starts_with("chr - ")).count(), 3);
        assert!(read("karyotype.txt").contains("chr - p1 ONE 0 400 p1"));
        let links = read("links.txt");
        assert!(links.contains("p1 10 10 p2 50 50 color=xl_inter,thickness="), "{links}");
        assert_eq!(links.lines().filter(|l| !l.starts_with('#')).count(), 3);
        let conf = read("circos.conf");
        for needed in ["karyotype = karyotype.txt", "<<include etc/housekeeping.conf>>", "file = links.txt", "xl_inter = 255,82,59", "p1 = 78,121,167"] {
            assert!(conf.contains(needed), "{needed}");
        }
        assert!(read("sites.txt").contains("p1 9 10 3"));
        let _ = std::fs::remove_dir_all(&dir);

        let table = crosslink_table(&data(), &CircosStyle::default());
        let mut lines = table.lines();
        assert_eq!(lines.next().unwrap(), "data\tONE\tTWO\tTHREE");
        assert_eq!(lines.next().unwrap(), "ONE\t1\t1\t-");
        assert_eq!(lines.next().unwrap(), "TWO\t-\t1\t-");
    }

    /// Tick spacing: chosen from the longest protein, or the style's
    #[test]
    fn ticks() {
        assert_eq!(tick_steps(1733, 0), (500, 100));
        assert_eq!(tick_steps(200, 0), (50, 10));
        assert_eq!(tick_steps(200, 25), (25, 5));
    }
}
