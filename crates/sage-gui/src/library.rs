//! Crosslinkers and modifications saved under a name for later searches: the
//! "Saved crosslinkers" panel of the Crosslinking tab and the "Saved
//! modifications" panel of the Modifications tab. They are kept on this
//! computer in `%APPDATA%\sage-gui\library.json`, apart from the search
//! settings, and written as soon as one is saved, edited or deleted.

use crate::settings::{self, Crosslink, ModKind, Modification, Reagent, Site};
use crate::theme;
use eframe::egui::{self, Color32, ComboBox, RichText, TextEdit, Ui};
use serde_json::{json, Map, Value};
use std::path::{Path, PathBuf};

/// A crosslinker's chemistry, saved under a name
#[derive(Clone, Debug, PartialEq)]
pub struct SavedLinker {
    pub name: String,
    /// Mass the linker adds to the pair of peptides (Da)
    pub mass: f64,
    /// MS-cleavable, leaving the stubs below
    pub cleavable: bool,
    pub stub_short: f64,
    pub stub_long: f64,
    /// Extra stub masses, comma separated
    pub stub_alt: String,
    /// Cleavable only: the second arm leaves its own stubs
    pub asymmetric: bool,
    pub stub2_short: f64,
    pub stub2_long: f64,
    pub stub2_alt: String,
    /// Each end reacts with a different group (`site1`, `site2`), else both
    /// with `residues` (and the protein N-terminus)
    pub two_sites: bool,
    pub residues: String,
    pub protein_nterm: bool,
    pub site1: Site,
    pub site2: Site,
}

impl SavedLinker {
    /// The crosslinker set up in `xl` (a preset or a custom one)
    pub fn from_crosslink(name: &str, xl: &Crosslink) -> Self {
        let cleavable = xl.is_cleavable();
        Self {
            name: name.trim().to_string(),
            mass: xl.linker_mass(),
            cleavable,
            stub_short: xl.stub_short,
            stub_long: xl.stub_long,
            stub_alt: xl.stub_alt.trim().to_string(),
            asymmetric: cleavable && xl.has_second_arm(),
            stub2_short: xl.stub2_short,
            stub2_long: xl.stub2_long,
            stub2_alt: xl.stub2_alt.trim().to_string(),
            two_sites: xl.two_sites,
            residues: xl.residues.trim().to_string(),
            protein_nterm: xl.protein_nterm,
            site1: xl.site1.clone(),
            site2: xl.site2.clone(),
        }
    }

    /// Set up `xl` with this crosslinker (as a custom reagent); the search
    /// settings (backend, peaks, loop-links, ...) are left as they are
    pub fn apply(&self, xl: &mut Crosslink) {
        xl.reagent = Reagent::Custom;
        xl.custom_mass = self.mass;
        xl.custom_cleavable = self.cleavable;
        xl.stub_short = self.stub_short;
        xl.stub_long = self.stub_long;
        xl.stub_alt = self.stub_alt.clone();
        xl.asymmetric = self.asymmetric;
        xl.stub2_short = self.stub2_short;
        xl.stub2_long = self.stub2_long;
        xl.stub2_alt = self.stub2_alt.clone();
        xl.two_sites = self.two_sites;
        xl.residues = self.residues.clone();
        xl.protein_nterm = self.protein_nterm;
        xl.site1 = self.site1.clone();
        xl.site2 = self.site2.clone();
    }

    /// Is `xl` set up with this crosslinker?
    pub fn matches(&self, xl: &Crosslink) -> bool {
        xl.reagent == Reagent::Custom && self.same_chemistry(&SavedLinker::from_crosslink(&self.name, xl))
    }

    fn same_chemistry(&self, o: &SavedLinker) -> bool {
        let same = |a: f64, b: f64| (a - b).abs() < 1e-6;
        let lists = |a: &str, b: &str| match (settings::parse_list(a), settings::parse_list(b)) {
            (Ok(x), Ok(y)) => x.len() == y.len() && x.iter().zip(&y).all(|(p, q)| same(*p, *q)),
            _ => a == b,
        };
        let arm1 = same(self.stub_short, o.stub_short) && same(self.stub_long, o.stub_long) && lists(&self.stub_alt, &o.stub_alt);
        let arm2 =
            same(self.stub2_short, o.stub2_short) && same(self.stub2_long, o.stub2_long) && lists(&self.stub2_alt, &o.stub2_alt);
        let stubs = !self.cleavable || (arm1 && (!self.asymmetric || arm2));
        let sites = if self.two_sites {
            same_site(&self.site1, &o.site1) && same_site(&self.site2, &o.site2)
        } else {
            same_residues(&self.residues, &o.residues) && self.protein_nterm == o.protein_nterm
        };
        same(self.mass, o.mass)
            && self.cleavable == o.cleavable
            && self.asymmetric == o.asymmetric
            && self.two_sites == o.two_sites
            && stubs
            && sites
    }

    /// The stub masses, e.g. `54.01056 / 103.9932 (+ 85.98264)`, or `non-cleavable`
    /// (an asymmetric linker: one line per arm)
    pub fn stubs_text(&self) -> String {
        if !self.cleavable {
            return "non-cleavable".into();
        }
        let arm = |short: f64, long: f64, alt: &str| {
            let alt = if alt.trim().is_empty() { String::new() } else { format!(" (+ {})", alt.trim()) };
            format!("{} / {}{alt}", number(short), number(long))
        };
        let first = arm(self.stub_short, self.stub_long, &self.stub_alt);
        if self.asymmetric {
            format!("arm 1: {first}\narm 2: {}", arm(self.stub2_short, self.stub2_long, &self.stub2_alt))
        } else {
            first
        }
    }

    /// The reactive sites, e.g. `K, protein N-term`, or one line per site
    pub fn sites_text(&self) -> String {
        if self.two_sites {
            format!("Site 1: {}\nSite 2: {}", site_text(&self.site1), site_text(&self.site2))
        } else {
            site_text(&Site { residues: self.residues.clone(), protein_nterm: self.protein_nterm, protein_cterm: false })
        }
    }

    /// Mass, stubs and sites in one line
    pub fn summary(&self) -> String {
        let flat = |t: String| t.replace('\n', "; ");
        format!("{} Da, {}; reacts with {}", number(self.mass), flat(self.stubs_text()), flat(self.sites_text()))
    }

    /// In the form of Sage's `crosslink` settings, with a name
    fn to_json(&self) -> Value {
        let mut m = Map::new();
        m.insert("name".into(), json!(self.name));
        m.insert("mass".into(), json!(self.mass));
        if self.cleavable {
            let arm = |short: f64, long: f64, alt: &str| {
                json!({"short": short, "long": long, "alt": settings::parse_list(alt).unwrap_or_default()})
            };
            let mut stubs = arm(self.stub_short, self.stub_long, &self.stub_alt);
            if self.asymmetric {
                stubs["other"] = arm(self.stub2_short, self.stub2_long, &self.stub2_alt);
            }
            m.insert("cleavable".into(), stubs);
        }
        if self.two_sites {
            m.insert("site1".into(), json!(self.site1.to_spec()));
            m.insert("site2".into(), json!(self.site2.to_spec()));
        } else {
            m.insert("residues".into(), json!(self.residues));
            m.insert("protein_nterm".into(), json!(self.protein_nterm));
        }
        Value::Object(m)
    }

    fn from_json(v: &Value) -> Option<Self> {
        let name = v.get("name")?.as_str()?.trim().to_string();
        let mass = v.get("mass")?.as_f64()?;
        if name.is_empty() || !mass.is_finite() {
            return None;
        }
        let d = Crosslink::default();
        let list = |v: Option<&Value>| {
            v.and_then(Value::as_array)
                .map(|a| a.iter().filter_map(Value::as_f64).map(|x| x.to_string()).collect::<Vec<_>>().join(", "))
                .unwrap_or_default()
        };
        let f = |v: Option<&Value>, default: f64| v.and_then(Value::as_f64).unwrap_or(default);
        let mut l = SavedLinker::from_crosslink(&name, &d);
        l.mass = mass;
        l.cleavable = false;
        l.asymmetric = false;
        if let Some(stubs) = v.get("cleavable").filter(|s| s.is_object()) {
            l.cleavable = true;
            l.stub_short = f(stubs.get("short"), d.stub_short);
            l.stub_long = f(stubs.get("long"), d.stub_long);
            l.stub_alt = list(stubs.get("alt"));
            if let Some(other) = stubs.get("other").filter(|s| s.is_object()) {
                l.asymmetric = true;
                l.stub2_short = f(other.get("short"), d.stub2_short);
                l.stub2_long = f(other.get("long"), d.stub2_long);
                l.stub2_alt = list(other.get("alt"));
            }
        }
        if let (Some(s1), Some(s2)) = (
            v.get("site1").and_then(Value::as_str),
            v.get("site2").and_then(Value::as_str),
        ) {
            l.two_sites = true;
            l.site1 = Site::from_spec(s1);
            l.site2 = Site::from_spec(s2);
        } else {
            l.two_sites = false;
            l.residues = v.get("residues").and_then(Value::as_str).unwrap_or("K").trim().to_string();
            l.protein_nterm = v.get("protein_nterm").and_then(Value::as_bool).unwrap_or(true);
        }
        Some(l)
    }
}

fn same_residues(a: &str, b: &str) -> bool {
    let set = |s: &str| {
        let mut v: Vec<char> = s.chars().filter(char::is_ascii_alphabetic).map(|c| c.to_ascii_uppercase()).collect();
        v.sort_unstable();
        v.dedup();
        v
    };
    set(a) == set(b)
}

fn same_site(a: &Site, b: &Site) -> bool {
    same_residues(&a.residues, &b.residues) && a.protein_nterm == b.protein_nterm && a.protein_cterm == b.protein_cterm
}

fn site_text(site: &Site) -> String {
    let mut parts = Vec::new();
    if !site.residues.trim().is_empty() {
        parts.push(site.residues.trim().to_string());
    }
    if site.protein_nterm {
        parts.push("protein N-term".into());
    }
    if site.protein_cterm {
        parts.push("protein C-term".into());
    }
    if parts.is_empty() {
        "nothing".into()
    } else {
        parts.join(", ")
    }
}

/// A mass without trailing zeros, e.g. `138.06808`, `-18.010565`
fn number(x: f64) -> String {
    let s = format!("{x:.6}");
    let s = s.trim_end_matches('0').trim_end_matches('.');
    if s == "-0" {
        "0".into()
    } else {
        s.into()
    }
}

/// A modification saved under a name
#[derive(Clone, Debug, PartialEq)]
pub struct SavedMod {
    pub name: String,
    pub target: String,
    pub mass: f64,
    /// ProForma name written to the results (e.g. `UNIMOD:35`), may be empty
    pub proforma: String,
    pub kind: ModKind,
    /// Variable mods only: maximum occurrences per peptide, 0 = no limit
    pub max_count: u32,
}

impl SavedMod {
    pub fn from_modification(name: &str, m: &Modification) -> Self {
        Self {
            name: name.trim().to_string(),
            target: m.target.trim().to_string(),
            mass: m.mass,
            proforma: m.name.trim().to_string(),
            kind: m.kind,
            max_count: m.max_count,
        }
    }

    pub fn to_modification(&self) -> Modification {
        Modification {
            target: self.target.clone(),
            mass: self.mass,
            name: self.proforma.clone(),
            kind: self.kind,
            max_count: self.max_count,
        }
    }

    /// Is this modification already in the search (same kind, target and mass)?
    pub fn in_search(&self, mods: &[Modification]) -> bool {
        mods.iter().any(|m| m.kind == self.kind && m.target.trim() == self.target && (m.mass - self.mass).abs() < 1e-6)
    }

    pub fn summary(&self) -> String {
        let kind = if self.kind == ModKind::Fixed { "fixed" } else { "variable" };
        let mut s = format!("{:+.6} Da on {}, {kind}", self.mass, self.target);
        if self.kind == ModKind::Variable && self.max_count > 0 {
            s += &format!(", at most {} per peptide", self.max_count);
        }
        if !self.proforma.is_empty() {
            s += &format!(" ({})", self.proforma);
        }
        s
    }

    fn to_json(&self) -> Value {
        json!({
            "name": self.name,
            "target": self.target,
            "mass": self.mass,
            "proforma": self.proforma,
            "kind": if self.kind == ModKind::Fixed { "fixed" } else { "variable" },
            "max_count": self.max_count,
        })
    }

    fn from_json(v: &Value) -> Option<Self> {
        let name = v.get("name")?.as_str()?.trim().to_string();
        let target = v.get("target")?.as_str()?.trim().to_string();
        let mass = v.get("mass")?.as_f64()?;
        if name.is_empty() || target.is_empty() || !mass.is_finite() {
            return None;
        }
        Some(Self {
            name,
            target,
            mass,
            proforma: v.get("proforma").and_then(Value::as_str).unwrap_or("").trim().to_string(),
            kind: if v.get("kind").and_then(Value::as_str) == Some("fixed") { ModKind::Fixed } else { ModKind::Variable },
            max_count: v.get("max_count").and_then(Value::as_u64).unwrap_or(0).min(u32::MAX as u64) as u32,
        })
    }
}

/// Where the library is kept
pub fn library_path() -> Option<PathBuf> {
    let base = std::env::var_os("APPDATA").or_else(|| std::env::var_os("HOME"))?;
    Some(PathBuf::from(base).join("sage-gui").join("library.json"))
}

/// The saved crosslinkers and modifications
#[derive(Default)]
pub struct Library {
    pub crosslinkers: Vec<SavedLinker>,
    pub modifications: Vec<SavedMod>,
    /// The file it is kept in; `None` keeps it in memory only
    path: Option<PathBuf>,
    /// The file could not be read: it is backed up before it is replaced
    unreadable: bool,
    /// The last problem reading or writing the file
    pub error: Option<String>,
}

impl Library {
    /// The library of this computer
    pub fn load() -> Self {
        library_path().map(|p| Self::load_from(&p)).unwrap_or_default()
    }

    pub fn load_from(path: &Path) -> Self {
        let mut lib = Library { path: Some(path.to_path_buf()), ..Default::default() };
        let parsed = match std::fs::read_to_string(path) {
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => return lib,
            Err(e) => Err(e.to_string()),
            Ok(text) => serde_json::from_str::<Value>(&text).map_err(|e| e.to_string()),
        };
        match parsed {
            Ok(v) if v.is_object() => {
                let all = |key: &str| v.get(key).and_then(Value::as_array).cloned().unwrap_or_default();
                lib.crosslinkers = all("crosslinkers").iter().filter_map(SavedLinker::from_json).collect();
                lib.modifications = all("modifications").iter().filter_map(SavedMod::from_json).collect();
            }
            Ok(_) => lib.unreadable(path, "it is not a library"),
            Err(e) => lib.unreadable(path, &e),
        }
        lib
    }

    fn unreadable(&mut self, path: &Path, why: &str) {
        self.unreadable = true;
        self.error = Some(format!(
            "The saved crosslinkers and modifications could not be read from {} ({why}). \
             It is kept as {} when something is saved.",
            path.display(),
            backup_path(path).display()
        ));
    }

    pub fn to_json(&self) -> Value {
        json!({
            "version": 1,
            "crosslinkers": self.crosslinkers.iter().map(SavedLinker::to_json).collect::<Vec<_>>(),
            "modifications": self.modifications.iter().map(SavedMod::to_json).collect::<Vec<_>>(),
        })
    }

    /// Write the library to its file (after every change)
    fn save(&mut self) -> Result<(), String> {
        let Some(path) = self.path.clone() else { return Ok(()) };
        let result = self.write(&path);
        self.error = result.clone().err();
        result
    }

    fn write(&mut self, path: &Path) -> Result<(), String> {
        let failed = |e: std::io::Error| format!("Could not save {}: {e}", path.display());
        if let Some(dir) = path.parent() {
            std::fs::create_dir_all(dir).map_err(failed)?;
        }
        if self.unreadable && path.exists() {
            std::fs::copy(path, backup_path(path)).map_err(failed)?;
        }
        self.unreadable = false;
        // written whole, then put in place: a crash never leaves half a file
        let tmp = path.with_extension("json.tmp");
        let text = serde_json::to_string_pretty(&self.to_json()).map_err(|e| e.to_string())?;
        std::fs::write(&tmp, text).map_err(failed)?;
        std::fs::rename(&tmp, path).map_err(failed)
    }

    /// A name for a new entry, or a renamed one (`except` = its own index)
    fn check_name<'a>(names: impl Iterator<Item = &'a str>, name: &str, except: Option<usize>, what: &str) -> Result<(), String> {
        let name = name.trim();
        if name.is_empty() {
            return Err(format!("Give the {what} a name."));
        }
        match names.enumerate().find(|(i, n)| Some(*i) != except && n.eq_ignore_ascii_case(name)) {
            Some(_) => Err(format!("A saved {what} is already called \"{name}\": choose another name, or use Edit to change it.")),
            None => Ok(()),
        }
    }

    /// Save a new crosslinker (`index` None) or replace a saved one
    pub fn put_linker(&mut self, linker: SavedLinker, index: Option<usize>) -> Result<(), String> {
        Self::check_name(self.crosslinkers.iter().map(|l| l.name.as_str()), &linker.name, index, "crosslinker")?;
        match index {
            Some(i) if i < self.crosslinkers.len() => self.crosslinkers[i] = linker,
            _ => self.crosslinkers.push(linker),
        }
        self.save()
    }

    pub fn delete_linker(&mut self, index: usize) -> Result<(), String> {
        if index < self.crosslinkers.len() {
            self.crosslinkers.remove(index);
        }
        self.save()
    }

    /// The saved crosslinker `xl` is set up with, if any
    pub fn matching_linker(&self, xl: &Crosslink) -> Option<&SavedLinker> {
        self.crosslinkers.iter().find(|l| l.matches(xl))
    }

    /// Save a new modification (`index` None) or replace a saved one
    pub fn put_mod(&mut self, m: SavedMod, index: Option<usize>) -> Result<(), String> {
        Self::check_name(self.modifications.iter().map(|m| m.name.as_str()), &m.name, index, "modification")?;
        settings::check_target(&m.target)?;
        if !m.mass.is_finite() {
            return Err("The mass shift is not a number.".into());
        }
        match index {
            Some(i) if i < self.modifications.len() => self.modifications[i] = m,
            _ => self.modifications.push(m),
        }
        self.save()
    }

    pub fn delete_mod(&mut self, index: usize) -> Result<(), String> {
        if index < self.modifications.len() {
            self.modifications.remove(index);
        }
        self.save()
    }
}

fn backup_path(path: &Path) -> PathBuf {
    path.with_extension("json.bak")
}

/// Copy the linker chemistry (reagent, mass, stubs, reactive sites) from one
/// crosslink setting to another, leaving the search settings alone
fn copy_chemistry(from: &Crosslink, to: &mut Crosslink) {
    to.reagent = from.reagent;
    to.custom_mass = from.custom_mass;
    to.custom_cleavable = from.custom_cleavable;
    to.stub_short = from.stub_short;
    to.stub_long = from.stub_long;
    to.stub_alt = from.stub_alt.clone();
    to.asymmetric = from.asymmetric;
    to.stub2_short = from.stub2_short;
    to.stub2_long = from.stub2_long;
    to.stub2_alt = from.stub2_alt.clone();
    to.two_sites = from.two_sites;
    to.residues = from.residues.clone();
    to.protein_nterm = from.protein_nterm;
    to.site1 = from.site1.clone();
    to.site2 = from.site2.clone();
}

// ---------------------------------------------------------------------------
// Panels

/// A saved crosslinker being edited
pub struct LinkerEdit {
    pub index: usize,
    /// The name, possibly changed
    pub name: String,
    /// The crosslinker set up before Edit was clicked, restored by Cancel
    pub before: Crosslink,
}

/// State of the two panels
#[derive(Default)]
pub struct Editor {
    /// Name typed for "Save the current crosslinker as"
    pub linker_name: String,
    pub linker_edit: Option<LinkerEdit>,
    /// A saved crosslinker whose Delete waits for confirmation
    pub linker_delete: Option<usize>,
    pub linker_message: Option<(bool, String)>,
    /// A row of the modifications table being saved: its index and the name typed
    pub mod_save: Option<(usize, String)>,
    /// A saved modification being edited: its index and the changed copy
    pub mod_edit: Option<(usize, SavedMod)>,
    pub mod_delete: Option<usize>,
    pub mod_message: Option<(bool, String)>,
}

fn message(ui: &mut Ui, m: &Option<(bool, String)>, lib: &Library) {
    if let Some(e) = &lib.error {
        ui.label(RichText::new(format!("✖ {e}")).color(Color32::from_rgb(210, 60, 60)));
    }
    if let Some((error, text)) = m {
        if *error {
            ui.label(RichText::new(format!("✖ {text}")).color(Color32::from_rgb(210, 60, 60)));
        } else {
            ui.label(RichText::new(format!("✔ {text}")).weak());
        }
    }
}

fn outcome(result: Result<(), String>, done: String) -> Option<(bool, String)> {
    Some(match result {
        Ok(()) => (false, done),
        Err(e) => (true, e),
    })
}

/// The Reagent menu's saved crosslinkers; picking one sets it up
pub fn reagent_menu(ui: &mut Ui, xl: &mut Crosslink, lib: &Library) {
    if lib.crosslinkers.is_empty() {
        return;
    }
    ui.separator();
    ui.label(RichText::new("Saved crosslinkers").weak().small());
    for l in &lib.crosslinkers {
        if ui.selectable_label(l.matches(xl), &l.name).on_hover_text(l.summary()).clicked() {
            l.apply(xl);
        }
    }
}

/// The text the Reagent menu shows: a saved crosslinker's name when `xl` is set up with one
pub fn reagent_text(xl: &Crosslink, lib: &Library) -> String {
    match lib.matching_linker(xl) {
        Some(l) => format!("{} (saved)", l.name),
        None => xl.reagent.label().to_string(),
    }
}

/// "Saved crosslinkers": save the crosslinker set up above under a name, and
/// use, edit or delete saved ones
pub fn linker_panel(ui: &mut Ui, xl: &mut Crosslink, lib: &mut Library, ed: &mut Editor) {
    crate::heading(ui, "Saved crosslinkers");
    crate::note(
        ui,
        "Save the crosslinker set up above under a name: its linker mass, whether it is MS-cleavable, its stub \
         masses and its reactive sites. Saved crosslinkers are kept on this computer for later searches and are \
         also listed in the Reagent menu.",
    );

    if let Some(edit) = &mut ed.linker_edit {
        let original = lib.crosslinkers.get(edit.index).map(|l| l.name.clone()).unwrap_or_default();
        let mut save = false;
        let mut cancel = false;
        ui.group(|ui| {
            ui.label(theme::bold(format!("Editing the saved crosslinker \"{original}\"")));
            ui.label(RichText::new("Change its settings above, then save the changes here.").weak());
            ui.horizontal(|ui| {
                ui.label("Name");
                ui.add(TextEdit::singleline(&mut edit.name).desired_width(200.0));
                save = ui.button("Save changes").on_hover_text("Replace the saved crosslinker with the settings above.").clicked();
                cancel = ui.button("Cancel").on_hover_text("Keep the saved crosslinker as it was, and go back to the crosslinker set up before Edit.").clicked();
            });
        });
        if save {
            let linker = SavedLinker::from_crosslink(&edit.name, xl);
            let name = linker.name.clone();
            let result = lib.put_linker(linker, Some(edit.index));
            let ok = result.is_ok();
            ed.linker_message = outcome(result, format!("Saved the changes to \"{name}\"."));
            if ok {
                ed.linker_edit = None;
            }
        } else if cancel {
            copy_chemistry(&edit.before, xl);
            ed.linker_edit = None;
            ed.linker_message = None;
        }
    } else {
        ui.horizontal(|ui| {
            ui.label("Save the current crosslinker as");
            let field = ui.add(TextEdit::singleline(&mut ed.linker_name).hint_text("name, e.g. DSSO-d8").desired_width(200.0));
            let enter = field.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
            let save = ui
                .button("Save")
                .on_hover_text("Save the reagent set up above (mass, cleavable or not, stub masses, reactive sites) under this name.");
            if save.clicked() || (enter && !ed.linker_name.trim().is_empty()) {
                let linker = SavedLinker::from_crosslink(&ed.linker_name, xl);
                let name = linker.name.clone();
                let result = lib.put_linker(linker, None);
                let ok = result.is_ok();
                ed.linker_message = outcome(result, format!("Saved \"{name}\"."));
                if ok {
                    ed.linker_name.clear();
                }
            }
        });
    }
    message(ui, &ed.linker_message, lib);

    if lib.crosslinkers.is_empty() {
        return;
    }
    ui.add_space(4.0);
    let editing = ed.linker_edit.is_some();
    let mut action = None;
    egui::Grid::new("saved_linkers").num_columns(5).spacing([14.0, 6.0]).striped(true).show(ui, |ui| {
        ui.label(theme::bold("Name"));
        ui.label(theme::bold("Linker mass"));
        ui.label(theme::bold("Stub masses")).on_hover_text("Short / long stub (Da), with any additional stubs in brackets.");
        ui.label(theme::bold("Reactive sites"));
        ui.label("");
        ui.end_row();
        for (i, l) in lib.crosslinkers.iter().enumerate() {
            if l.matches(xl) {
                ui.label(RichText::new(&l.name).strong()).on_hover_text("The crosslinker set up above.");
            } else {
                ui.label(&l.name);
            }
            ui.label(format!("{} Da", number(l.mass)));
            ui.label(l.stubs_text());
            ui.label(l.sites_text());
            ui.horizontal(|ui| {
                if ed.linker_delete == Some(i) {
                    ui.label(format!("Delete \"{}\"?", l.name));
                    if ui.small_button("Delete").clicked() {
                        action = Some(("delete", i));
                    }
                    if ui.small_button("Keep").clicked() {
                        action = Some(("keep", i));
                    }
                } else {
                    let editing_this = ed.linker_edit.as_ref().is_some_and(|e| e.index == i);
                    if ui
                        .add_enabled(!editing, egui::Button::new("Use").small())
                        .on_hover_text("Set up this crosslinker for the search.")
                        .clicked()
                    {
                        action = Some(("use", i));
                    }
                    if ui
                        .add_enabled(!editing, egui::Button::new(if editing_this { "Editing…" } else { "Edit" }).small())
                        .on_hover_text("Set this crosslinker up above to change it, then Save changes.")
                        .clicked()
                    {
                        action = Some(("edit", i));
                    }
                    if ui
                        .add_enabled(!editing, egui::Button::new("Delete").small())
                        .on_hover_text("Remove this saved crosslinker (you are asked to confirm).")
                        .clicked()
                    {
                        action = Some(("ask", i));
                    }
                }
            });
            ui.end_row();
        }
    });
    let Some((what, i)) = action else { return };
    let name = lib.crosslinkers.get(i).map(|l| l.name.clone()).unwrap_or_default();
    match what {
        "use" => {
            lib.crosslinkers[i].apply(xl);
            ed.linker_message = Some((false, format!("Set up \"{name}\".")));
        }
        "edit" => {
            ed.linker_delete = None;
            ed.linker_edit = Some(LinkerEdit { index: i, name: name.clone(), before: xl.clone() });
            lib.crosslinkers[i].apply(xl);
            ed.linker_message = None;
        }
        "ask" => ed.linker_delete = Some(i),
        "keep" => ed.linker_delete = None,
        "delete" => {
            ed.linker_delete = None;
            ed.linker_message = outcome(lib.delete_linker(i), format!("Deleted \"{name}\"."));
        }
        _ => {}
    }
}

/// The "Add saved modification…" menu next to "Add common modification…"
pub fn mod_menu(ui: &mut Ui, mods: &mut Vec<Modification>, lib: &Library) {
    if lib.modifications.is_empty() {
        return;
    }
    ComboBox::from_id_salt("add_saved_mod")
        .selected_text("Add saved modification…")
        .width(240.0)
        .show_ui(ui, |ui| {
            for m in &lib.modifications {
                let added = m.in_search(mods);
                let label = if added { format!("{}  ✔", m.name) } else { m.name.clone() };
                if ui.selectable_label(false, label).on_hover_text(m.summary()).clicked() && !added {
                    mods.push(m.to_modification());
                }
            }
        })
        .response
        .on_hover_text("Add a modification you saved under a name (Saved modifications, below).");
}

/// A name for a modification of the search table: its Unimod title, or its mass and target
fn default_mod_name(m: &Modification) -> String {
    let unimod = crate::unimod::parse_accession(&m.name)
        .filter(|_| m.name.to_ascii_uppercase().starts_with("UNIMOD"))
        .and_then(|a| crate::unimod::bundled().get(&a));
    match unimod {
        Some(e) => format!("{} ({})", e.title, m.target.trim()),
        None => format!("{:+.4} Da ({})", m.mass, m.target.trim()),
    }
}

/// The Save button of a row of the modifications table
pub fn mod_save_button(ui: &mut Ui, i: usize, m: &Modification, ed: &mut Editor) {
    if ui.small_button("Save…").on_hover_text("Save this modification under a name, for later searches.").clicked() {
        ed.mod_save = Some((i, default_mod_name(m)));
        ed.mod_message = None;
    }
}

/// After a row of the modifications table was removed
pub fn mod_row_removed(ed: &mut Editor, i: usize) {
    let Some(j) = ed.mod_save.as_ref().map(|(j, _)| *j) else { return };
    if j == i {
        ed.mod_save = None;
    } else if let (true, Some((j, _))) = (j > i, ed.mod_save.as_mut()) {
        *j -= 1;
    }
}

/// Naming a row of the modifications table being saved, and the saved modifications
pub fn mod_panel(ui: &mut Ui, mods: &mut Vec<Modification>, lib: &mut Library, ed: &mut Editor) {
    if let Some((i, name)) = &mut ed.mod_save {
        let i = *i;
        let Some(m) = mods.get(i).cloned() else {
            ed.mod_save = None;
            return;
        };
        let mut save = false;
        let mut cancel = false;
        ui.add_space(4.0);
        ui.horizontal(|ui| {
            ui.label(format!("Save {:+.6} Da on {} as", m.mass, m.target.trim()));
            let field = ui.add(TextEdit::singleline(name).desired_width(220.0));
            let enter = field.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
            save = ui.button("Save").clicked() || enter;
            cancel = ui.button("Cancel").clicked();
        });
        if save {
            let saved = SavedMod::from_modification(name, &m);
            let shown = saved.name.clone();
            let result = lib.put_mod(saved, None);
            let ok = result.is_ok();
            ed.mod_message = outcome(result, format!("Saved \"{shown}\"."));
            if ok {
                ed.mod_save = None;
            }
        } else if cancel {
            ed.mod_save = None;
            ed.mod_message = None;
        }
    }

    crate::heading(ui, "Saved modifications");
    crate::note(
        ui,
        "Modifications saved under a name (Save…, in the table above) are kept on this computer for later \
         searches. Add one to this search with Add, or with \"Add saved modification…\" at the top.",
    );
    message(ui, &ed.mod_message, lib);
    if lib.modifications.is_empty() {
        ui.label(RichText::new("No saved modifications yet.").weak());
        return;
    }

    let editing = ed.mod_edit.is_some();
    let mut action = None;
    egui::Grid::new("saved_mods").num_columns(7).spacing([12.0, 6.0]).striped(true).show(ui, |ui| {
        ui.label(theme::bold("Name"));
        ui.label(theme::bold("Type"));
        ui.label(theme::bold("Target")).on_hover_text(crate::TARGET_TIP);
        ui.label(theme::bold("Mass shift"));
        ui.label(theme::bold("Max per peptide"));
        ui.label(theme::bold("ProForma name"));
        ui.label("");
        ui.end_row();
        for (i, saved) in lib.modifications.iter().enumerate() {
            match &mut ed.mod_edit {
                Some((j, draft)) if *j == i => {
                    ui.add(TextEdit::singleline(&mut draft.name).desired_width(150.0));
                    ComboBox::from_id_salt(("saved_mod_kind", i))
                        .width(90.0)
                        .selected_text(if draft.kind == ModKind::Fixed { "Fixed" } else { "Variable" })
                        .show_ui(ui, |ui| {
                            ui.selectable_value(&mut draft.kind, ModKind::Fixed, "Fixed");
                            ui.selectable_value(&mut draft.kind, ModKind::Variable, "Variable");
                        });
                    let bad = settings::check_target(&draft.target).is_err();
                    ui.add(
                        TextEdit::singleline(&mut draft.target)
                            .desired_width(50.0)
                            .text_color_opt(bad.then_some(Color32::RED)),
                    )
                    .on_hover_text(crate::TARGET_TIP);
                    crate::mass(ui, &mut draft.mass);
                    ui.add_enabled_ui(draft.kind == ModKind::Variable, |ui| {
                        theme::number(ui, &mut draft.max_count, 0..=10, 1.0, |d| d)
                    });
                    ui.add(TextEdit::singleline(&mut draft.proforma).hint_text("UNIMOD:…").desired_width(110.0));
                    ui.horizontal(|ui| {
                        if ui.small_button("Save").on_hover_text("Keep these changes.").clicked() {
                            action = Some(("save", i));
                        }
                        if ui.small_button("Cancel").on_hover_text("Leave the saved modification as it was.").clicked() {
                            action = Some(("cancel", i));
                        }
                    });
                }
                _ => {
                    ui.label(&saved.name);
                    ui.label(if saved.kind == ModKind::Fixed { "Fixed" } else { "Variable" });
                    ui.monospace(&saved.target);
                    ui.label(format!("{:+.6} Da", saved.mass));
                    ui.label(match (saved.kind, saved.max_count) {
                        (ModKind::Fixed, _) => "-".to_string(),
                        (_, 0) => "no limit".to_string(),
                        (_, n) => n.to_string(),
                    });
                    ui.label(if saved.proforma.is_empty() { "-" } else { saved.proforma.as_str() });
                    ui.horizontal(|ui| {
                        if ed.mod_delete == Some(i) {
                            ui.label(format!("Delete \"{}\"?", saved.name));
                            if ui.small_button("Delete").clicked() {
                                action = Some(("delete", i));
                            }
                            if ui.small_button("Keep").clicked() {
                                action = Some(("keep", i));
                            }
                        } else {
                            let added = saved.in_search(mods);
                            if ui
                                .add_enabled(!editing && !added, egui::Button::new(if added { "✔ added" } else { "Add" }).small())
                                .on_hover_text("Add this modification to the search.")
                                .clicked()
                            {
                                action = Some(("add", i));
                            }
                            if ui
                                .add_enabled(!editing, egui::Button::new("Edit").small())
                                .on_hover_text("Change this saved modification.")
                                .clicked()
                            {
                                action = Some(("edit", i));
                            }
                            if ui
                                .add_enabled(!editing, egui::Button::new("Delete").small())
                                .on_hover_text("Remove this saved modification (you are asked to confirm).")
                                .clicked()
                            {
                                action = Some(("ask", i));
                            }
                        }
                    });
                }
            }
            ui.end_row();
        }
    });
    let Some((what, i)) = action else { return };
    let name = lib.modifications.get(i).map(|m| m.name.clone()).unwrap_or_default();
    match what {
        "add" => {
            mods.push(lib.modifications[i].to_modification());
            ed.mod_message = Some((false, format!("Added \"{name}\" to the search.")));
        }
        "edit" => {
            ed.mod_delete = None;
            ed.mod_edit = Some((i, lib.modifications[i].clone()));
            ed.mod_message = None;
        }
        "save" => {
            if let Some((_, draft)) = ed.mod_edit.clone() {
                let shown = draft.name.trim().to_string();
                let result = lib.put_mod(draft, Some(i));
                let ok = result.is_ok();
                ed.mod_message = outcome(result, format!("Saved the changes to \"{shown}\"."));
                if ok {
                    ed.mod_edit = None;
                }
            }
        }
        "cancel" => {
            ed.mod_edit = None;
            ed.mod_message = None;
        }
        "ask" => ed.mod_delete = Some(i),
        "keep" => ed.mod_delete = None,
        "delete" => {
            ed.mod_delete = None;
            ed.mod_message = outcome(lib.delete_mod(i), format!("Deleted \"{name}\"."));
        }
        _ => {}
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn temp_file(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("xnova_library_{name}_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir.join("library.json")
    }

    fn with(reagent: Reagent) -> Crosslink {
        let mut xl = Crosslink { enabled: true, reagent, ..Default::default() };
        if let Some((s1, s2)) = reagent.sites() {
            xl.two_sites = true;
            xl.site1 = s1;
            xl.site2 = s2;
        }
        xl.reset_stubs();
        xl
    }

    #[test]
    fn presets_save_as_their_chemistry() {
        // DSSO: cleavable, its stubs, K + protein N-term
        let dsso = SavedLinker::from_crosslink(" My DSSO ", &with(Reagent::Dsso));
        assert_eq!(dsso.name, "My DSSO");
        assert!((dsso.mass - 158.00376).abs() < 1e-9 && dsso.cleavable && !dsso.asymmetric);
        assert_eq!(dsso.stubs_text(), "54.01056 / 103.9932 (+ 85.98264)");
        assert_eq!(dsso.sites_text(), "K, protein N-term");
        // set up again, it is a custom reagent with the same mass and stubs
        let mut xl = with(Reagent::Dss);
        dsso.apply(&mut xl);
        assert_eq!(xl.reagent, Reagent::Custom);
        assert!(xl.is_cleavable() && (xl.linker_mass() - 158.00376).abs() < 1e-9);
        assert!(dsso.matches(&xl));
        xl.custom_mass += 1.0;
        assert!(!dsso.matches(&xl), "a changed mass is no longer the saved crosslinker");

        // cliXlink: asymmetric
        let clix = SavedLinker::from_crosslink("clix", &with(Reagent::Clixlink));
        assert!(clix.asymmetric);
        assert!(clix.stubs_text().starts_with("arm 1: 54.01056 / 103.9932"), "{}", clix.stubs_text());
        // EDC: zero-length, two sites
        let edc = SavedLinker::from_crosslink("edc", &with(Reagent::Edc));
        assert!(!edc.cleavable && edc.two_sites);
        assert_eq!(edc.sites_text(), "Site 1: K, protein N-term\nSite 2: DE, protein C-term");
        assert!(edc.summary().ends_with("reacts with Site 1: K, protein N-term; Site 2: DE, protein C-term"));
        assert_eq!(edc.stubs_text(), "non-cleavable");
        // a preset is not "the saved one", even with the same chemistry
        assert!(!SavedLinker::from_crosslink("dss", &with(Reagent::Dss)).matches(&with(Reagent::Dss)));
    }

    #[test]
    fn the_search_settings_of_a_saved_linker_are_its_chemistry() {
        let mut custom = with(Reagent::Custom);
        custom.custom_mass = 200.1;
        custom.custom_cleavable = true;
        custom.stub_short = 50.0;
        custom.stub_long = 100.0;
        custom.stub_alt = "70.5".into();
        custom.residues = "KSTY".into();
        let saved = SavedLinker::from_crosslink("x", &custom);
        let mut xl = with(Reagent::Dss);
        xl.candidates = 77;
        xl.multivalent = true;
        saved.apply(&mut xl);
        // the Sage settings written are those of the custom linker
        let mut a = Settings::default();
        a.crosslink = custom.clone();
        let mut b = Settings::default();
        b.crosslink = xl.clone();
        let (ja, jb) = (a.to_json(&[]), b.to_json(&[]));
        for key in ["mass", "cleavable", "residues", "protein_nterm"] {
            assert_eq!(ja["crosslink"][key], jb["crosslink"][key], "{key}");
        }
        // and the search settings stay
        assert_eq!((xl.candidates, xl.multivalent), (77, true));
    }

    use crate::settings::Settings;

    #[test]
    fn library_round_trips_through_its_file() {
        let path = temp_file("round_trip");
        let mut lib = Library::load_from(&path);
        assert!(lib.crosslinkers.is_empty() && lib.error.is_none());
        lib.put_linker(SavedLinker::from_crosslink("DSSO-d8", &with(Reagent::Dsso)), None).unwrap();
        lib.put_linker(SavedLinker::from_crosslink("clix", &with(Reagent::Clixlink)), None).unwrap();
        lib.put_linker(SavedLinker::from_crosslink("EDC", &with(Reagent::Edc)), None).unwrap();
        let m = Modification { target: "K".into(), mass: 42.010565, name: "UNIMOD:1".into(), kind: ModKind::Variable, max_count: 2 };
        lib.put_mod(SavedMod::from_modification("Acetyl (K)", &m), None).unwrap();

        let back = Library::load_from(&path);
        assert_eq!(back.crosslinkers, lib.crosslinkers);
        assert_eq!(back.modifications, lib.modifications);
        let m2 = back.modifications[0].to_modification();
        assert_eq!((m2.target.as_str(), m2.mass, m2.name.as_str(), m2.kind, m2.max_count), ("K", 42.010565, "UNIMOD:1", ModKind::Variable, 2));

        // edit and delete are written at once
        let mut lib = back;
        let mut changed = lib.crosslinkers[0].clone();
        changed.name = "DSSO heavy".into();
        changed.mass = 166.05;
        lib.put_linker(changed, Some(0)).unwrap();
        lib.delete_linker(1).unwrap();
        lib.delete_mod(0).unwrap();
        let back = Library::load_from(&path);
        let names: Vec<_> = back.crosslinkers.iter().map(|l| l.name.as_str()).collect();
        assert_eq!(names, ["DSSO heavy", "EDC"]);
        assert!((back.crosslinkers[0].mass - 166.05).abs() < 1e-9);
        assert!(back.modifications.is_empty());
        let _ = std::fs::remove_dir_all(path.parent().unwrap());
    }

    #[test]
    fn names_are_required_and_unique() {
        let mut lib = Library::default();
        let dss = SavedLinker::from_crosslink("DSS", &with(Reagent::Dss));
        lib.put_linker(dss.clone(), None).unwrap();
        let mut empty = dss.clone();
        empty.name = "  ".into();
        assert!(lib.put_linker(empty, None).unwrap_err().contains("name"));
        let mut same = dss.clone();
        same.name = "dss".into();
        assert!(lib.put_linker(same.clone(), None).unwrap_err().contains("already called"));
        // its own name is fine when it is edited
        lib.put_linker(same, Some(0)).unwrap();
        assert_eq!(lib.crosslinkers.len(), 1);
        // modifications: the target is checked too
        let m = Modification { target: "K".into(), mass: 1.0, name: String::new(), kind: ModKind::Fixed, max_count: 0 };
        let mut bad = SavedMod::from_modification("x", &m);
        bad.target = "KK".into();
        assert!(lib.put_mod(bad, None).is_err());
        assert!(lib.modifications.is_empty());
    }

    #[test]
    fn an_unreadable_file_is_kept_as_a_backup() {
        let path = temp_file("unreadable");
        std::fs::write(&path, "{ not json").unwrap();
        let mut lib = Library::load_from(&path);
        assert!(lib.error.as_deref().unwrap().contains("could not be read"));
        assert!(lib.crosslinkers.is_empty());
        lib.put_linker(SavedLinker::from_crosslink("DSS", &with(Reagent::Dss)), None).unwrap();
        assert_eq!(std::fs::read_to_string(backup_path(&path)).unwrap(), "{ not json");
        assert_eq!(Library::load_from(&path).crosslinkers.len(), 1);
        assert!(lib.error.is_none());
        let _ = std::fs::remove_dir_all(path.parent().unwrap());
    }

    #[test]
    fn unimod_names_for_saved_mods() {
        let m = |name: &str, target: &str, mass: f64| Modification {
            target: target.into(),
            mass,
            name: name.into(),
            kind: ModKind::Variable,
            max_count: 0,
        };
        assert_eq!(default_mod_name(&m("UNIMOD:35", "M", 15.994915)), "Oxidation (M)");
        assert_eq!(default_mod_name(&m("", "K", 156.078644)), "+156.0786 Da (K)");
    }
}
