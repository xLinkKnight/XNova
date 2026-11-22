//! The end of a search, written for people: its settings as a list, its
//! spectrum and output files, and the quality-control tables. results.json
//! holds the same as JSON, for programs.
//!
//! Every setting is listed: the ones known here in words, any other one (a
//! new setting, an optional block) under its own name, so nothing in
//! results.json is left out.

use sage_core::crosslink::{takes_any_residue, CrosslinkSettings, Reagent};
use serde_json::{Map, Value};

/// Width of the settings' labels
const LABEL: usize = 26;

/// Log the summary of a finished search
pub fn to_log(params: &crate::input::Search) {
    // through JSON text, as results.json has it: f32 values keep their short
    // form (0.7, not 0.699999988079071)
    let value = serde_json::to_string(params).ok().and_then(|s| serde_json::from_str::<Value>(&s).ok());
    if let Some(value) = value {
        for line in lines(&value) {
            log::info!("{line}");
        }
    }
}

/// The summary's lines, from the parameters as written to results.json
pub fn lines(params: &Value) -> Vec<String> {
    let mut p = params.as_object().cloned().unwrap_or_default();
    let mut out = Vec::new();
    out.push("=== summary of the search (results.json has the same as JSON) ===".to_string());

    let spectra = strings(take(&mut p, "mzml_paths"));
    let outputs = strings(take(&mut p, "output_paths"));
    let qc = take(&mut p, "qc");

    out.push("Settings".into());
    let mut list = List::default();
    if let Some(v) = take(&mut p, "version") {
        list.item("XNova version", text(&v));
    }
    if let Some(Value::Object(db)) = take(&mut p, "database") {
        database(&mut list, db);
    }
    matching(&mut list, &mut p);
    if let Some(Value::Object(quant)) = take(&mut p, "quant") {
        quantification(&mut list, quant);
    }
    if let Some(Value::Object(xl)) = take(&mut p, "crosslink") {
        crosslinker(&mut list, xl);
    }
    let bruker = take(&mut p, "bruker_config");
    if spectra.iter().any(|s| is_bruker(s)) {
        if let Some(Value::Object(b)) = bruker {
            timstof(&mut list, b);
        }
    }
    if let Some(Value::Object(lib)) = take(&mut p, "spectral_library") {
        libraries(&mut list, lib);
    }
    if !p.is_empty() {
        list.heading("Other settings");
        rest(&mut list, p);
    }
    out.extend(list.lines);

    files(&mut out, "Spectrum files", &spectra, false);
    files(&mut out, "Output files", &outputs, true);

    if let Some(Value::Object(qc)) = qc {
        checks(&mut out, qc);
    }
    out
}

/// What one end of a linker reacts with, in words: `KSTY[` is "KSTY and
/// protein N-termini", every amino acid "any residue"
pub fn site_text(spec: &str) -> String {
    let residues: String = spec.chars().filter(char::is_ascii_uppercase).collect();
    let mut text = if takes_any_residue(&residues) {
        "any residue".to_string()
    } else {
        residues
    };
    if spec.contains('[') {
        text.push_str(" and protein N-termini");
    }
    if spec.contains(']') {
        text.push_str(" and protein C-termini");
    }
    text
}

/// The settings, as labels and values under headings
#[derive(Default)]
struct List {
    lines: Vec<String>,
}

impl List {
    fn heading(&mut self, text: &str) {
        self.lines.push(format!("  {text}"));
    }

    fn item(&mut self, label: &str, value: String) {
        self.items(label, vec![value]);
    }

    /// One value per line, the first beside the label
    fn items(&mut self, label: &str, values: Vec<String>) {
        for (i, value) in values.into_iter().enumerate() {
            let label = if i == 0 { label } else { "" };
            self.lines.push(format!("    {label:<LABEL$} {value}").trim_end().to_string());
        }
    }
}

fn take(m: &mut Map<String, Value>, key: &str) -> Option<Value> {
    m.remove(key)
}

/// A key as a label: `min_ion_index` is "min ion index"
fn label(key: &str) -> String {
    key.replace('_', " ")
}

/// A number as JSON wrote it, without a trailing `.0`
fn num(v: &Value) -> String {
    match v {
        Value::Number(n) => {
            let s = n.to_string();
            s.strip_suffix(".0").map(str::to_string).unwrap_or(s)
        }
        other => text(other),
    }
}

/// A mass shift with its sign
fn signed(v: &Value) -> String {
    let s = num(v);
    if s.starts_with('-') {
        s
    } else {
        format!("+{s}")
    }
}

/// A fraction as a percentage: 0.01 is "1%"
fn percent(v: &Value) -> String {
    match v.as_f64() {
        Some(x) => percent_text(x),
        None => text(v),
    }
}

/// A fraction as a percentage: 0.01 is "1%", 0.005 is "0.5%"
pub fn percent_text(x: impl Into<f64>) -> String {
    format!("{}%", trim(&format!("{:.2}", x.into() * 100.0)))
}

/// A mass held as f32, to 4 decimals as the run's header gives it
fn mass4(m: f32) -> String {
    trim(&format!("{m:.4}"))
}

/// `5.00` is "5", `1.50` is "1.5"
fn trim(s: &str) -> String {
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        s.to_string()
    }
}

fn yes(v: Option<Value>) -> String {
    match v {
        Some(Value::Bool(true)) => "yes".into(),
        Some(Value::Bool(false)) => "no".into(),
        Some(other) => text(&other),
        None => "no".into(),
    }
}

/// Any value on one line: no braces, quotes or indents
fn text(v: &Value) -> String {
    match v {
        Value::Null => "none".into(),
        Value::Bool(b) => if *b { "yes" } else { "no" }.into(),
        Value::Number(_) => num(v),
        Value::String(s) => path(s),
        Value::Array(a) if a.is_empty() => "none".into(),
        Value::Array(a) => {
            let sep = if a.iter().all(|x| !x.is_object() && !x.is_array()) { ", " } else { "; " };
            a.iter().map(text).collect::<Vec<_>>().join(sep)
        }
        Value::Object(m) => m
            .iter()
            .map(|(k, v)| match v {
                Value::Object(_) => format!("{} ({})", label(k), text(v)),
                _ => format!("{} {}", label(k), text(v)),
            })
            .collect::<Vec<_>>()
            .join(", "),
    }
}

/// A `file://` URL as a path; anything else as it is
fn path(s: &str) -> String {
    url::Url::parse(s)
        .ok()
        .filter(|u| u.scheme() == "file")
        .and_then(|u| u.to_file_path().ok())
        .map(|p| p.display().to_string())
        .unwrap_or_else(|| s.to_string())
}

fn strings(v: Option<Value>) -> Vec<String> {
    match v {
        Some(Value::Array(a)) => a.iter().map(|x| x.as_str().map(str::to_string).unwrap_or_else(|| text(x))).collect(),
        _ => Vec::new(),
    }
}

fn is_bruker(url: &str) -> bool {
    let u = url.trim_end_matches('/').to_ascii_lowercase();
    u.ends_with(".d") || u.ends_with(".tdf")
}

/// The leftover keys of a block, each under its own name
fn rest(list: &mut List, m: Map<String, Value>) {
    for (k, v) in m {
        list.item(&label(&k), text(&v));
    }
}

/// `{"ppm": [-15, 15]}` is "-15 to 15 ppm"
fn tolerance(v: &Value) -> String {
    if let Some((unit, Value::Array(range))) = v.as_object().and_then(|m| m.iter().next()) {
        if range.len() == 2 {
            let unit = match unit.as_str() {
                "ppm" => "ppm",
                "da" => "Da",
                "pct" => "%",
                other => other,
            };
            return format!("{} to {} {unit}", num(&range[0]), num(&range[1]));
        }
    }
    text(v)
}

/// `[2, 6]` is "2 to 6"
fn range(v: &Value) -> String {
    match v.as_array() {
        Some(a) if a.len() == 2 => format!("{} to {}", num(&a[0]), num(&a[1])),
        _ => text(v),
    }
}

/// A modification target in words: `^Q` is "peptide N-term Q"
fn target(t: &str) -> String {
    let mut chars = t.chars();
    let end = match chars.next() {
        Some('^') => "peptide N-term",
        Some('$') => "peptide C-term",
        Some('[') => "protein N-term",
        Some(']') => "protein C-term",
        _ => return t.to_string(),
    };
    let residue: String = chars.collect();
    if residue.is_empty() {
        end.to_string()
    } else {
        format!("{end} {residue}")
    }
}

/// One modification: its target, mass, name, and anything else it has
fn modification(t: &str, v: &Value, variable: bool) -> String {
    let Some(m) = v.as_object() else {
        return format!("{} {}", target(t), signed(v));
    };
    let mut m = m.clone();
    let mut s = target(t);
    if let Some(mass) = m.remove("mass") {
        s.push(' ');
        s.push_str(&signed(&mass));
    }
    match m.remove("name") {
        Some(Value::String(name)) if !name.is_empty() => s.push_str(&format!(" ({name})")),
        _ => {}
    }
    if let Some(count) = m.remove("max_count") {
        if variable && count.as_u64().is_some_and(|c| c > 0) {
            s.push_str(&format!(", up to {} per peptide", num(&count)));
        }
    }
    for (k, v) in m {
        s.push_str(&format!(", {} {}", label(&k), text(&v)));
    }
    s
}

fn modifications(v: Option<Value>, variable: bool) -> Vec<String> {
    let Some(Value::Object(m)) = v else { return vec!["none".into()] };
    let mut all = Vec::new();
    for (t, mods) in &m {
        match mods {
            Value::Array(list) => all.extend(list.iter().map(|x| modification(t, x, variable))),
            other => all.push(modification(t, other, variable)),
        }
    }
    if all.is_empty() {
        all.push("none".into());
    }
    all
}

fn database(list: &mut List, mut db: Map<String, Value>) {
    list.heading("Database");
    if let Some(v) = take(&mut db, "fasta") {
        list.item("FASTA", text(&v));
    }
    if let Some(Value::Object(mut e)) = take(&mut db, "enzyme") {
        let mut parts = Vec::new();
        match take(&mut e, "enzymes") {
            Some(Value::Array(names)) if !names.is_empty() => {
                parts.push(names.iter().map(text).collect::<Vec<_>>().join(" + "));
            }
            _ => {}
        }
        if let Some(Value::String(at)) = take(&mut e, "cleave_at") {
            let side = if e.get("c_terminal").and_then(Value::as_bool) == Some(false) { "before" } else { "after" };
            parts.push(if at.is_empty() { "no cleavage rule".into() } else { format!("cleaves {side} {at}") });
        }
        if let Some(Value::String(r)) = take(&mut e, "restrict") {
            parts.push(format!("not next to {r}"));
        }
        take(&mut e, "c_terminal");
        if let Some(n) = take(&mut e, "missed_cleavages") {
            parts.push(format!("up to {} missed cleavages", num(&n)));
        }
        let (lo, hi) = (take(&mut e, "min_len"), take(&mut e, "max_len"));
        if let (Some(lo), Some(hi)) = (&lo, &hi) {
            parts.push(format!("{} to {} residues", num(lo), num(hi)));
        }
        if let Some(Value::Bool(semi)) = take(&mut e, "semi_enzymatic") {
            if semi {
                parts.push("semi-enzymatic".into());
            }
        }
        for (k, v) in e {
            if !v.is_null() {
                parts.push(format!("{} {}", label(&k), text(&v)));
            }
        }
        list.item("enzyme", parts.join(", "));
    }
    if let (Some(lo), Some(hi)) = (take(&mut db, "peptide_min_mass"), take(&mut db, "peptide_max_mass")) {
        list.item("peptide mass", format!("{} to {} Da", num(&lo), num(&hi)));
    }
    if let Some(kinds) = take(&mut db, "ion_kinds") {
        let from = take(&mut db, "min_ion_index").map(|i| format!(", from ion {}", num(&i))).unwrap_or_default();
        list.item("fragment ions", format!("{}{from}", text(&kinds)));
    }
    list.items("fixed modifications", modifications(take(&mut db, "static_mods"), false));
    list.items("variable modifications", modifications(take(&mut db, "variable_mods"), true));
    if let Some(n) = take(&mut db, "max_variable_mods") {
        let mut s = format!("up to {} per peptide", num(&n));
        if let Some(c) = take(&mut db, "max_combinations").filter(|c| !c.is_null()) {
            s.push_str(&format!(", at most {} combinations", num(&c)));
        }
        list.item("variable mods", s);
    }
    let tag = take(&mut db, "decoy_tag").map(|t| text(&t)).unwrap_or_default();
    match take(&mut db, "generate_decoys") {
        Some(Value::Bool(true)) => list.item("decoys", format!("reversed proteins, tag {tag}")),
        Some(_) => list.item("decoys", format!("from the FASTA, tag {tag}")),
        None => {}
    }
    if let Some(v) = take(&mut db, "clip_n_term_met") {
        list.item("clip N-terminal Met", yes(Some(v)));
    }
    let (pre, chunk, low) = (take(&mut db, "prefilter"), take(&mut db, "prefilter_chunk_size"), take(&mut db, "prefilter_low_memory"));
    match pre {
        Some(Value::Bool(true)) => list.item(
            "prefilter",
            format!(
                "on, {} proteins per chunk, low memory {}",
                chunk.as_ref().map(num).unwrap_or_default(),
                yes(low)
            ),
        ),
        Some(_) => list.item("prefilter", "off".into()),
        None => {}
    }
    if let Some(v) = take(&mut db, "bucket_size") {
        list.item("fragment bucket size", num(&v));
    }
    rest(list, db);
}

/// The settings of the spectra and of matching them
fn matching(list: &mut List, p: &mut Map<String, Value>) {
    list.heading("Spectra and matching");
    if let Some(v) = take(p, "precursor_tol") {
        list.item("precursor tolerance", tolerance(&v));
    }
    if let Some(v) = take(p, "precursor_charge") {
        let override_ = take(p, "override_precursor_charge").and_then(|v| v.as_bool()).unwrap_or(false);
        let when = if override_ { "for every spectrum (the files' charges ignored)" } else { "where a spectrum has no charge" };
        list.item("precursor charge", format!("{}, {when}", range(&v)));
    }
    if let Some(v) = take(p, "isotope_errors") {
        let none = v.as_array().is_some_and(|a| a.iter().all(|x| x.as_f64() == Some(0.0)));
        list.item("isotope errors", if none { "none".into() } else { range(&v) });
    }
    if let Some(v) = take(p, "fragment_tol") {
        list.item("fragment tolerance", tolerance(&v));
    }
    if let Some(v) = take(p, "max_fragment_charge") {
        list.item("fragment charge", if v.is_null() { "up to the precursor's - 1".into() } else { format!("up to {}", num(&v)) });
    }
    if let (Some(lo), Some(hi)) = (take(p, "min_peaks"), take(p, "max_peaks")) {
        list.item("peaks per spectrum", format!("{} to {} (the most intense)", num(&lo), num(&hi)));
    }
    if let Some(v) = take(p, "min_matched_peaks") {
        list.item("matched peaks", format!("at least {}", num(&v)));
    }
    for (key, name) in [("deisotope", "deisotoping"), ("chimera", "chimeric spectra"), ("wide_window", "wide window")] {
        if let Some(v) = take(p, key) {
            list.item(name, yes(Some(v)));
        }
    }
    if let Some(v) = take(p, "report_psms") {
        list.item("matches reported", format!("{} per spectrum", num(&v)));
    }
    if let Some(v) = take(p, "score_type") {
        list.item("score", text(&v));
    }
    if let Some(v) = take(p, "predict_rt") {
        list.item("retention time prediction", yes(Some(v)));
    }
    if let Some(v) = take(p, "protein_grouping") {
        list.heading("Proteins");
        let fdr = take(p, "protein_grouping_peptide_fdr").map(|f| format!(", from peptides at {} FDR", percent(&f)));
        let on = v.as_bool() == Some(true);
        list.item("protein groups", if on { format!("yes{}", fdr.unwrap_or_default()) } else { "no".into() });
    }
}

fn quantification(list: &mut List, mut q: Map<String, Value>) {
    list.heading("Quantification");
    let settings = take(&mut q, "tmt_settings");
    match take(&mut q, "tmt") {
        Some(Value::Null) | None => list.item("TMT", "off".into()),
        Some(tmt) => {
            let s = settings.map(|s| format!(" ({})", text(&s))).unwrap_or_default();
            list.item("TMT", format!("{}{s}", text(&tmt)));
        }
    }
    let settings = take(&mut q, "lfq_settings");
    match take(&mut q, "lfq") {
        Some(Value::Bool(true)) => list.item("LFQ", settings.map(|s| text(&s)).unwrap_or_else(|| "on".into())),
        _ => list.item("LFQ", "off".into()),
    }
    rest(list, q);
}

fn crosslinker(list: &mut List, mut xl: Map<String, Value>) {
    list.heading("Crosslinker");
    let typed: Option<CrosslinkSettings> = serde_json::from_value(Value::Object(xl.clone())).ok();
    let reagent: Option<Reagent> = take(&mut xl, "reagent").and_then(|r| serde_json::from_value(r).ok());
    take(&mut xl, "mass");
    let cleavable = take(&mut xl, "cleavable");
    if let Some(s) = &typed {
        let name = reagent.map(|r| r.name()).unwrap_or("custom");
        let kind = if s.is_cleavable() { "MS-cleavable" } else { "non-cleavable" };
        list.item("reagent", format!("{name}, {} Da, {kind}", mass4(s.linker_mass())));
        if let Some(stubs) = s.stubs() {
            let arm = |short: f32, long: f32, alt: &[f32]| {
                let mut t = format!("{} and {} Da", mass4(short), mass4(long));
                if !alt.is_empty() {
                    let alt: Vec<String> = alt.iter().map(|&m| mass4(m)).collect();
                    t.push_str(&format!(", also {}", alt.join(", ")));
                }
                t
            };
            let mut stub_lines = vec![arm(stubs.short, stubs.long, &stubs.alt)];
            if let Some(o) = &stubs.other {
                stub_lines[0] = format!("first arm {}", stub_lines[0]);
                stub_lines.push(format!("second arm {}", arm(o.short, o.long, &o.alt)));
            }
            list.items("stubs", stub_lines);
        }
        let sites = match s.hetero_sites() {
            Some((s1, s2)) => format!("site 1: {}; site 2: {}", site_text(&s1), site_text(&s2)),
            None => site_text(&format!("{}{}", s.residues, if s.protein_nterm { "[" } else { "" })),
        };
        list.item("reacts with", sites);
    } else {
        if let Some(r) = reagent {
            list.item("reagent", r.name().into());
        }
        if let Some(c) = cleavable {
            list.item("stubs", text(&c));
        }
    }
    for key in ["residues", "protein_nterm", "site1", "site2"] {
        take(&mut xl, key);
    }
    let flag = |xl: &mut Map<String, Value>, key: &str| take(xl, key).map(|v| yes(Some(v)));
    if let Some(v) = flag(&mut xl, "allow_cterm_residue") {
        list.item("peptide C-terminal links", v);
    }
    if let Some(v) = flag(&mut xl, "loop_links") {
        list.item("loop-links", v);
    }
    if let Some(v) = take(&mut xl, "multivalent") {
        list.item("multivalent", if v.is_null() { "no".into() } else { format!("up to {} linkers", num(&v)) });
    }
    if let Some(v) = take(&mut xl, "candidates") {
        list.item("seed peptides", format!("{} per spectrum", num(&v)));
    }
    if let Some(v) = take(&mut xl, "min_chain_peaks") {
        list.item("peaks per chain", format!("at least {}", num(&v)));
    }
    if let Some(v) = take(&mut xl, "min_linked_peaks") {
        list.item("linked ions per chain", format!("at least {} (non-cleavable)", num(&v)));
    }
    if let Some(v) = flag(&mut xl, "require_complete_doublets") {
        list.item("complete doublets", v);
    }
    if let Some(v) = flag(&mut xl, "precursor_window") {
        list.item("precursor mass window", v);
    }
    if let Some(v) = take(&mut xl, "backend") {
        let survivors = take(&mut xl, "chain_survivors").map(|s| format!(" ({} chain survivors on reflect/gpu)", num(&s)));
        list.item("backend", format!("{}{}", text(&v), survivors.unwrap_or_default()));
    }
    if let Some(v) = flag(&mut xl, "linear_search") {
        list.item("linear search too", v);
    }
    let fdr = take(&mut xl, "linear_protein_fdr");
    if let Some(v) = take(&mut xl, "linear_proteins") {
        let narrowed = v.as_bool() == Some(true);
        let cutoff = fdr.as_ref().map(percent).unwrap_or_else(|| "1%".into());
        list.item(
            "crosslink database",
            if narrowed {
                format!("the proteins of the linear search's peptides at {cutoff} FDR (LinearProteins.fasta)")
            } else {
                "the whole FASTA".into()
            },
        );
    }
    rest(list, xl);
}

fn timstof(list: &mut List, mut b: Map<String, Value>) {
    list.heading("Bruker timsTOF");
    if let Some(v) = take(&mut b, "ion_mobility_scale") {
        list.item("ion mobility", text(&v));
    }
    if let Some(Value::Object(mut ms2)) = take(&mut b, "ms2") {
        if let Some(v) = take(&mut ms2, "spectrum_processing_params") {
            list.item("MS2 processing", text(&v));
        }
        if let Some(v) = take(&mut ms2, "frame_splitting_params") {
            list.item("MS2 frame splitting", text(&v));
        }
        for (k, v) in ms2 {
            list.item(&format!("MS2 {}", label(&k)), text(&v));
        }
    }
    if let Some(Value::Object(mut ms1)) = take(&mut b, "ms1") {
        let mut parts = Vec::new();
        if let Some(v) = take(&mut ms1, "mz_ppm") {
            parts.push(format!("m/z within {} ppm", num(&v)));
        }
        if let Some(v) = take(&mut ms1, "ims_pct") {
            parts.push(format!("ion mobility within {}%", num(&v)));
        }
        parts.extend(ms1.iter().map(|(k, v)| format!("{} {}", label(k), text(v))));
        list.item("MS1 peak merging", parts.join(", "));
    }
    rest(list, b);
}

fn libraries(list: &mut List, mut lib: Map<String, Value>) {
    list.heading("Spectral libraries");
    if let Some(v) = take(&mut lib, "formats") {
        list.item("formats", text(&v));
    }
    let mut fdr = Vec::new();
    for (key, name) in [("psm_fdr", "PSMs"), ("peptide_fdr", "peptides"), ("protein_fdr", "proteins"), ("crosslink_fdr", "crosslinks")] {
        if let Some(v) = take(&mut lib, key) {
            fdr.push(format!("{name} {}", percent(&v)));
        }
    }
    if !fdr.is_empty() {
        list.item("FDR", fdr.join(", "));
    }
    if let (Some(lo), Some(hi)) = (take(&mut lib, "min_fragments"), take(&mut lib, "max_fragments")) {
        list.item("fragments", format!("{} to {} per precursor", num(&lo), num(&hi)));
    }
    for (key, name) in [("proteotypic", "proteotypic only"), ("consensus", "consensus spectra"), ("irt", "iRT"), ("crosslink_decoys", "crosslink decoys")] {
        if let Some(v) = take(&mut lib, key) {
            list.item(name, yes(Some(v)));
        }
    }
    rest(list, lib);
}

/// A file list: "Spectrum files: 21, in <folder>", then one name per line
fn files(out: &mut Vec<String>, title: &str, urls: &[String], sizes: bool) {
    if urls.is_empty() {
        return;
    }
    let paths: Vec<std::path::PathBuf> = urls.iter().map(|u| std::path::PathBuf::from(path(u))).collect();
    let folder = paths[0].parent().map(|p| p.to_path_buf());
    let same = folder.is_some() && paths.iter().all(|p| p.parent() == folder.as_deref());
    match (&folder, same) {
        (Some(f), true) => out.push(format!("{title}: {}, in {}", paths.len(), f.display())),
        _ => out.push(format!("{title}: {}", paths.len())),
    }
    let names: Vec<String> = paths
        .iter()
        .map(|p| match (same, p.file_name()) {
            (true, Some(n)) => n.to_string_lossy().into_owned(),
            _ => p.display().to_string(),
        })
        .collect();
    let width = names.iter().map(|n| n.chars().count()).max().unwrap_or(0);
    let digits = names.len().to_string().len();
    for (i, (name, p)) in names.iter().zip(&paths).enumerate() {
        let size = if sizes { std::fs::metadata(p).ok().filter(|m| m.is_file()).map(|m| bytes(m.len())) } else { None };
        let line = match size {
            Some(size) => format!("  {:>digits$}  {name:<width$}  {size:>9}", i + 1),
            None => format!("  {:>digits$}  {name}", i + 1),
        };
        out.push(line);
    }
}

fn bytes(n: u64) -> String {
    let n = n as f64;
    if n >= 1e9 {
        format!("{:.2} GB", n / 1e9)
    } else if n >= 1e6 {
        format!("{:.1} MB", n / 1e6)
    } else if n >= 1e3 {
        format!("{:.1} kB", n / 1e3)
    } else {
        format!("{n} B")
    }
}

/// A table: the first column left-aligned, the others right-aligned, a rule
/// under the head and, with `total`, above the last row
fn table(out: &mut Vec<String>, head: &[String], rows: &[Vec<String>], total: bool) {
    let mut width: Vec<usize> = head.iter().map(|h| h.chars().count()).collect();
    for row in rows {
        for (w, cell) in width.iter_mut().zip(row) {
            *w = (*w).max(cell.chars().count());
        }
    }
    let line = |cells: &[String]| {
        let mut s = String::from("  ");
        for (i, (cell, w)) in cells.iter().zip(&width).enumerate() {
            if i == 0 {
                s.push_str(&format!("{cell:<w$}"));
            } else {
                s.push_str(&format!("  {cell:>w$}"));
            }
        }
        s.trim_end().to_string()
    };
    let rule = format!("  {}", "-".repeat(width.iter().sum::<usize>() + 2 * width.len().saturating_sub(1)));
    out.push(line(head));
    out.push(rule.clone());
    for (i, row) in rows.iter().enumerate() {
        if total && i + 1 == rows.len() && rows.len() > 1 {
            out.push(rule.clone());
        }
        out.push(line(row));
    }
}

fn count(row: &Map<String, Value>, key: &str) -> u64 {
    row.get(key).and_then(Value::as_u64).unwrap_or(0)
}

fn pct1(row: &Map<String, Value>, key: &str) -> String {
    format!("{:.1}", row.get(key).and_then(Value::as_f64).unwrap_or(0.0))
}

/// The quality-control tables
fn checks(out: &mut Vec<String>, mut qc: Map<String, Value>) {
    if let Some(Value::Array(rows)) = take(&mut qc, "digestion") {
        let rows: Vec<Map<String, Value>> = rows.into_iter().filter_map(|r| r.as_object().cloned()).collect();
        if !rows.is_empty() {
            out.push(
                "Digestion of the peptides at 1% FDR (peptides: targets less decoys; percentages of the peptides)".into(),
            );
            let any = |key: &str| rows.iter().any(|r| count(r, key) > 0);
            let semi = any("semi_n") || any("semi_c") || any("non_enzymatic");
            let mut head: Vec<String> = ["file", "targets", "decoys", "peptides", "0 missed", "1 missed", "2+ missed", "missed %"]
                .map(String::from)
                .to_vec();
            if semi {
                head.extend(["semi N", "semi C", "non-specific", "semi N %", "semi C %", "non-specific %"].map(String::from));
            }
            let body: Vec<Vec<String>> = rows
                .iter()
                .map(|r| {
                    let mut cells = vec![r.get("file").map(text).unwrap_or_default()];
                    for key in ["target_peptides", "decoy_peptides", "peptides", "missed_cleavages_0", "missed_cleavages_1", "missed_cleavages_2_plus"] {
                        cells.push(count(r, key).to_string());
                    }
                    cells.push(pct1(r, "missed_cleavage_pct"));
                    if semi {
                        for key in ["semi_n", "semi_c", "non_enzymatic"] {
                            cells.push(count(r, key).to_string());
                        }
                        for key in ["semi_n_pct", "semi_c_pct", "non_enzymatic_pct"] {
                            cells.push(pct1(r, key));
                        }
                    }
                    cells
                })
                .collect();
            let total = rows.last().and_then(|r| r.get("file")).and_then(Value::as_str) == Some("total");
            table(out, &head, &body, total);
        }
    }
    if let Some(Value::Object(files)) = take(&mut qc, "polymers") {
        if let Some(Value::Object(first)) = files.values().next() {
            let names: Vec<String> = first
                .get("polymers")
                .and_then(Value::as_array)
                .map(|a| a.iter().filter_map(|p| p.get("name").map(text)).collect())
                .unwrap_or_default();
            out.push("MS1 polymer signal (percent of the MS1 ion current)".into());
            let mut head = vec!["file".to_string(), "MS1 spectra".to_string()];
            head.extend(names.iter().map(|n| format!("{n} %")));
            let body: Vec<Vec<String>> = files
                .iter()
                .map(|(file, stats)| {
                    let mut cells = vec![file.clone(), stats.get("ms1_spectra").map(num).unwrap_or_default()];
                    let shares = stats.get("polymers").and_then(Value::as_array).cloned().unwrap_or_default();
                    for n in &names {
                        let pct = shares
                            .iter()
                            .find(|p| p.get("name").map(text).as_deref() == Some(n.as_str()))
                            .and_then(|p| p.get("tic_pct"))
                            .and_then(Value::as_f64);
                        cells.push(pct.map(|x| format!("{x:.2}")).unwrap_or_default());
                    }
                    cells
                })
                .collect();
            table(out, &head, &body, false);
        }
    }
    if let Some(Value::Array(rows)) = take(&mut qc, "diagnostic_ions") {
        if !rows.is_empty() {
            out.push("Diagnostic ions (MS2 spectra with the ion)".into());
            let head = ["file", "ion", "spectra", "%"].map(String::from).to_vec();
            let body: Vec<Vec<String>> = rows
                .iter()
                .filter_map(Value::as_object)
                .map(|r| {
                    vec![
                        r.get("file").map(text).unwrap_or_default(),
                        r.get("ion").map(text).unwrap_or_default(),
                        r.get("spectra").map(num).unwrap_or_default(),
                        r.get("pct").and_then(Value::as_f64).map(|x| format!("{x:.2}")).unwrap_or_default(),
                    ]
                })
                .collect();
            table(out, &head, &body, false);
        }
    }
    if let Some(Value::Array(rows)) = take(&mut qc, "entrapment") {
        if !rows.is_empty() {
            out.push("Entrapment: false discoveries at the reported FDR, from the entrapment proteins".into());
            let head = ["level", "FDR", "targets", "entrapments", "size ratio", "at least false", "estimated false"]
                .map(String::from)
                .to_vec();
            let share = |r: &Map<String, Value>, k: &str| {
                r.get(k).and_then(Value::as_f64).map(|x| format!("{:.2}%", x * 100.0)).unwrap_or_default()
            };
            let body: Vec<Vec<String>> = rows
                .iter()
                .filter_map(Value::as_object)
                .map(|r| {
                    vec![
                        r.get("level").map(text).unwrap_or_default(),
                        r.get("fdr").map(percent).unwrap_or_default(),
                        r.get("targets").map(num).unwrap_or_default(),
                        r.get("entrapments").map(num).unwrap_or_default(),
                        r.get("ratio").and_then(Value::as_f64).map(|x| format!("{x:.2}")).unwrap_or_default(),
                        share(r, "lower_bound"),
                        share(r, "combined"),
                    ]
                })
                .collect();
            table(out, &head, &body, false);
        }
    }
    if !qc.is_empty() {
        out.push("Other checks".into());
        for (k, v) in qc {
            out.push(format!("  {:<LABEL$} {}", label(&k), text(&v)).trim_end().to_string());
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use serde_json::json;

    /// The end of the user's DSSO search of RNA polymerase II, as results.json has it
    fn polii() -> Value {
        serde_json::from_str(
            r#"{
            "version": "1.9.4",
            "database": {
                "bucket_size": 8192,
                "enzyme": {"missed_cleavages": 2, "min_len": 5, "max_len": 50, "cleave_at": null, "restrict": null,
                           "c_terminal": null, "semi_enzymatic": false, "enzymes": ["trypsin"]},
                "peptide_min_mass": 350.0, "peptide_max_mass": 5000.0,
                "ion_kinds": ["b", "y"], "min_ion_index": 2,
                "static_mods": {"C": {"mass": 57.021465, "name": "UNIMOD:4"}},
                "variable_mods": {"^Q": [{"mass": -17.026548, "max_count": 2, "name": "UNIMOD:28"}],
                                  "M": [{"mass": 15.994915, "max_count": 2, "name": "UNIMOD:35"}]},
                "max_variable_mods": 3, "max_combinations": null, "decoy_tag": "rev_", "generate_decoys": true,
                "fasta": "D:\\AaronI\\XNova\\PolII and PRP40_FixedHeader.fasta",
                "prefilter_chunk_size": 0, "prefilter": false, "prefilter_low_memory": true, "clip_n_term_met": false
            },
            "quant": {"tmt": null, "tmt_settings": {"level": 3, "sn": false}, "lfq": false,
                      "lfq_settings": {"peak_scoring": "Hybrid", "integration": "Sum", "spectral_angle": 0.7}},
            "precursor_tol": {"ppm": [-15.0, 15.0]},
            "fragment_tol": {"ppm": [-20.0, 20.0]},
            "precursor_charge": [2, 6], "override_precursor_charge": false, "isotope_errors": [0, 0],
            "deisotope": true, "chimera": false, "wide_window": false, "min_peaks": 15, "max_peaks": 300,
            "max_fragment_charge": null, "min_matched_peaks": 4, "report_psms": 1, "predict_rt": true,
            "mzml_paths": ["file:///D:/AaronI/XNova/Bruker%20Data/20250305_AI_Rui_Xlink_F1_Slot1-10_1_9607.d",
                           "file:///D:/AaronI/XNova/Bruker%20Data/20250305_AI_Rui_Xlink_Unfrac_Slot1-30_1_9631.d"],
            "output_paths": ["file:///D:/AaronI/XNova/XNova_BrukerResults194/results.sage.tsv",
                             "file:///D:/AaronI/XNova/XNova_BrukerResults194/crosslinks.sage.tsv"],
            "bruker_config": {"ms2": {"spectrum_processing_params": {"smoothing_window": 1, "centroiding_window": 1,
                                                                     "calibration_tolerance": 0.1, "calibrate": false},
                                      "frame_splitting_params": {"Quadrupole": {"Even": 1}}},
                              "ms1": {"mz_ppm": 5.0, "ims_pct": 3.0}, "ion_mobility_scale": "calibrated"},
            "protein_grouping": true, "protein_grouping_peptide_fdr": 0.01,
            "crosslink": {"reagent": "dsso", "mass": null, "cleavable": null, "residues": "KSTY", "protein_nterm": true,
                          "allow_cterm_residue": false, "loop_links": true, "multivalent": 2, "candidates": 30,
                          "min_chain_peaks": 3, "require_complete_doublets": false, "precursor_window": true,
                          "backend": "cpu", "chain_survivors": 20000, "linear_search": true, "linear_proteins": false,
                          "linear_protein_fdr": 0.01, "min_linked_peaks": 2},
            "score_type": "SageHyperScore",
            "spectral_library": {"formats": ["openswath", "spectronaut"], "psm_fdr": 0.01, "peptide_fdr": 0.01,
                                 "protein_fdr": 0.01, "crosslink_fdr": 0.01, "proteotypic": true, "consensus": true,
                                 "min_fragments": 3, "max_fragments": 20, "irt": false, "crosslink_decoys": true},
            "qc": {"digestion": [
                {"file": "20250305_AI_Rui_Xlink_F1_Slot1-10_1_9607.d", "target_peptides": 0, "decoy_peptides": 0,
                 "peptides": 0, "missed_cleavages_0": 0, "missed_cleavages_1": 0, "missed_cleavages_2_plus": 0,
                 "semi_n": 0, "semi_c": 0, "non_enzymatic": 0, "missed_cleavage_pct": 0.0, "semi_n_pct": 0.0,
                 "semi_c_pct": 0.0, "non_enzymatic_pct": 0.0},
                {"file": "total", "target_peptides": 844, "decoy_peptides": 9, "peptides": 835,
                 "missed_cleavages_0": 289, "missed_cleavages_1": 331, "missed_cleavages_2_plus": 215,
                 "semi_n": 0, "semi_c": 0, "non_enzymatic": 0, "missed_cleavage_pct": 65.38922155688623,
                 "semi_n_pct": 0.0, "semi_c_pct": 0.0, "non_enzymatic_pct": 0.0}]}
        }"#,
        )
        .unwrap()
    }

    fn find<'a>(lines: &'a [String], start: &str) -> &'a str {
        lines
            .iter()
            .find(|l| l.trim_start().starts_with(start))
            .unwrap_or_else(|| panic!("no line starting {start:?} in\n{}", lines.join("\n")))
    }

    fn value(line: &str) -> &str {
        line[4 + LABEL..].trim()
    }

    #[test]
    fn settings_read_as_words() {
        let lines = lines(&polii());
        assert!(lines[0].starts_with("=== summary"));
        assert_eq!(value(find(&lines, "XNova version")), "1.9.4");
        assert_eq!(value(find(&lines, "enzyme")), "trypsin, up to 2 missed cleavages, 5 to 50 residues");
        assert_eq!(value(find(&lines, "peptide mass")), "350 to 5000 Da");
        assert_eq!(value(find(&lines, "fragment ions")), "b, y, from ion 2");
        assert_eq!(value(find(&lines, "fixed modifications")), "C +57.021465 (UNIMOD:4)");
        // one variable modification per line
        let i = lines.iter().position(|l| l.trim_start().starts_with("variable modifications")).unwrap();
        assert_eq!(value(&lines[i]), "M +15.994915 (UNIMOD:35), up to 2 per peptide");
        assert_eq!(value(&lines[i + 1]), "peptide N-term Q -17.026548 (UNIMOD:28), up to 2 per peptide");
        assert_eq!(value(find(&lines, "decoys")), "reversed proteins, tag rev_");
        assert_eq!(value(find(&lines, "precursor tolerance")), "-15 to 15 ppm");
        assert_eq!(value(find(&lines, "precursor charge")), "2 to 6, where a spectrum has no charge");
        assert_eq!(value(find(&lines, "isotope errors")), "none");
        assert_eq!(value(find(&lines, "protein groups")), "yes, from peptides at 1% FDR");
        assert_eq!(value(find(&lines, "TMT")), "off");
        assert_eq!(value(find(&lines, "LFQ")), "off");
        assert_eq!(value(find(&lines, "reagent")), "DSSO, 158.0038 Da, MS-cleavable");
        assert_eq!(value(find(&lines, "stubs")), "54.0106 and 103.9932 Da, also 85.9826");
        assert_eq!(value(find(&lines, "reacts with")), "KSTY and protein N-termini");
        assert_eq!(value(find(&lines, "multivalent")), "up to 2 linkers");
        assert_eq!(value(find(&lines, "ion mobility")), "calibrated");
        assert_eq!(value(find(&lines, "MS2 frame splitting")), "Quadrupole (Even 1)");
        assert_eq!(value(find(&lines, "MS1 peak merging")), "m/z within 5 ppm, ion mobility within 3%");
        assert_eq!(value(find(&lines, "FDR")), "PSMs 1%, peptides 1%, proteins 1%, crosslinks 1%");
        // no JSON left
        assert!(lines.iter().all(|l| !l.contains(['{', '}', '"'])), "{}", lines.join("\n"));
        // nothing the Progress tab reads as a stage of the search
        assert!(lines.iter().all(|l| !l.contains("model") && !l.contains("discriminant")));
    }

    #[test]
    fn files_are_listed_by_name_under_their_folder() {
        let lines = lines(&polii());
        let i = lines.iter().position(|l| l.starts_with("Spectrum files: 2")).unwrap();
        assert!(lines[i].contains("Bruker Data"), "{}", lines[i]);
        assert_eq!(lines[i + 1], "  1  20250305_AI_Rui_Xlink_F1_Slot1-10_1_9607.d");
        assert_eq!(lines[i + 2], "  2  20250305_AI_Rui_Xlink_Unfrac_Slot1-30_1_9631.d");
        let o = lines.iter().position(|l| l.starts_with("Output files: 2")).unwrap();
        assert_eq!(lines[o + 1], "  1  results.sage.tsv");
    }

    #[test]
    fn digestion_is_a_table_with_its_total_last() {
        let lines = lines(&polii());
        let i = lines.iter().position(|l| l.starts_with("Digestion of the peptides")).unwrap();
        let head = &lines[i + 1];
        assert!(head.starts_with("  file"), "{head}");
        assert!(head.ends_with("missed %"), "semi-enzymatic columns only when there are any: {head}");
        assert!(lines[i + 2].trim_start().starts_with("---"));
        assert!(lines[i + 3].starts_with("  20250305_AI_Rui_Xlink_F1"));
        assert!(lines[i + 4].trim_start().starts_with("---"), "a rule above the total");
        let total: Vec<&str> = lines[i + 5].split_whitespace().collect();
        assert_eq!(total, ["total", "844", "9", "835", "289", "331", "215", "65.4"]);
        // the columns line up
        let end = |l: &str| l.len();
        assert_eq!(end(head), end(&lines[i + 3]));
        assert_eq!(end(head), end(&lines[i + 5]));
    }

    #[test]
    fn unknown_settings_are_still_listed() {
        let mut v = polii();
        v["max_memory_gb"] = json!(32.0);
        v["entrapment"] = json!({"tag": "ENTRAP_"});
        v["database"]["new_setting"] = json!([1, 2]);
        let lines = lines(&v);
        assert_eq!(value(find(&lines, "max memory gb")), "32");
        assert_eq!(value(find(&lines, "entrapment")), "tag ENTRAP_");
        assert_eq!(value(find(&lines, "new setting")), "1, 2");
    }

    #[test]
    fn two_site_linkers_name_each_end() {
        assert_eq!(site_text("KSTY["), "KSTY and protein N-termini");
        assert_eq!(site_text("DE]"), "DE and protein C-termini");
        assert_eq!(site_text(sage_core::crosslink::SDA_CARBENE), "any residue");
        let mut v = polii();
        v["crosslink"] = json!({"reagent": "sda", "site1": "KSTY[", "site2": sage_core::crosslink::SDA_CARBENE});
        let lines = lines(&v);
        assert_eq!(value(find(&lines, "reagent")), "SDA, 82.0419 Da, non-cleavable");
        assert_eq!(value(find(&lines, "reacts with")), "site 1: KSTY and protein N-termini; site 2: any residue");
    }

    #[test]
    fn the_crosslink_database_is_named() {
        let mut v = polii();
        assert_eq!(value(find(&lines(&v), "crosslink database")), "the whole FASTA");
        v["crosslink"]["linear_proteins"] = json!(true);
        v["crosslink"]["linear_protein_fdr"] = json!(0.01);
        assert_eq!(
            value(find(&lines(&v), "crosslink database")),
            "the proteins of the linear search's peptides at 1% FDR (LinearProteins.fasta)"
        );
        assert_eq!(percent_text(0.005), "0.5%");
    }

    #[test]
    fn bruker_settings_only_for_bruker_files() {
        let mut v = polii();
        v["mzml_paths"] = json!(["file:///data/run1.mzML"]);
        let lines = lines(&v);
        assert!(lines.iter().all(|l| !l.contains("timsTOF") && !l.contains("ion mobility")));
    }
}

#[cfg(test)]
mod preview {
    /// Print the summary of a results.json: `XNOVA_SUMMARY_JSON=<file>`
    #[test]
    #[ignore]
    fn print_summary() {
        let path = std::env::var("XNOVA_SUMMARY_JSON").expect("set XNOVA_SUMMARY_JSON");
        let value: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
        for line in super::lines(&value) {
            println!("{line}");
        }
    }
}
