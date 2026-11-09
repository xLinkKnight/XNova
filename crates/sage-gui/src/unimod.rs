//! Look up modifications in Unimod (https://www.unimod.org) by accession number.
//!
//! XNova bundles Unimod's XML export (`assets/unimod/unimod.xml`, unchanged,
//! under Unimod's Design Science License; see `assets/unimod/README.md`), and
//! sage-gui compiles it in. Lookups need no internet connection and send no
//! traffic to unimod.org. unimod.org is contacted only when the user asks for
//! an accession that is not in the bundled copy, and then at most once a day.
//! That download is kept next to the GUI state and read first while it is newer
//! than the bundled copy.

use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::mpsc::{channel, Receiver};
use std::sync::OnceLock;
use std::time::{Duration, SystemTime};

pub const XML_URL: &str = "https://www.unimod.org/xml/unimod.xml";

/// Unimod's XML export, as published at [`XML_URL`]
pub const BUNDLED: &str = include_str!("../assets/unimod/unimod.xml");
/// When unimod.org last changed the bundled copy
pub const BUNDLED_DATE: &str = "17 February 2026";
/// The same, as the file's `Last-Modified` time (seconds since 1970)
const BUNDLED_MODIFIED: u64 = 1_771_329_288;
/// unimod.org is asked for its current export at most this often
const DOWNLOAD_INTERVAL: Duration = Duration::from_secs(24 * 3600);

/// Unimod's copyright notice, shown with the data
pub const NOTICE: &str = "Copyright (C) 2002-2006 Unimod; this information may be copied, distributed \
    and/or modified under certain conditions, but it comes WITHOUT ANY WARRANTY; see the accompanying \
    Design Science License for more details";

/// One allowed position of a modification
#[derive(Clone, Debug, PartialEq)]
pub struct Specificity {
    /// Residue letter, or `N-term` / `C-term`
    pub site: String,
    /// `Anywhere`, `Any N-term`, `Any C-term`, `Protein N-term`, `Protein C-term`
    pub position: String,
    pub classification: String,
    /// Unimod hides rare or unverified specificities from its default view
    pub hidden: bool,
    /// The same position as a Sage modification target, if Sage can express it
    pub target: Option<String>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Entry {
    pub accession: u32,
    /// PSI-MS name, e.g. `Oxidation`
    pub title: String,
    pub full_name: String,
    pub mono_mass: f64,
    pub composition: String,
    pub specificities: Vec<Specificity>,
}

impl Entry {
    /// ProForma name written into the Sage configuration
    pub fn proforma(&self) -> String {
        format!("UNIMOD:{}", self.accession)
    }
}

/// Accept `35`, `#35`, `UNIMOD:35`, `unimod 35`
pub fn parse_accession(text: &str) -> Option<u32> {
    let t = text.trim().to_ascii_lowercase();
    let t = t.strip_prefix("unimod").unwrap_or(&t);
    let t = t.trim_start_matches([':', '#', ' ']);
    t.parse().ok().filter(|&n| n > 0)
}

/// Unimod site + position as a Sage target: `K`, `^`, `^Q`, `$`, `[`, `[M`, `]`
fn sage_target(site: &str, position: &str) -> Option<String> {
    let terminus = match position {
        "Anywhere" => "",
        "Any N-term" => "^",
        "Any C-term" => "$",
        "Protein N-term" => "[",
        "Protein C-term" => "]",
        _ => return None,
    };
    match site {
        "N-term" | "C-term" if !terminus.is_empty() => Some(terminus.into()),
        aa if aa.len() == 1 && crate::settings::VALID_AA.contains(aa) => Some(format!("{terminus}{aa}")),
        _ => None,
    }
}

fn parse(xml: &str) -> Result<roxmltree::Document<'_>, String> {
    roxmltree::Document::parse(xml).map_err(|e| format!("Unimod XML is unreadable: {e}"))
}

/// One `<umod:mod>` record
fn entry(node: roxmltree::Node, accession: u32) -> Result<Entry, String> {
    let attr = |n: roxmltree::Node, k: &str| n.attribute(k).unwrap_or("").to_string();
    let delta = node
        .children()
        .find(|n| n.has_tag_name("delta"))
        .ok_or(format!("Unimod record {accession} has no mass"))?;
    let mono_mass = delta
        .attribute("mono_mass")
        .and_then(|m| m.parse().ok())
        .ok_or(format!("Unimod record {accession} has no monoisotopic mass"))?;
    let mut specificities: Vec<Specificity> = node
        .children()
        .filter(|n| n.has_tag_name("specificity"))
        .map(|s| {
            let (site, position) = (attr(s, "site"), attr(s, "position"));
            Specificity {
                target: sage_target(&site, &position),
                site,
                position,
                classification: attr(s, "classification"),
                hidden: s.attribute("hidden") == Some("1"),
            }
        })
        .collect();
    // Shown specificities first, then by site, as on unimod.org
    specificities.sort_by(|a, b| (a.hidden, &a.site, &a.position).cmp(&(b.hidden, &b.site, &b.position)));
    Ok(Entry {
        accession,
        title: attr(node, "title"),
        full_name: attr(node, "full_name"),
        mono_mass,
        composition: attr(delta, "composition"),
        specificities,
    })
}

/// Find one record in a Unimod XML. `Ok(None)` when it is not there.
pub fn find(xml: &str, accession: u32) -> Result<Option<Entry>, String> {
    let doc = parse(xml)?;
    let id = accession.to_string();
    match doc
        .descendants()
        .find(|n| n.has_tag_name("mod") && n.attribute("record_id") == Some(id.as_str()))
    {
        Some(node) => entry(node, accession).map(Some),
        None => Ok(None),
    }
}

/// Every record of a Unimod XML, by accession
pub fn index(xml: &str) -> Result<HashMap<u32, Entry>, String> {
    let doc = parse(xml)?;
    Ok(doc
        .descendants()
        .filter(|n| n.has_tag_name("mod"))
        .filter_map(|n| {
            let accession = n.attribute("record_id")?.parse().ok()?;
            entry(n, accession).ok().map(|e| (accession, e))
        })
        .collect())
}

/// The bundled copy, read once
pub fn bundled() -> &'static HashMap<u32, Entry> {
    static INDEX: OnceLock<HashMap<u32, Entry>> = OnceLock::new();
    INDEX.get_or_init(|| index(BUNDLED).unwrap_or_default())
}

fn cache_path() -> Option<PathBuf> {
    let base = std::env::var_os("APPDATA").or_else(|| std::env::var_os("HOME"))?;
    Some(PathBuf::from(base).join("sage-gui").join("unimod.xml"))
}

/// When the downloaded copy was saved, if there is one
fn downloaded_at() -> Option<SystemTime> {
    std::fs::metadata(cache_path()?).ok()?.modified().ok()
}

fn newer_than_bundled(t: SystemTime) -> bool {
    t > SystemTime::UNIX_EPOCH + Duration::from_secs(BUNDLED_MODIFIED)
}

/// Look up an accession on this computer: in a downloaded copy while it is
/// newer than the bundled one, else in the bundled copy (then an older
/// download). `Ok(None)` when no local copy has it. Never uses the network.
pub fn lookup(accession: u32) -> Option<Entry> {
    let downloaded = cache_path().and_then(|p| std::fs::read_to_string(p).ok());
    let newer = downloaded_at().is_some_and(newer_than_bundled);
    lookup_in(downloaded.as_deref(), newer, accession)
}

fn lookup_in(downloaded: Option<&str>, newer: bool, accession: u32) -> Option<Entry> {
    // an unreadable download is passed over: the bundled copy still answers
    let from_download = || downloaded.and_then(|xml| find(xml, accession).ok().flatten());
    if newer {
        from_download().or_else(|| bundled().get(&accession).cloned())
    } else {
        bundled().get(&accession).cloned().or_else(from_download)
    }
}

/// Whether unimod.org may be asked again: no download in the last day
pub fn may_download() -> bool {
    downloaded_at()
        .and_then(|t| t.elapsed().ok())
        .map_or(true, |age| age >= DOWNLOAD_INTERVAL)
}

/// Download Unimod's current export, only when the user asks, and at most
/// once a day. `Ok(false)` when a download from the last day is already kept.
pub fn download() -> Result<bool, String> {
    if !may_download() {
        return Ok(false);
    }
    let xml = ureq::get(XML_URL)
        .timeout(Duration::from_secs(60))
        .call()
        .map_err(|e| format!("Could not download Unimod ({XML_URL}): {e}"))?
        .into_string()
        .map_err(|e| format!("Could not read the Unimod download: {e}"))?;
    if !xml.contains("<umod:unimod") || index(&xml).map_or(true, |records| records.is_empty()) {
        return Err(format!("{XML_URL} did not return the Unimod database"));
    }
    let path = cache_path().ok_or("There is no folder to keep the Unimod download in")?;
    if let Some(dir) = path.parent() {
        let _ = std::fs::create_dir_all(dir);
    }
    std::fs::write(&path, &xml).map_err(|e| format!("Could not save the Unimod download: {e}"))?;
    Ok(true)
}

/// What a lookup found
#[derive(Clone, Debug, PartialEq)]
pub enum Found {
    Entry(Entry),
    /// In no copy on this computer
    Missing {
        accession: u32,
        message: String,
        /// unimod.org may be asked (no download in the last day)
        may_download: bool,
    },
    Error(String),
}

fn missing(accession: u32, downloaded: bool) -> Found {
    let may_download = may_download();
    let message = if downloaded || !may_download {
        format!("Unimod has no record with accession #{accession} (unimod.org was checked within the last day).")
    } else {
        format!(
            "#{accession} is not in the copy of Unimod bundled with XNova (Unimod of {BUNDLED_DATE}). \
             It may be a newer record."
        )
    };
    Found::Missing { accession, message, may_download }
}

/// State of the Unimod panel on the Modifications tab
#[derive(Default)]
pub struct Lookup {
    pub input: String,
    pub pending: Option<Receiver<Found>>,
    /// The pending work is a download from unimod.org
    pub downloading: bool,
    pub result: Option<Found>,
}

impl Lookup {
    /// Look up the typed accession on this computer, on a background thread
    pub fn start(&mut self, ctx: &eframe::egui::Context) {
        let Some(accession) = parse_accession(&self.input) else {
            self.result = Some(Found::Error(format!(
                "`{}` is not a Unimod accession number; enter e.g. 35 or UNIMOD:35",
                self.input.trim()
            )));
            return;
        };
        self.spawn(ctx, false, move || match lookup(accession) {
            Some(entry) => Found::Entry(entry),
            None => missing(accession, false),
        });
    }

    /// The user asked to check unimod.org for an accession no local copy has
    pub fn check_online(&mut self, ctx: &eframe::egui::Context, accession: u32) {
        self.spawn(ctx, true, move || match download() {
            Err(e) => Found::Error(e),
            Ok(downloaded) => match lookup(accession) {
                Some(entry) => Found::Entry(entry),
                None => missing(accession, downloaded),
            },
        });
    }

    fn spawn(&mut self, ctx: &eframe::egui::Context, downloading: bool, work: impl FnOnce() -> Found + Send + 'static) {
        let (tx, rx) = channel();
        let ctx = ctx.clone();
        std::thread::spawn(move || {
            let _ = tx.send(work());
            ctx.request_repaint();
        });
        self.pending = Some(rx);
        self.downloading = downloading;
        self.result = None;
    }

    pub fn busy(&mut self) -> bool {
        if let Some(rx) = &self.pending {
            if let Ok(result) = rx.try_recv() {
                self.result = Some(result);
                self.pending = None;
            }
        }
        self.pending.is_some()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    const XML: &str = r#"<?xml version="1.0" encoding="utf-8"?>
<umod:unimod xmlns:umod="http://www.unimod.org/xmlns/schema/unimod_2">
  <umod:modifications>
    <umod:mod title="Oxidation" full_name="Oxidation or Hydroxylation" record_id="35">
      <umod:specificity hidden="1" site="C" position="Anywhere" classification="Post-translational" spec_group="9"/>
      <umod:specificity hidden="0" site="M" position="Anywhere" classification="Artefact" spec_group="8"/>
      <umod:specificity hidden="1" site="G" position="Any C-term" classification="Pre-translational" spec_group="11"/>
      <umod:delta mono_mass="15.994915" avge_mass="15.9994" composition="O">
        <umod:element symbol="O" number="1"/>
      </umod:delta>
    </umod:mod>
    <umod:mod title="Acetyl" full_name="Acetylation" record_id="1">
      <umod:specificity hidden="0" site="N-term" position="Protein N-term" classification="Post-translational" spec_group="5"/>
      <umod:specificity hidden="0" site="N-term" position="Any N-term" classification="Chemical derivative" spec_group="2"/>
      <umod:delta mono_mass="42.010565" avge_mass="42.0367" composition="H(2) C(2) O"/>
    </umod:mod>
    <umod:mod title="Future" full_name="A record newer than the bundled copy" record_id="99999">
      <umod:specificity hidden="0" site="K" position="Anywhere" classification="Other" spec_group="1"/>
      <umod:delta mono_mass="1.5" avge_mass="1.5" composition="H"/>
    </umod:mod>
  </umod:modifications>
</umod:unimod>"#;

    #[test]
    fn accessions() {
        for (text, n) in [("35", Some(35)), (" #35 ", Some(35)), ("UNIMOD:35", Some(35)), ("unimod 4", Some(4))] {
            assert_eq!(parse_accession(text), n, "{text}");
        }
        for bad in ["", "0", "ox", "UNIMOD:", "-3"] {
            assert_eq!(parse_accession(bad), None, "{bad}");
        }
    }

    #[test]
    fn records_map_to_sage_targets() {
        let ox = find(XML, 35).unwrap().unwrap();
        assert_eq!((ox.title.as_str(), ox.mono_mass, ox.composition.as_str()), ("Oxidation", 15.994915, "O"));
        assert_eq!(ox.proforma(), "UNIMOD:35");
        // shown specificity first
        let targets: Vec<_> = ox.specificities.iter().map(|s| (s.target.clone().unwrap(), s.hidden)).collect();
        assert_eq!(targets, [("M".into(), false), ("C".into(), true), ("$G".into(), true)]);

        let ac = find(XML, 1).unwrap().unwrap();
        let targets: Vec<_> = ac.specificities.iter().filter_map(|s| s.target.clone()).collect();
        assert_eq!(targets, ["^", "["]);

        assert_eq!(find(XML, 999).unwrap(), None);
        // the index agrees with find
        let all = index(XML).unwrap();
        assert_eq!(all.len(), 3);
        assert_eq!(all[&35], ox);
    }

    /// The bundled copy: complete, unchanged, with its notice; no network
    #[test]
    fn bundled_unimod() {
        assert!(BUNDLED.contains(&NOTICE[..60]), "copyright notice kept");
        assert!(BUNDLED.contains("Design Science License"));
        assert!(BUNDLED.trim_end().ends_with("</umod:unimod>"));
        let all = bundled();
        assert!(all.len() > 1500, "{} records", all.len());

        let cam = &all[&4];
        assert_eq!(cam.title, "Carbamidomethyl");
        assert!((cam.mono_mass - 57.021464).abs() < 1e-6);
        assert!(cam.specificities.iter().any(|s| s.target.as_deref() == Some("C") && !s.hidden));

        // common modifications reach the targets the GUI presets use
        let shown = |acc: u32| -> Vec<String> {
            all[&acc].specificities.iter().filter(|s| !s.hidden).filter_map(|s| s.target.clone()).collect()
        };
        for (acc, targets) in [(21, &["S", "T", "Y"][..]), (1, &["[", "K"]), (737, &["K", "^"]), (28, &["^Q"])] {
            let got = shown(acc);
            for t in targets {
                assert!(got.contains(&t.to_string()), "UNIMOD:{acc}: {t} not in {got:?}");
            }
        }
        // the license travels with the data
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("assets/unimod");
        let license = std::fs::read_to_string(dir.join("LICENSE-unimod.txt")).unwrap();
        assert!(license.contains("DESIGN SCIENCE LICENSE"));
    }

    #[test]
    fn newer_downloads_answer_first() {
        let fake = XML.replace("Oxidation or Hydroxylation", "Changed upstream");
        // a newer download answers first, and adds records the bundled copy lacks
        assert_eq!(lookup_in(Some(&fake), true, 35).unwrap().full_name, "Changed upstream");
        assert_eq!(lookup_in(Some(&fake), true, 99999).unwrap().title, "Future");
        // an older download only fills in what the bundled copy lacks
        assert_eq!(lookup_in(Some(&fake), false, 35).unwrap().full_name, "Oxidation or Hydroxylation");
        assert_eq!(lookup_in(Some(&fake), false, 99999).unwrap().title, "Future");
        // no download, or an unreadable one: the bundled copy
        assert_eq!(lookup_in(None, false, 4).unwrap().title, "Carbamidomethyl");
        assert_eq!(lookup_in(Some("not xml"), true, 4).unwrap().title, "Carbamidomethyl");
        assert_eq!(lookup_in(None, false, 99999), None);

        let t = |s| SystemTime::UNIX_EPOCH + Duration::from_secs(s);
        assert!(newer_than_bundled(t(BUNDLED_MODIFIED + 1)));
        assert!(!newer_than_bundled(t(BUNDLED_MODIFIED)));
    }
}
