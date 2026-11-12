//! What is known about the proteins of a search, looked up by accession:
//! UniProt (name, gene, organism, function, subcellular location, GO terms,
//! ComplexPortal complexes, Reactome pathways and keywords) and STRING (the
//! known and predicted interactions between them, with the evidence of each).
//!
//! Looked up only when the user asks: the accessions (nothing else) go to
//! rest.uniprot.org and string-db.org. What comes back is kept in the results
//! folder (`network/protein_info.json`), so that a network can be drawn again
//! offline, and is the same however often it is drawn.

use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::io::Read;
use std::path::{Path, PathBuf};
use std::time::Duration;

/// The file the lookups are kept in, in the network folder of the results
pub const FILE: &str = "protein_info.json";
const UNIPROT: &str = "https://rest.uniprot.org/uniprotkb/accessions";
const UNIPROT_FIELDS: &str = "accession,id,reviewed,gene_primary,gene_oln,protein_name,organism_name,organism_id,\
                              cc_function,cc_subcellular_location,go,xref_complexportal,xref_reactome,keyword";
const STRING: &str = "https://string-db.org/api/json";
/// Who is asking, as STRING asks of its users (`caller_identity`)
const CALLER: &str = "XNova";
/// Accessions per UniProt request
const UNIPROT_CHUNK: usize = 100;
/// Proteins per STRING request (STRING's own limit is 2000)
const STRING_MAX: usize = 2000;
/// The lowest STRING score asked for: every interaction STRING has
/// (its "low confidence"); the network's threshold is applied when drawing
const STRING_FLOOR: u32 = 150;

/// A term of an ontology or database: `GO:0005665`, `CPX-2662`, `R-SCE-674695`
/// (locations and keywords: the name, no id)
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(default)]
pub struct Term {
    pub id: String,
    pub name: String,
}

impl Term {
    fn new(id: &str, name: &str) -> Term {
        Term { id: id.to_string(), name: name.to_string() }
    }

    /// What identifies it: the id, else the name
    pub fn key(&self) -> &str {
        if self.id.is_empty() {
            &self.name
        } else {
            &self.id
        }
    }
}

/// A protein's UniProt entry, as far as the network needs it
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct ProteinInfo {
    /// The UniProt accession looked up (`P04050`)
    pub accession: String,
    /// `RPB1_YEAST`
    pub entry: String,
    /// Swiss-Prot (reviewed) rather than TrEMBL
    pub reviewed: bool,
    pub gene: String,
    pub name: String,
    pub organism: String,
    pub taxon: u32,
    /// UniProt's "Function" text
    pub function: String,
    /// UniProt's subcellular locations (`Nucleus`, `Cytoplasm`)
    pub locations: Vec<String>,
    pub go_process: Vec<Term>,
    pub go_function: Vec<Term>,
    pub go_component: Vec<Term>,
    /// ComplexPortal's complexes
    pub complexes: Vec<Term>,
    /// Reactome's pathways
    pub pathways: Vec<Term>,
    /// UniProt's keywords, but for its technical ones (`3D-structure`)
    pub keywords: Vec<Term>,
    /// STRING's id (`4932.YDL140C`), once looked up
    pub string_id: String,
}

/// An interaction STRING knows between two proteins of the results (their
/// accessions as in the results, `a` < `b`), with its scores (0 to 1): the
/// combined score and the evidence channels
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct Interaction {
    pub a: String,
    pub b: String,
    pub score: f32,
    pub experiments: f32,
    pub databases: f32,
    pub textmining: f32,
    pub coexpression: f32,
}

/// Everything looked up for a results folder
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct Annotations {
    /// By accession as in the results (`sp|P04050|RPB1_YEAST`)
    pub proteins: BTreeMap<String, ProteinInfo>,
    /// Accessions UniProt did not know, or that are not UniProt accessions
    pub not_found: BTreeSet<String>,
    pub interactions: Vec<Interaction>,
    /// Accessions whose interactions with each other STRING was asked for
    pub string_queried: BTreeSet<String>,
    /// `2026_04`
    pub uniprot_release: String,
    /// `12.5`
    pub string_version: String,
    /// When last looked up (`2026-11-12`)
    pub looked_up: String,
}

impl Annotations {
    /// Whether anything was looked up yet
    pub fn is_empty(&self) -> bool {
        self.proteins.is_empty() && self.interactions.is_empty() && self.not_found.is_empty()
    }

    /// The interaction STRING has between `a` and `b`; `None` when STRING
    /// was not asked about both (`Some(None)`: asked, and it has none)
    pub fn interaction(&self, a: &str, b: &str) -> Option<Option<&Interaction>> {
        if !self.string_queried.contains(a) || !self.string_queried.contains(b) {
            return None;
        }
        let (a, b) = if a <= b { (a, b) } else { (b, a) };
        Some(self.interactions.iter().find(|i| i.a == a && i.b == b))
    }

    /// The accessions of `wanted` not looked up yet (in UniProt or STRING)
    pub fn missing<'a>(&self, wanted: &'a [String]) -> Vec<&'a String> {
        wanted
            .iter()
            .filter(|a| {
                let uniprot = self.proteins.contains_key(*a) || self.not_found.contains(*a);
                let string = self.string_queried.contains(*a) || !self.proteins.contains_key(*a);
                !uniprot || !string
            })
            .collect()
    }

    /// The saved lookups of a results folder, if any
    pub fn load(folder: &Path) -> Option<Annotations> {
        let text = std::fs::read_to_string(path(folder)).ok()?;
        serde_json::from_str(&text).ok()
    }

    /// Keep them in the results folder; returns the file
    pub fn save(&self, folder: &Path) -> Result<PathBuf, String> {
        let file = path(folder);
        if let Some(dir) = file.parent() {
            std::fs::create_dir_all(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
        }
        let text = serde_json::to_string_pretty(self).map_err(|e| e.to_string())?;
        std::fs::write(&file, text).map_err(|e| format!("{}: {e}", file.display()))?;
        Ok(file)
    }
}

/// Where a results folder's lookups are kept
pub fn path(folder: &Path) -> PathBuf {
    folder.join(super::network::FOLDER).join(FILE)
}

/// Whether `a` is a UniProt accession (`P04050`, `A0A023GPI8`)
fn is_uniprot(a: &str) -> bool {
    let b = a.as_bytes();
    let digit = |c: u8| c.is_ascii_digit();
    let upper = |c: u8| c.is_ascii_uppercase();
    let alnum = |c: u8| digit(c) || upper(c);
    let opq = |c: u8| matches!(c, b'O' | b'P' | b'Q');
    match b.len() {
        6 if opq(b[0]) => digit(b[1]) && alnum(b[2]) && alnum(b[3]) && alnum(b[4]) && digit(b[5]),
        6 | 10 if upper(b[0]) && !opq(b[0]) => {
            digit(b[1]) && b[2..].chunks(4).all(|c| upper(c[0]) && alnum(c[1]) && alnum(c[2]) && digit(c[3]))
        }
        _ => false,
    }
}

/// The UniProt accession of an accession as in a FASTA file:
/// `sp|P04050|RPB1_YEAST` and `tr|…|…` (UniProt's headers), `P04050` and
/// isoforms (`P04050-2`: the entry's); `None` for anything else
pub fn uniprot_accession(accession: &str) -> Option<String> {
    let mut parts = accession.split('|');
    let first = parts.next().unwrap_or("");
    let candidate = match (first, parts.next()) {
        ("sp" | "tr", Some(middle)) => middle,
        _ => first,
    };
    let candidate = candidate.split('-').next().unwrap_or("").trim();
    is_uniprot(candidate).then(|| candidate.to_string())
}

/// Today's date (UTC), `2026-11-12`
pub fn today() -> String {
    let days = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_secs() / 86_400).unwrap_or(0) as i64;
    // civil from days (H. Hinnant)
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe + era * 400 + i64::from(m <= 2);
    format!("{y:04}-{m:02}-{d:02}")
}

// ---------------------------------------------------------------- reading the answers

fn text_at<'a>(v: &'a Value, path: &[&str]) -> &'a str {
    let mut v = v;
    for key in path {
        match v.get(*key) {
            Some(next) => v = next,
            None => return "",
        }
    }
    v.as_str().unwrap_or("")
}

/// The entries of a UniProt answer (`/uniprotkb/accessions`, JSON)
pub fn parse_uniprot(json: &str) -> Result<Vec<ProteinInfo>, String> {
    let v: Value = serde_json::from_str(json).map_err(|e| format!("UniProt's answer is not JSON: {e}"))?;
    let Some(results) = v.get("results").and_then(Value::as_array) else {
        let messages = v.get("messages").map(|m| m.to_string()).unwrap_or_default();
        return Err(format!("UniProt answered without entries {messages}"));
    };
    let mut out = Vec::new();
    for r in results {
        let mut p = ProteinInfo {
            accession: text_at(r, &["primaryAccession"]).to_string(),
            entry: text_at(r, &["uniProtkbId"]).to_string(),
            reviewed: text_at(r, &["entryType"]).contains("Swiss-Prot"),
            organism: text_at(r, &["organism", "scientificName"]).to_string(),
            taxon: r.pointer("/organism/taxonId").and_then(Value::as_u64).unwrap_or(0) as u32,
            ..Default::default()
        };
        p.name = match text_at(r, &["proteinDescription", "recommendedName", "fullName", "value"]) {
            "" => r
                .pointer("/proteinDescription/submissionNames/0/fullName/value")
                .and_then(Value::as_str)
                .unwrap_or("")
                .to_string(),
            name => name.to_string(),
        };
        if let Some(gene) = r.pointer("/genes/0") {
            p.gene = match text_at(gene, &["geneName", "value"]) {
                "" => gene.pointer("/orderedLocusNames/0/value").and_then(Value::as_str).unwrap_or("").to_string(),
                name => name.to_string(),
            };
        }
        for c in r.get("comments").and_then(Value::as_array).into_iter().flatten() {
            match text_at(c, &["commentType"]) {
                "FUNCTION" => {
                    for t in c.get("texts").and_then(Value::as_array).into_iter().flatten() {
                        let text = text_at(t, &["value"]).trim();
                        if !text.is_empty() {
                            if !p.function.is_empty() {
                                p.function.push(' ');
                            }
                            p.function.push_str(text);
                        }
                    }
                }
                "SUBCELLULAR LOCATION" => {
                    for l in c.get("subcellularLocations").and_then(Value::as_array).into_iter().flatten() {
                        let location = text_at(l, &["location", "value"]).trim().to_string();
                        if !location.is_empty() && !p.locations.contains(&location) {
                            p.locations.push(location);
                        }
                    }
                }
                _ => {}
            }
        }
        for x in r.get("uniProtKBCrossReferences").and_then(Value::as_array).into_iter().flatten() {
            let id = text_at(x, &["id"]);
            let property = |key: &str| {
                x.get("properties")
                    .and_then(Value::as_array)
                    .into_iter()
                    .flatten()
                    .find(|p| text_at(p, &["key"]) == key)
                    .map(|p| text_at(p, &["value"]).to_string())
                    .unwrap_or_default()
            };
            match text_at(x, &["database"]) {
                "GO" => {
                    let term = property("GoTerm");
                    let (aspect, name) = term.split_once(':').unwrap_or(("", term.as_str()));
                    let list = match aspect {
                        "P" => &mut p.go_process,
                        "F" => &mut p.go_function,
                        "C" => &mut p.go_component,
                        _ => continue,
                    };
                    list.push(Term::new(id, name));
                }
                "ComplexPortal" => p.complexes.push(Term::new(id, &property("EntryName"))),
                "Reactome" => p.pathways.push(Term::new(id, &property("PathwayName"))),
                _ => {}
            }
        }
        for k in r.get("keywords").and_then(Value::as_array).into_iter().flatten() {
            let name = text_at(k, &["name"]);
            if !name.is_empty() && text_at(k, &["category"]) != "Technical term" {
                p.keywords.push(Term::new(text_at(k, &["id"]), name));
            }
        }
        if !p.accession.is_empty() {
            out.push(p);
        }
    }
    Ok(out)
}

/// STRING's ids (`get_string_ids`): `(query, string id, taxon)`
pub fn parse_string_ids(json: &str) -> Result<Vec<(String, String, u32)>, String> {
    let v: Value = serde_json::from_str(json).map_err(|e| format!("STRING's answer is not JSON: {e}"))?;
    let list = v.as_array().ok_or_else(|| format!("STRING answered: {}", short(json)))?;
    Ok(list
        .iter()
        .filter_map(|x| {
            let query = text_at(x, &["queryItem"]);
            let id = text_at(x, &["stringId"]);
            let taxon = x.get("ncbiTaxonId").and_then(|t| t.as_u64().or_else(|| t.as_str()?.parse().ok())).unwrap_or(0) as u32;
            (!query.is_empty() && !id.is_empty()).then(|| (query.to_string(), id.to_string(), taxon))
        })
        .collect())
}

/// STRING's interactions (`network`): `(string id a, string id b, scores)`
/// with the scores as in [`Interaction`] (`a`, `b` empty)
pub fn parse_string_network(json: &str) -> Result<Vec<(String, String, Interaction)>, String> {
    let v: Value = serde_json::from_str(json).map_err(|e| format!("STRING's answer is not JSON: {e}"))?;
    let list = v.as_array().ok_or_else(|| format!("STRING answered: {}", short(json)))?;
    let score = |x: &Value, key: &str| x.get(key).and_then(Value::as_f64).unwrap_or(0.0) as f32;
    Ok(list
        .iter()
        .filter_map(|x| {
            let (a, b) = (text_at(x, &["stringId_A"]), text_at(x, &["stringId_B"]));
            (!a.is_empty() && !b.is_empty()).then(|| {
                let scores = Interaction {
                    score: score(x, "score"),
                    experiments: score(x, "escore"),
                    databases: score(x, "dscore"),
                    textmining: score(x, "tscore"),
                    coexpression: score(x, "ascore"),
                    ..Default::default()
                };
                (a.to_string(), b.to_string(), scores)
            })
        })
        .collect())
}

fn short(text: &str) -> String {
    text.chars().take(200).collect()
}

// ---------------------------------------------------------------- asking

fn agent() -> ureq::Agent {
    ureq::AgentBuilder::new()
        .timeout(Duration::from_secs(120))
        .try_proxy_from_env(true)
        .user_agent(concat!("XNova/", env!("CARGO_PKG_VERSION"), " (sage-gui)"))
        .build()
}

/// An answer's text (larger than ureq's `into_string` takes)
fn body(response: ureq::Response) -> Result<String, String> {
    let mut bytes = Vec::new();
    response.into_reader().take(512 << 20).read_to_end(&mut bytes).map_err(|e| e.to_string())?;
    text_of(bytes)
}

/// Bytes as text, gzip undone (UniProt gzips its answers twice when asked
/// for gzip, and ureq undoes only the outer layer)
fn text_of(mut bytes: Vec<u8>) -> Result<String, String> {
    for _ in 0..3 {
        if !bytes.starts_with(&[0x1f, 0x8b]) {
            break;
        }
        let mut plain = Vec::new();
        flate2::read::MultiGzDecoder::new(bytes.as_slice()).take(512 << 20).read_to_end(&mut plain).map_err(|e| format!("a gzipped answer that does not unpack: {e}"))?;
        bytes = plain;
    }
    String::from_utf8(bytes).map_err(|_| "an answer that is not text (UTF-8)".to_string())
}

fn failed(service: &str, e: ureq::Error) -> String {
    match e {
        ureq::Error::Status(code, response) => {
            let text = response.into_string().unwrap_or_default();
            format!("{service} answered {code}: {}", short(&text))
        }
        e => format!("Could not reach {service}: {e}"),
    }
}

/// Look `accessions` (as in the results) up in UniProt and STRING, adding to
/// `have` (all of them again if `again`). `progress` hears what is being
/// done. Returns what is known now and the problems met (what could be
/// looked up is kept when the rest could not).
pub fn look_up(accessions: &[String], have: &Annotations, again: bool, progress: &dyn Fn(String)) -> (Annotations, Vec<String>) {
    let mut ann = have.clone();
    let mut problems = Vec::new();
    let agent = agent();
    let wanted: Vec<&String> = if again { accessions.iter().collect() } else { ann.missing(accessions) };

    // UniProt: the entries not looked up yet
    let mut by_uniprot: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for a in accessions {
        match uniprot_accession(a) {
            Some(u) => by_uniprot.entry(u).or_default().push(a.clone()),
            None => {
                ann.not_found.insert(a.clone());
            }
        }
    }
    let ask: Vec<String> = by_uniprot
        .iter()
        .filter(|(_, own)| own.iter().any(|a| wanted.contains(&a) && (again || !ann.proteins.contains_key(a))))
        .map(|(u, _)| u.clone())
        .collect();
    for (i, chunk) in ask.chunks(UNIPROT_CHUNK).enumerate() {
        progress(format!(
            "UniProt: proteins {} to {} of {}…",
            i * UNIPROT_CHUNK + 1,
            (i * UNIPROT_CHUNK + chunk.len()),
            ask.len()
        ));
        let answer = agent
            .get(UNIPROT)
            .query("accessions", &chunk.join(","))
            .query("fields", UNIPROT_FIELDS)
            .query("format", "json")
            .call()
            .map_err(|e| failed("UniProt", e))
            .and_then(|r| {
                let release = r.header("X-UniProt-Release").unwrap_or("").to_string();
                body(r).map(|text| (release, text))
            })
            .and_then(|(release, text)| parse_uniprot(&text).map(|entries| (release, entries)));
        match answer {
            Ok((release, entries)) => {
                if !release.is_empty() {
                    ann.uniprot_release = release;
                }
                let found: BTreeSet<&str> = entries.iter().map(|e| e.accession.as_str()).collect();
                for u in chunk {
                    if !found.contains(u.as_str()) {
                        for a in &by_uniprot[u] {
                            ann.proteins.remove(a);
                            ann.not_found.insert(a.clone());
                        }
                    }
                }
                for entry in entries {
                    for a in by_uniprot.get(&entry.accession).into_iter().flatten() {
                        ann.not_found.remove(a);
                        // a STRING id found before stays
                        let string_id = ann.proteins.get(a).map(|p| p.string_id.clone()).unwrap_or_default();
                        ann.proteins.insert(a.clone(), ProteinInfo { string_id, ..entry.clone() });
                    }
                }
            }
            Err(e) => {
                problems.push(e);
                break;
            }
        }
    }

    // STRING: the interactions between all of them, once any is new
    let in_uniprot: Vec<&String> = accessions.iter().filter(|a| ann.proteins.contains_key(*a)).collect();
    let new = again || in_uniprot.iter().any(|a| !ann.string_queried.contains(*a));
    if new && in_uniprot.len() >= 2 {
        match string_network(&agent, &mut ann, &in_uniprot, progress) {
            Ok(()) => {}
            Err(e) => problems.push(e),
        }
    }
    ann.looked_up = today();
    (ann, problems)
}

/// STRING's ids for the proteins, then the interactions between them, species by species
fn string_network(agent: &ureq::Agent, ann: &mut Annotations, accessions: &[&String], progress: &dyn Fn(String)) -> Result<(), String> {
    let mut by_uniprot: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for a in accessions {
        by_uniprot.entry(ann.proteins[*a].accession.clone()).or_default().push((*a).clone());
    }
    let uniprot: Vec<&String> = by_uniprot.keys().collect();
    if uniprot.len() > STRING_MAX {
        return Err(format!(
            "STRING takes at most {STRING_MAX} proteins at a time; the network has {}: its interactions were not looked up",
            uniprot.len()
        ));
    }
    progress(format!("STRING: finding the {} proteins…", uniprot.len()));
    let ids = agent
        .post(&format!("{STRING}/get_string_ids"))
        .send_form(&[
            ("identifiers", &uniprot.iter().map(|s| s.as_str()).collect::<Vec<_>>().join("\r")),
            ("echo_query", "1"),
            ("limit", "1"),
            ("caller_identity", CALLER),
        ])
        .map_err(|e| failed("STRING", e))
        .and_then(body)
        .and_then(|text| parse_string_ids(&text))?;
    // STRING's id → the results' accessions; species → STRING's ids
    let mut accessions_of: HashMap<String, Vec<String>> = HashMap::new();
    let mut species: BTreeMap<u32, BTreeSet<String>> = BTreeMap::new();
    for (query, id, taxon) in ids {
        for a in by_uniprot.get(&query).into_iter().flatten() {
            if let Some(p) = ann.proteins.get_mut(a) {
                p.string_id = id.clone();
            }
            accessions_of.entry(id.clone()).or_default().push(a.clone());
        }
        species.entry(taxon).or_default().insert(id);
    }
    let mut found = Vec::new();
    for (taxon, ids) in species.iter().filter(|(_, ids)| ids.len() >= 2) {
        // STRING asks for a second between requests
        std::thread::sleep(Duration::from_millis(1100));
        progress(format!("STRING: interactions between {} proteins (taxon {taxon})…", ids.len()));
        let edges = agent
            .post(&format!("{STRING}/network"))
            .send_form(&[
                ("identifiers", &ids.iter().map(String::as_str).collect::<Vec<_>>().join("\r")),
                ("species", &taxon.to_string()),
                ("required_score", &STRING_FLOOR.to_string()),
                ("caller_identity", CALLER),
            ])
            .map_err(|e| failed("STRING", e))
            .and_then(body)
            .and_then(|text| parse_string_network(&text))?;
        found.extend(edges);
    }
    std::thread::sleep(Duration::from_millis(1100));
    if let Ok(version) = agent.get(&format!("{STRING}/version")).call().map_err(|e| failed("STRING", e)).and_then(body) {
        if let Ok(v) = serde_json::from_str::<Value>(&version) {
            let v = text_at(&v[0], &["string_version"]);
            if !v.is_empty() {
                ann.string_version = v.to_string();
            }
        }
    }
    // the interactions between the proteins asked about, replaced
    let asked: BTreeSet<&str> = accessions.iter().map(|a| a.as_str()).collect();
    ann.interactions.retain(|i| !(asked.contains(i.a.as_str()) && asked.contains(i.b.as_str())));
    let mut merged: BTreeMap<(String, String), Interaction> = BTreeMap::new();
    for (id_a, id_b, scores) in found {
        for a in accessions_of.get(&id_a).into_iter().flatten() {
            for b in accessions_of.get(&id_b).into_iter().flatten() {
                if a == b {
                    continue;
                }
                let (a, b) = if a <= b { (a.clone(), b.clone()) } else { (b.clone(), a.clone()) };
                let entry = merged.entry((a.clone(), b.clone())).or_insert(Interaction { a, b, ..Default::default() });
                if scores.score > entry.score {
                    *entry = Interaction { a: entry.a.clone(), b: entry.b.clone(), ..scores.clone() };
                }
            }
        }
    }
    ann.interactions.extend(merged.into_values());
    ann.interactions.sort_by(|x, y| (&x.a, &x.b).cmp(&(&y.a, &y.b)));
    ann.string_queried.extend(accessions.iter().map(|a| (*a).clone()));
    Ok(())
}

#[cfg(test)]
pub(crate) mod test {
    use super::*;

    /// Two entries as UniProt gives them (shortened)
    pub(crate) const UNIPROT_JSON: &str = r#"{"results":[
      {"entryType":"UniProtKB reviewed (Swiss-Prot)","primaryAccession":"P04050","uniProtkbId":"RPB1_YEAST",
       "organism":{"scientificName":"Saccharomyces cerevisiae (strain ATCC 204508 / S288c)","taxonId":559292},
       "proteinDescription":{"recommendedName":{"fullName":{"value":"DNA-directed RNA polymerase II subunit RPB1"}}},
       "genes":[{"geneName":{"value":"RPO21"},"orderedLocusNames":[{"value":"YDL140C"}]}],
       "comments":[{"commentType":"FUNCTION","texts":[{"value":"DNA-dependent RNA polymerase."}]},
                   {"commentType":"SUBCELLULAR LOCATION","subcellularLocations":[{"location":{"value":"Nucleus","id":"SL-0191"}}]}],
       "keywords":[{"id":"KW-0539","category":"Cellular component","name":"Nucleus"},{"id":"KW-1185","category":"Technical term","name":"Reference proteome"}],
       "uniProtKBCrossReferences":[
         {"database":"ComplexPortal","id":"CPX-2662","properties":[{"key":"EntryName","value":"DNA-directed RNA polymerase II complex"}]},
         {"database":"KEGG","id":"sce:YDL140C","properties":[{"key":"Description","value":"-"}]},
         {"database":"Reactome","id":"R-SCE-674695","properties":[{"key":"PathwayName","value":"RNA Polymerase II Pre-transcription Events"}]},
         {"database":"GO","id":"GO:0005665","properties":[{"key":"GoTerm","value":"C:RNA polymerase II, core complex"},{"key":"GoEvidenceType","value":"IDA:SGD"}]},
         {"database":"GO","id":"GO:0003899","properties":[{"key":"GoTerm","value":"F:DNA-directed RNA polymerase activity"}]},
         {"database":"GO","id":"GO:0006366","properties":[{"key":"GoTerm","value":"P:transcription by RNA polymerase II"}]}]},
      {"entryType":"UniProtKB unreviewed (TrEMBL)","primaryAccession":"A0A023GPI8","uniProtkbId":"A0A023GPI8_CANAX",
       "organism":{"scientificName":"Canavalia lineata","taxonId":28957},
       "proteinDescription":{"submissionNames":[{"fullName":{"value":"Lectin alpha chain"}}]}}
    ]}"#;

    #[test]
    fn accessions() {
        assert_eq!(uniprot_accession("sp|P04050|RPB1_YEAST").as_deref(), Some("P04050"));
        assert_eq!(uniprot_accession("tr|A0A023GPI8|A0A023GPI8_CANAX").as_deref(), Some("A0A023GPI8"));
        assert_eq!(uniprot_accession("P04050").as_deref(), Some("P04050"));
        assert_eq!(uniprot_accession("P04050-2").as_deref(), Some("P04050"));
        assert_eq!(uniprot_accession("Q9Y6K9").as_deref(), Some("Q9Y6K9"));
        for not in ["ENSP00000354587", "gene1", "sp|nothing|X", "PO4050", "", "rev_sp"] {
            assert_eq!(uniprot_accession(not), None, "{not}");
        }
        let date = today();
        assert_eq!(date.len(), 10);
        assert!(date.starts_with("20"), "{date}");
    }

    #[test]
    fn uniprot_entries() {
        let entries = parse_uniprot(UNIPROT_JSON).unwrap();
        assert_eq!(entries.len(), 2);
        let rpb1 = &entries[0];
        assert_eq!((rpb1.accession.as_str(), rpb1.entry.as_str(), rpb1.gene.as_str()), ("P04050", "RPB1_YEAST", "RPO21"));
        assert!(rpb1.reviewed);
        assert_eq!(rpb1.taxon, 559292);
        assert_eq!(rpb1.name, "DNA-directed RNA polymerase II subunit RPB1");
        assert_eq!(rpb1.function, "DNA-dependent RNA polymerase.");
        assert_eq!(rpb1.locations, ["Nucleus"]);
        assert_eq!(rpb1.complexes, [Term::new("CPX-2662", "DNA-directed RNA polymerase II complex")]);
        assert_eq!(rpb1.pathways, [Term::new("R-SCE-674695", "RNA Polymerase II Pre-transcription Events")]);
        assert_eq!(rpb1.go_component, [Term::new("GO:0005665", "RNA polymerase II, core complex")]);
        assert_eq!(rpb1.go_function.len(), 1);
        assert_eq!(rpb1.go_process[0].name, "transcription by RNA polymerase II");
        assert_eq!(rpb1.keywords, [Term::new("KW-0539", "Nucleus")], "no technical keywords");
        // TrEMBL: the submitted name, no gene
        let lectin = &entries[1];
        assert!(!lectin.reviewed);
        assert_eq!(lectin.name, "Lectin alpha chain");
        assert_eq!(lectin.gene, "");
        assert!(parse_uniprot(r#"{"messages":["Invalid fields parameter value 'x'"]}"#).unwrap_err().contains("Invalid fields"));
    }

    #[test]
    fn gzipped_answers() {
        use std::io::Write;
        let gzip = |bytes: &[u8]| {
            let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
            e.write_all(bytes).unwrap();
            e.finish().unwrap()
        };
        assert_eq!(text_of(b"{}".to_vec()).unwrap(), "{}");
        assert_eq!(text_of(gzip(b"{\"results\":[]}")).unwrap(), "{\"results\":[]}");
        assert_eq!(text_of(gzip(&gzip(b"twice"))).unwrap(), "twice");
        assert!(text_of(vec![0xff, 0xfe, 0x00]).is_err());
    }

    #[test]
    fn string_answers() {
        let ids = parse_string_ids(
            r#"[{"queryIndex":0,"queryItem":"P04050","stringId":"4932.YDL140C","ncbiTaxonId":4932,"preferredName":"RPO21"},
                {"queryIndex":1,"queryItem":"P08518","stringId":"4932.YOR151C","ncbiTaxonId":"4932"}]"#,
        )
        .unwrap();
        assert_eq!(ids, [("P04050".into(), "4932.YDL140C".into(), 4932), ("P08518".into(), "4932.YOR151C".into(), 4932)]);
        let edges = parse_string_network(
            r#"[{"stringId_A":"4932.YDL140C","stringId_B":"4932.YOR151C","score":0.999,"escore":0.999,"dscore":0.9,"tscore":0.969,"ascore":0.1}]"#,
        )
        .unwrap();
        assert_eq!(edges.len(), 1);
        assert_eq!(edges[0].2.score, 0.999);
        assert_eq!(edges[0].2.databases, 0.9);
        assert!(parse_string_ids(r#"{"Error":"not found"}"#).is_err());
    }

    #[test]
    fn interactions_and_what_is_missing() {
        let mut ann = Annotations::default();
        let all: Vec<String> = ["sp|A|X", "sp|B|Y", "custom"].map(String::from).to_vec();
        assert_eq!(ann.missing(&all).len(), 3);
        ann.proteins.insert("sp|A|X".into(), ProteinInfo::default());
        ann.proteins.insert("sp|B|Y".into(), ProteinInfo::default());
        ann.not_found.insert("custom".into());
        // in UniProt, not yet asked of STRING
        assert_eq!(ann.missing(&all), [&all[0], &all[1]]);
        assert_eq!(ann.interaction("sp|A|X", "sp|B|Y"), None);
        ann.string_queried.extend(["sp|A|X".to_string(), "sp|B|Y".to_string()]);
        assert!(ann.missing(&all).is_empty());
        assert_eq!(ann.interaction("sp|B|Y", "sp|A|X"), Some(None));
        ann.interactions.push(Interaction { a: "sp|A|X".into(), b: "sp|B|Y".into(), score: 0.9, ..Default::default() });
        assert_eq!(ann.interaction("sp|B|Y", "sp|A|X").unwrap().unwrap().score, 0.9);
        // kept in, and read from, the results folder
        let dir = std::env::temp_dir().join(format!("xnova_protein_info_{}", std::process::id()));
        let file = ann.save(&dir).unwrap();
        assert!(file.ends_with(Path::new("network").join(FILE)));
        assert_eq!(Annotations::load(&dir), Some(ann));
        let _ = std::fs::remove_dir_all(&dir);
    }
}
