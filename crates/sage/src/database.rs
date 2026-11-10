use crate::enzyme::{group_digests, DigestGroup, Enzyme, EnzymeParameters};
use crate::fasta::Fasta;
use crate::ion_series::{IonSeries, Kind};
use crate::mass::Tolerance;
use crate::modification::{
    DetailedMod, ModEntry, ModId, ModificationError, ModificationSpecificity, Modifications, SearchMode,
    SiteMode,
};
use crate::peptide::{Peptide, Site};
use crate::scoring::Feature;
use crate::site_rules::{SiteContext, SiteRule};
use std::borrow::Cow;
use dashmap::DashSet;
use fnv::FnvBuildHasher;
use itertools::Itertools;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap};
use std::hash::Hash;
use std::sync::Arc;

#[derive(Deserialize, Serialize, Clone, Debug)]
pub struct EnzymeBuilder {
    /// How many missed cleavages to use
    pub missed_cleavages: Option<u8>,
    /// Minimum peptide length that will be fragmented
    pub min_len: Option<usize>,
    /// Maximum peptide length that will be fragmented
    pub max_len: Option<usize>,
    pub cleave_at: Option<String>,
    pub restrict: Option<String>,
    pub c_terminal: Option<bool>,
    pub semi_enzymatic: Option<bool>,
    /// Several enzymes used together (co-digestion), each a protease name
    /// (`"trypsin"`, `"asp-n"`, see [`crate::enzyme::PROTEASES`]) or a rule
    /// `{"cleave_at": "D", "restrict": "", "c_terminal": false}`. Replaces
    /// `cleave_at` / `restrict` / `c_terminal` when set.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub enzymes: Option<Vec<EnzymeRule>>,
}

/// One enzyme of [`EnzymeBuilder::enzymes`]
#[derive(Deserialize, Serialize, Clone, Debug, PartialEq)]
#[serde(untagged)]
pub enum EnzymeRule {
    Named(String),
    Custom {
        cleave_at: String,
        #[serde(default)]
        restrict: String,
        #[serde(default = "default_c_terminal")]
        c_terminal: bool,
    },
}

fn default_c_terminal() -> bool {
    true
}

impl EnzymeRule {
    /// `(cleave_at, restrict, c_terminal)`
    pub fn resolve(&self) -> Result<(String, String, bool), String> {
        let valid = |s: &str| s.bytes().all(|b| crate::mass::VALID_AA.contains(&b));
        match self {
            EnzymeRule::Named(name) => crate::enzyme::protease(name)
                .map(|(c, r, t)| (c.into(), r.into(), t))
                .ok_or_else(|| {
                    format!(
                        "unknown enzyme `{name}`; use one of {} or a rule {{\"cleave_at\": .., \"restrict\": .., \"c_terminal\": ..}}",
                        crate::enzyme::PROTEASES.join(", ")
                    )
                }),
            EnzymeRule::Custom {
                cleave_at,
                restrict,
                c_terminal,
            } => {
                if cleave_at.is_empty() || !valid(cleave_at) {
                    return Err(format!(
                        "enzyme rule `cleave_at` must be amino acid letters, got `{cleave_at}`"
                    ));
                }
                if !valid(restrict) {
                    return Err(format!(
                        "enzyme rule `restrict` must be amino acid letters, got `{restrict}`"
                    ));
                }
                Ok((cleave_at.clone(), restrict.clone(), *c_terminal))
            }
        }
    }
}

impl EnzymeBuilder {
    /// Check the `enzymes` list (names must be known, rules well-formed)
    pub fn validate(&self) -> Result<(), String> {
        if let Some(rules) = &self.enzymes {
            if rules.is_empty() {
                return Err("`database.enzyme.enzymes` must list at least one enzyme".into());
            }
            for rule in rules {
                rule.resolve()?;
            }
        }
        Ok(())
    }
}

impl Default for EnzymeBuilder {
    fn default() -> Self {
        Self {
            missed_cleavages: Some(0),
            min_len: Some(5),
            max_len: Some(50),
            cleave_at: Some("KR".into()),
            restrict: Some("P".into()),
            c_terminal: Some(true),
            semi_enzymatic: Some(false),
            enzymes: None,
        }
    }
}

impl From<EnzymeBuilder> for EnzymeParameters {
    fn from(en: EnzymeBuilder) -> EnzymeParameters {
        let semi_enzymatic = en.semi_enzymatic.unwrap_or(false);
        let enzyme = match &en.enzymes {
            Some(rules) => Enzyme::combine(
                rules
                    .iter()
                    .map(|rule| {
                        let (cleave, restrict, c_terminal) =
                            rule.resolve().unwrap_or_else(|e| panic!("{e}"));
                        Enzyme::new(&cleave, &restrict, c_terminal, semi_enzymatic)
                            .expect("enzyme rules always cleave")
                    })
                    .collect(),
            ),
            None => Enzyme::new(
                &en.cleave_at.unwrap_or_else(|| "KR".into()),
                &en.restrict.unwrap_or_else(|| "".into()),
                en.c_terminal.unwrap_or(true),
                semi_enzymatic,
            ),
        };
        EnzymeParameters {
            missed_cleavages: en.missed_cleavages.unwrap_or(1),
            min_len: en.min_len.unwrap_or(5),
            max_len: en.max_len.unwrap_or(50),
            enzyme,
        }
    }
}

#[derive(Deserialize, Default)]
/// Parameters used for generating the fragment database
pub struct Builder {
    /// This parameter allows tuning of the internal search structure
    pub bucket_size: Option<usize>,

    pub enzyme: Option<EnzymeBuilder>,
    /// Minimum peptide monoisotopic mass that will be fragmented
    pub peptide_min_mass: Option<f32>,
    /// Maximum peptide monoisotopic mass that will be fragmented
    pub peptide_max_mass: Option<f32>,
    /// Which kind of fragment ions to generate (a, b, c, x, y, z, z_radical)
    pub ion_kinds: Option<Vec<Kind>>,
    /// Minimum ion index to be generated: 1 will remove b1/y1 ions
    /// 2 will remove b1/b2/y1/y2 ions, etc
    pub min_ion_index: Option<usize>,
    /// Static modifications to add to matching amino acids.
    /// Each entry is either a bare mass (`57.0215`) or an object with `mass` and
    /// an optional ProForma `name` (`{"mass": 57.0215, "name": "UNIMOD:4"}`).
    /// Keys: a residue or `^ $ [ ]` symbol, or a named site (see
    /// [`crate::site_rules`]).
    #[serde(default, deserialize_with = "crate::site_rules::named_static_mods")]
    pub static_mods: Option<HashMap<ModificationSpecificity, ModEntry>>,
    /// Variable modifications to add to matching amino acids.
    /// Each entry is either a bare mass (`15.9949`) or an object with `mass` and
    /// optional `max_count` and `name` fields
    /// (`{"mass": 15.9949, "max_count": 1, "name": "UNIMOD:35"}`), and
    /// optionally `position`, `motif` and `neutral_losses`.
    #[serde(default, deserialize_with = "crate::site_rules::named_variable_mods")]
    pub variable_mods: Option<HashMap<ModificationSpecificity, Vec<ModEntry>>>,
    /// Limit number of variable modifications on a peptide
    pub max_variable_mods: Option<usize>,
    /// Hard cap on the total peptide variants generated per input peptide,
    /// including its unmodified form. Values below 1 are normalized to 1.
    /// Variants with fewer PTMs are preferred (generated first).
    pub max_combinations: Option<usize>,
    /// Use this prefix for decoy proteins
    pub decoy_tag: Option<String>,

    pub generate_decoys: Option<bool>,
    /// Path to fasta database
    pub fasta: Option<String>,
    /// Number of sequences to handle simultaneously when pre-filtering the db
    pub prefilter_chunk_size: Option<usize>,
    /// Pre-filter the database to minimize memory usage
    pub prefilter: Option<bool>,
    /// Pre-filter the database with a minimal amount of memory at the cost of speed
    pub prefilter_low_memory: Option<bool>,
    /// Also digest each protein without its initiator Met when methionine
    /// aminopeptidase removes it (second residue G, A, S, T, C, P or V);
    /// those peptides are protein N-terminal (default: false)
    pub clip_n_term_met: Option<bool>,
    /// SILAC: labelled channels; every peptide is also searched in each
    pub silac: Option<crate::silac::SilacSettings>,
    /// Proteins uniformly labelled with 15N and/or 13C
    pub labelled_proteins: Option<crate::labelled::LabelledProteins>,
    /// Two FASTA files, searched as one for linear peptides; crosslinks only
    /// between a protein of one and a protein of the other. Replaces `fasta`.
    pub crosslink_databases: Option<CrosslinkDatabases>,
    /// Known modification sites, for variable modifications with
    /// `site_mode: library` (see [`crate::site_library`])
    pub site_library: Option<crate::site_library::SiteLibrarySettings>,
}

/// Two protein databases of a crosslink search between them, e.g. two
/// complexes: only crosslinks joining a protein of `database_1` to one of
/// `database_2` are searched. Linear peptides are searched in both, as one.
#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct CrosslinkDatabases {
    pub database_1: String,
    pub database_2: String,
}

/// Peptide of a protein in Database 1 / Database 2 (bits of a group)
pub const DATABASE_1: u8 = 1;
pub const DATABASE_2: u8 = 2;

impl Builder {
    pub fn make_parameters(self) -> Result<Parameters, ModificationError> {
        let bucket_size = self.bucket_size.unwrap_or(8192).next_power_of_two();
        let mut static_mods = self.static_mods.unwrap_or_default();
        // chemical labelling (dimethyl): the light channel's labels are static
        if let Some(silac) = &self.silac {
            silac.validate().map_err(ModificationError::Silac)?;
            for (specificity, entry) in silac.light_static_mods() {
                match static_mods.get(&specificity) {
                    // e.g. a configuration written back by a search
                    Some(existing) if (existing.mass() - entry.mass()).abs() < 1e-4 => {}
                    Some(_) => {
                        return Err(ModificationError::Silac(format!(
                            "silac: the light label on {specificity} replaces its static modification; remove that"
                        )))
                    }
                    None => {
                        static_mods.insert(specificity, entry);
                    }
                }
            }
        }
        let mut modifications =
            ResolvedModifications::new(&static_mods, &self.variable_mods.unwrap_or_default())?;
        if modifications.has_library_mods() && self.site_library.is_none() {
            return Err(ModificationError::Rule(
                "a modification with site_mode library needs database.site_library".into(),
            ));
        }
        if modifications.has_library_mods() && !self.generate_decoys.unwrap_or(true) {
            return Err(ModificationError::Rule(
                "site_mode library needs generate_decoys: decoys must carry their targets' library sites".into(),
            ));
        }
        if !modifications.offset_mods.is_empty() && (self.silac.is_some() || self.labelled_proteins.is_some()) {
            return Err(ModificationError::Rule(
                "mass offsets (search_mode mass_offset) cannot be combined with SILAC or labelled proteins".into(),
            ));
        }
        let mut configured: Vec<(ModificationSpecificity, ModId)> = modifications
            .static_mods
            .iter()
            .copied()
            .chain(modifications.variable_mods.iter().map(|&(s, id, _)| (s, id)))
            .collect();
        // uniformly labelled proteins first: SILAC labels go on top of theirs
        let protein_labels = match &self.labelled_proteins {
            Some(labelled) => {
                labelled.validate().map_err(ModificationError::Silac)?;
                let labels = crate::labelled::ProteinLabels::register(labelled, &mut modifications.registry, &configured)?;
                configured.extend(labels.all_ids().into_iter().map(|(r, id)| (ModificationSpecificity::Residue(r), id)));
                Some(Arc::new(labels))
            }
            None => None,
        };
        let labels = match &self.silac {
            Some(silac) => {
                silac.validate().map_err(ModificationError::Silac)?;
                Some(Arc::new(crate::silac::Labels::register(
                    silac,
                    &mut modifications.registry,
                    &configured,
                )?))
            }
            None => None,
        };
        Ok(Parameters {
            bucket_size,
            peptide_min_mass: self.peptide_min_mass.unwrap_or(500.0),
            peptide_max_mass: self.peptide_max_mass.unwrap_or(5000.0),
            ion_kinds: self.ion_kinds.unwrap_or(vec![Kind::B, Kind::Y]),
            min_ion_index: self.min_ion_index.unwrap_or(2),
            decoy_tag: self.decoy_tag.unwrap_or_else(|| "rev_".into()),
            enzyme: self.enzyme.unwrap_or_default(),
            modifications,
            max_variable_mods: self.max_variable_mods.map(|x| x.max(1)).unwrap_or(2),
            max_combinations: self.max_combinations.map(|x| x.max(1)),
            generate_decoys: self.generate_decoys.unwrap_or(true),
            // two databases replace the single FASTA
            fasta: match &self.crosslink_databases {
                Some(_) => String::new(),
                None => self.fasta.expect("A fasta file must be provided!"),
            },
            prefilter_chunk_size: self.prefilter_chunk_size.unwrap_or(0),
            prefilter: self.prefilter.unwrap_or(false),
            prefilter_low_memory: self.prefilter_low_memory.unwrap_or(true),
            clip_n_term_met: self.clip_n_term_met.unwrap_or(false),
            silac: self.silac,
            labels,
            labelled_proteins: self.labelled_proteins,
            protein_labels,
            crosslink_databases: self.crosslink_databases,
            database_groups: None,
            site_library: self.site_library,
            library: None,
        })
    }

    pub fn update_fasta(&mut self, fasta: String) {
        self.fasta = Some(fasta)
    }
}

#[derive(Serialize, Clone, Debug)]
pub struct Parameters {
    pub bucket_size: usize,
    pub enzyme: EnzymeBuilder,
    pub peptide_min_mass: f32,
    pub peptide_max_mass: f32,
    pub ion_kinds: Vec<Kind>,
    pub min_ion_index: usize,
    /// Serialized as the `static_mods` and `variable_mods` it was built from
    #[serde(flatten)]
    pub modifications: ResolvedModifications,
    pub max_variable_mods: usize,
    pub max_combinations: Option<usize>,
    pub decoy_tag: String,
    pub generate_decoys: bool,
    pub fasta: String,
    pub prefilter_chunk_size: usize,
    pub prefilter: bool,
    pub prefilter_low_memory: bool,
    /// Also search proteins without their initiator Met (see [`Builder::clip_n_term_met`])
    pub clip_n_term_met: bool,
    /// SILAC channels (`None`: no labels)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub silac: Option<crate::silac::SilacSettings>,
    /// The SILAC labels as modification ids
    #[serde(skip)]
    pub labels: Option<Arc<crate::silac::Labels>>,
    /// Uniformly 15N / 13C labelled proteins
    #[serde(skip_serializing_if = "Option::is_none")]
    pub labelled_proteins: Option<crate::labelled::LabelledProteins>,
    /// Their labels as modification ids
    #[serde(skip)]
    pub protein_labels: Option<Arc<crate::labelled::ProteinLabels>>,
    /// Two databases for a crosslink search between them (`fasta` is empty)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub crosslink_databases: Option<CrosslinkDatabases>,
    /// With two databases: which of them each protein is in (accession ->
    /// [`DATABASE_1`] and/or [`DATABASE_2`]), set by [`Parameters::merge_databases`]
    #[serde(skip)]
    pub database_groups: Option<Arc<HashMap<Arc<str>, u8>>>,
    /// Known modification sites (`site_mode: library`)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub site_library: Option<crate::site_library::SiteLibrarySettings>,
    /// The site library, resolved to the modifications by
    /// [`Parameters::load_site_library`]
    #[serde(skip)]
    pub library: Option<Arc<crate::site_library::ResolvedLibrary>>,
}

impl Parameters {
    /// Resolve the rows of a site library (`text`, see
    /// [`crate::site_library`]) to the `site_mode: library` modifications
    pub fn load_site_library(&mut self, text: &str) -> Result<&crate::site_library::ResolvedLibrary, String> {
        let max_q = self.site_library.as_ref().map_or(0.01, |s| s.max_q);
        let rows = crate::site_library::parse(text, max_q)?;
        let mut library = crate::site_library::ResolvedLibrary { rows: rows.len(), ..Default::default() };
        let mods = &self.modifications;
        for row in rows {
            let ids: Vec<ModId> = mods
                .variable_mods
                .iter()
                .zip(&mods.variable_rules)
                .filter(|((specificity, id, _), rule)| {
                    rule.library == Some(*id)
                        && *specificity == ModificationSpecificity::Residue(row.residue)
                        && crate::site_library::names(mods.registry.get(*id), &row.modification)
                })
                .map(|((_, id, _), _)| *id)
                .unique()
                .collect();
            if ids.is_empty() {
                continue;
            }
            library.used += 1;
            let entries = library.sites.entry(row.protein).or_default();
            for id in ids {
                entries.push((row.position - 1, row.residue, id));
            }
        }
        for entries in library.sites.values_mut() {
            entries.sort_unstable();
            entries.dedup();
        }
        self.library = Some(Arc::new(library));
        Ok(self.library.as_deref().expect("just set"))
    }

    /// The FASTA files searched: `fasta`, or the two crosslink databases
    pub fn fasta_files(&self) -> Vec<&str> {
        match &self.crosslink_databases {
            Some(d) => vec![d.database_1.as_str(), d.database_2.as_str()],
            None => vec![self.fasta.as_str()],
        }
    }

    /// Two crosslink databases read: one FASTA of both (for the linear
    /// search), remembering which database each protein is in
    pub fn merge_databases(&mut self, database_1: Fasta, database_2: Fasta) -> Fasta {
        let mut groups: HashMap<Arc<str>, u8> = HashMap::new();
        for (accession, _) in &database_1.targets {
            *groups.entry(accession.clone()).or_default() |= DATABASE_1;
        }
        for (accession, _) in &database_2.targets {
            *groups.entry(accession.clone()).or_default() |= DATABASE_2;
        }
        self.database_groups = Some(Arc::new(groups));
        database_1.merge(database_2)
    }

    pub fn auto_calculate_prefilter_chunk_size(&mut self, fasta: &Fasta) {
        const MAX_PEPS_PER_CHUNK: usize = 2usize.pow(23);
        self.prefilter_chunk_size = match self.prefilter_chunk_size {
            0 => {
                let enzyme = self.enzyme.clone().into();
                let total_unmodified_pep_count: usize = fasta.digest_clipping(&enzyme, self.clip_n_term_met).len();
                let combination_factor = if self.max_variable_mods >= usize::BITS as usize {
                    usize::MAX
                } else {
                    1usize << self.max_variable_mods
                };
                let mut mod_count_estimate = (self.modifications.variable_specificities() + 1)
                    .saturating_mul(combination_factor);
                if let Some(max_combinations) = self.max_combinations {
                    mod_count_estimate = mod_count_estimate.min(max_combinations);
                }
                let chunk_count = mod_count_estimate.saturating_mul(total_unmodified_pep_count)
                    / MAX_PEPS_PER_CHUNK;
                fasta
                    .targets
                    .len()
                    .checked_div(chunk_count)
                    .unwrap_or(fasta.targets.len())
            }
            x => x,
        };
    }

    pub fn digest(&self, fasta: &Fasta) -> Vec<Peptide> {
        log::trace!("digesting fasta");
        let enzyme = self.enzyme.clone().into();
        // Generate all tryptic peptide sequences, including reversed (decoy)
        // and missed cleavages, if applicable.
        let digests = fasta.digest_clipping(&enzyme, self.clip_n_term_met);

        log::trace!("grouping digests");
        let start_num = digests.len();
        let digests = group_digests(digests);
        log::trace!(
            "grouped {} digests into {} groups",
            start_num,
            digests.len()
        );

        let targets: DashSet<_, FnvBuildHasher> = DashSet::default();
        digests
            .par_iter()
            .filter(|digest| !digest.reference.decoy)
            .for_each(|digest| {
                targets.insert(digest.reference.sequence.clone().into_bytes());
            });

        log::trace!("modifying peptides");
        let proteins = self.motif_proteins(fasta);
        let mut target_decoys = digests
            .into_par_iter()
            .filter_map(|group| self.with_context(&proteins, group))
            .flat_map_iter(|(peptide, context)| self.forms(peptide, context, &targets))
            .collect::<Vec<_>>();

        Self::reorder_peptides(&mut target_decoys);

        target_decoys
    }

    /// Protein sequences by accession when a motif or the site library
    /// needs them (else empty)
    fn motif_proteins<'f>(&self, fasta: &'f Fasta) -> HashMap<&'f str, &'f [u8]> {
        match self.modifications.has_motifs() || self.library.is_some() {
            true => fasta.targets.iter().map(|(acc, seq)| (acc.as_ref(), seq.as_bytes())).collect(),
            false => HashMap::new(),
        }
    }

    /// A digested peptide, and where it sits in its (first) protein when
    /// motifs need to look past it; with a site library, the library's
    /// sites on the peptide in each of its proteins
    fn with_context(
        &self,
        proteins: &HashMap<&str, &[u8]>,
        group: DigestGroup,
    ) -> Option<(Peptide, Option<SiteContext>)> {
        let sequence = group.reference.sequence.as_bytes();
        let mut context = proteins
            .get(group.reference.protein.as_ref())
            .and_then(|protein| SiteContext::find(protein, sequence));
        if let Some(library) = &self.library {
            let mut sites: Vec<(usize, ModId)> = Vec::new();
            for accession in &group.proteins {
                let (Some(listed), Some(protein)) =
                    (library.sites.get(accession.as_ref()), proteins.get(accession.as_ref()))
                else {
                    continue;
                };
                let Some(start) = protein.windows(sequence.len().max(1)).position(|w| w == sequence) else {
                    continue;
                };
                for &(position, residue, id) in listed {
                    if (start..start + sequence.len()).contains(&position) && sequence[position - start] == residue {
                        sites.push((position - start, id));
                    }
                }
            }
            if !sites.is_empty() {
                sites.sort_unstable();
                sites.dedup();
                context
                    .get_or_insert_with(|| SiteContext { sequence: sequence.to_vec(), offset: 0, library: Vec::new() })
                    .library = sites;
            }
        }
        Peptide::try_from(group).ok().map(|peptide| (peptide, context))
    }

    /// Every form of one digested peptide the search indexes: its
    /// modifications, labelled states and SILAC channels within the mass
    /// range, each with its decoy (unless the decoy is a target sequence)
    fn forms<'a>(
        &'a self,
        peptide: Peptide,
        context: Option<crate::site_rules::SiteContext>,
        targets: &'a DashSet<Vec<u8>, FnvBuildHasher>,
    ) -> impl Iterator<Item = Peptide> + 'a {
        let mods = &self.modifications;
        peptide
            .apply_with_rules(
                &mods.variable_mods,
                &mods.variable_rules,
                &mods.static_mods,
                &mods.static_rules,
                &mods.registry,
                self.max_variable_mods,
                self.max_combinations,
                context.as_ref(),
            )
            .into_iter()
            // uniformly labelled proteins: their peptides labelled
            .flat_map(move |peptide| match &self.protein_labels {
                Some(labels) => labels.states(peptide),
                None => vec![peptide],
            })
            // SILAC: each peptide light and in every labelled channel
            .flat_map(move |peptide| match &self.labels {
                Some(labels) => labels.states(peptide),
                None => vec![peptide],
            })
            .filter(move |peptide| {
                peptide.monoisotopic >= self.peptide_min_mass && peptide.monoisotopic <= self.peptide_max_mass
            })
            .flat_map(move |peptide| {
                if self.generate_decoys {
                    vec![peptide.reverse(), peptide].into_iter()
                } else {
                    vec![peptide].into_iter()
                }
            })
            .filter(move |peptide| !peptide.decoy || !targets.contains(&(peptide.sequence[..])))
    }

    /// The size of the database a search would build, counted without
    /// building it: each digested peptide's forms are generated and counted
    /// one at a time (the same forms the search indexes), so only one
    /// peptide's forms are held per thread. Peptides found in several
    /// proteins are merged by the search, so the counts are an upper bound.
    pub fn estimate(&self, fasta: &Fasta) -> DatabaseEstimate {
        let enzyme = self.enzyme.clone().into();
        let digests = fasta.digest_clipping(&enzyme, self.clip_n_term_met);
        if digests.is_empty() {
            return DatabaseEstimate { proteins: fasta.targets.len(), ..Default::default() };
        }
        let digests = group_digests(digests);
        let targets: DashSet<_, FnvBuildHasher> = DashSet::default();
        digests.par_iter().filter(|digest| !digest.reference.decoy).for_each(|digest| {
            targets.insert(digest.reference.sequence.clone().into_bytes());
        });
        let per_peptide = std::mem::size_of::<Peptide>() as u64;
        let kinds = self.ion_kinds.len();
        let sum = |a: [u64; 4], b: [u64; 4]| [a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]];
        let proteins = self.motif_proteins(fasta);
        let [peptides, decoys, fragments, bytes] = digests
            .into_par_iter()
            .filter_map(|group| self.with_context(&proteins, group))
            .map(|(peptide, context)| {
                self.forms(peptide, context, &targets).fold([0u64; 4], |acc, p| {
                    let len = p.sequence.len();
                    let fragments = (kinds * len.saturating_sub(1 + self.min_ion_index)) as u64;
                    // the peptide, its sequence, modification masses and protein list; 8 bytes per fragment
                    let bytes = per_peptide + 5 * len as u64 + 16 + 8 * fragments;
                    sum(acc, [1, p.decoy as u64, fragments, bytes])
                })
            })
            .reduce(|| [0; 4], sum);
        DatabaseEstimate {
            proteins: fasta.targets.len(),
            peptides: peptides as usize,
            decoys: decoys as usize,
            fragments: fragments as usize,
            bytes,
        }
    }

    pub fn reorder_peptides(target_decoys: &mut Vec<Peptide>) {
        log::trace!("sorting and deduplicating peptides");

        let init_size = target_decoys.len();
        // This is equivalent to a stable sort
        target_decoys.par_sort_unstable_by(|a, b| {
            a.monoisotopic
                .total_cmp(&b.monoisotopic)
                .then_with(|| a.initial_sort(b))
        });
        target_decoys.dedup_by(|remove, keep| {
            if remove.monoisotopic == keep.monoisotopic
                && remove.sequence == keep.sequence
                && remove.modifications == keep.modifications
                && remove.nterm == keep.nterm
                && remove.cterm == keep.cterm
            {
                keep.proteins.extend(remove.proteins.iter().cloned());
                // When merging peptides from different Fastas,
                // decoys in one fasta might be targets in another
                keep.decoy &= remove.decoy;
                true
            } else {
                false
            }
        });

        target_decoys
            .par_iter_mut()
            .for_each(|peptide| peptide.proteins.sort_unstable());

        let num_dropped = init_size - target_decoys.len();
        log::trace!(
            "dropped {} t/d pairs, remaining {}",
            num_dropped,
            target_decoys.len(),
        );
    }

    pub fn build(self, fasta: Fasta) -> IndexedDatabase {
        let target_decoys = self.digest(&fasta);
        self.build_from_peptides(target_decoys)
    }

    pub fn build_from_peptides(self, target_decoys: Vec<Peptide>) -> IndexedDatabase {
        log::trace!("generating fragments");

        // Finally, perform in silico digest for our target sequences
        // Note that multiple charge states are actually handled by
        // [`SpectrumProcessor`] or during scoring - all theoretical
        // fragments are monoisotopic/uncharged
        let mut fragments = target_decoys
            .par_iter()
            .enumerate()
            .flat_map_iter(|(idx, peptide)| {
                // Generate both B and Y ions, then filter down to make sure that
                // theoretical fragments are within the search space
                self.ion_kinds
                    .iter()
                    .flat_map(|kind| {
                        IonSeries::new(peptide, *kind, &self.modifications.registry).enumerate()
                    })
                    .filter(|(ion_idx, ion)| {
                        // Don't store b1, b2, y1, y2 ions for preliminary scoring

                        match ion.kind {
                            Kind::A | Kind::B | Kind::C => (ion_idx + 1) > self.min_ion_index,
                            Kind::X | Kind::Y | Kind::Z | Kind::ZRadical => {
                                peptide.sequence.len().saturating_sub(1) - ion_idx
                                    > self.min_ion_index
                            }
                        }
                    })
                    .map(move |(_, ion)| Theoretical {
                        peptide_index: PeptideIx(idx as u32),
                        fragment_mz: ion.monoisotopic_mass,
                    })
            })
            .collect::<Vec<_>>();
        log::trace!("finalizing index");

        // Sort all of our theoretical fragments by m/z, from low to high
        fragments.par_sort_unstable_by(|a, b| a.fragment_mz.total_cmp(&b.fragment_mz));

        // Now, we bucket all of our theoretical fragments, and within each bucket
        // sort by precursor m/z - and save the minimum *fragment* m/z in a separate
        // vector so that we can perform an efficient binary search to reduce
        // the number of in silico fragments we evaluate
        //
        // Imagine our theoretical fragments look like this
        //
        // Fragment        A      B       C       D       E       F       G       H
        // Fragment m/z [ 1.0    1.2     1.3     2.5     2.5     2.6     3.5     4.0 ]
        // Parent m/z   [ 500    439     291     800     142     515     517     232 ]
        //
        // If we apply a bucket size of 4 we will end up with the following:
        //
        // Fragment        C      B       A       D       E       H       F       G
        // Fragment m/z [ 1.3    1.2     1.0     2.5     2.5     4.0     3.5     2.6 ]
        // Parent m/z   [ 291    439     500     800     142     232     515     517 ]
        //              |___________________________|   |____________________________|
        //               Bucket 1: min m/z 1.0          Bucket 2: min m/z 2.5
        //
        // * Example query: Fragment m/z 1.3 - 1.9 & Precursor m/z: 450 - 900
        // 1) Perform a binary search to narrow down our window to Bucket 1 only
        //      * Bucket 2 has a min m/z outside of our query range - nothing here can match
        //
        // Fragment        C      B       A       D
        // Fragment m/z [ 1.3    1.2     1.0     2.5
        // Parent m/z   [ 291    439     500     800
        //                            |_____________|
        //                                    ^
        //                                    |
        // Window with matching precursors ___|

        // and within Bucket 1, we can perform another binary search to find fragments
        // matching our desired precursor m/z tolerance

        let min_value = fragments
            .par_chunks_mut(self.bucket_size)
            .map(|chunk| {
                // There should always be at least one item in the chunk!
                //  we know the chunk is already sorted by fragment_mz too, so this is minimum value
                let min = chunk[0].fragment_mz;
                chunk.par_sort_unstable_by(|a, b| a.peptide_index.cmp(&b.peptide_index));
                min
            })
            .collect::<Vec<_>>();

        // SILAC: the channel of every peptide, for pairing crosslinks and quantifying
        let silac_channels = match &self.labels {
            Some(labels) => target_decoys
                .par_iter()
                .map(|p| labels.channel(p).unwrap_or(NO_CHANNEL))
                .collect(),
            None => Vec::new(),
        };

        // two databases: which of them each peptide's proteins are in
        let database_groups = match &self.database_groups {
            Some(groups) => target_decoys
                .par_iter()
                .map(|p| p.proteins.iter().fold(0, |g, acc| g | groups.get(acc).copied().unwrap_or(0)))
                .collect(),
            None => Vec::new(),
        };

        // mass offsets: one per modification, with all of its sites
        let mut mass_offsets: Vec<MassOffset> = Vec::new();
        for (specificity, id, rule) in &self.modifications.offset_mods {
            match mass_offsets.iter_mut().find(|o| o.id == *id) {
                Some(offset) => offset.sites.push((*specificity, rule.clone())),
                None => mass_offsets.push(MassOffset {
                    id: *id,
                    mass: self.modifications.registry.get(*id).mass,
                    sites: vec![(*specificity, rule.clone())],
                }),
            }
        }

        IndexedDatabase {
            database_groups,
            mass_offsets,
            offset_peptides: Vec::new(),
            peptides: target_decoys,
            fragments,
            min_value,
            bucket_size: self.bucket_size,
            ion_kinds: self.ion_kinds,
            generate_decoys: self.generate_decoys,
            modifications: self.modifications.registry,
            decoy_tag: self.decoy_tag,
            silac: self.labels,
            silac_channels,
            protein_labels: self.protein_labels,
        }
    }
}

/// Configured modifications, as ids into `registry`
#[derive(Clone, Debug, Default)]
pub struct ResolvedModifications {
    pub registry: Modifications,
    pub static_mods: Vec<(ModificationSpecificity, ModId)>,
    /// `(specificity, modification, per_mod_limit)`, in a stable order
    pub variable_mods: Vec<(ModificationSpecificity, ModId, Option<usize>)>,
    /// Where each static / variable modification may sit (same order)
    pub static_rules: Vec<crate::site_rules::SiteRule>,
    pub variable_rules: Vec<crate::site_rules::SiteRule>,
    /// Variable modifications searched as mass offsets
    /// (`search_mode: mass_offset`), not indexed
    pub offset_mods: Vec<(ModificationSpecificity, ModId, SiteRule)>,
}

/// Check a configured modification's rule and neutral losses, and register
/// the losses
fn register_rule(
    registry: &mut Modifications,
    specificity: ModificationSpecificity,
    id: ModId,
    entry: &ModEntry,
) -> Result<crate::site_rules::SiteRule, ModificationError> {
    use ModificationSpecificity::*;
    let rule = entry.rule();
    let residue = match specificity {
        Residue(r) => Some(r),
        PeptideN(r) | PeptideC(r) | ProteinN(r) | ProteinC(r) => r,
    };
    if !rule.is_any() {
        let Some(r) = residue else {
            return Err(ModificationError::Rule(format!(
                "modification on {specificity}: a position or motif needs a residue"
            )));
        };
        if !rule.position.is_anywhere() && !matches!(specificity, Residue(_)) {
            return Err(ModificationError::Rule(format!(
                "modification on {specificity}: a terminal residue has no other position"
            )));
        }
        if let Some(motif) = &rule.motif {
            if !motif.site_residues().contains(&r) {
                return Err(ModificationError::Rule(format!(
                    "modification on {specificity}: the site of motif {} is never {}",
                    motif.pattern(),
                    r as char
                )));
            }
        }
    }
    for &loss in entry.neutral_losses() {
        if !loss.is_finite() || loss <= 0.0 {
            return Err(ModificationError::Rule(format!(
                "modification on {specificity}: neutral losses must be positive masses"
            )));
        }
        registry.add_neutral_loss(id, residue.unwrap_or(0), loss);
    }
    let mut rule = rule;
    if entry.site_mode() == SiteMode::Library {
        if !matches!(specificity, Residue(_)) {
            return Err(ModificationError::Rule(format!(
                "modification on {specificity}: site_mode library is for residues, not termini"
            )));
        }
        rule.library = Some(id);
    }
    Ok(rule)
}

impl ResolvedModifications {
    /// Does a configured modification have a motif (which needs the protein)?
    pub fn has_motifs(&self) -> bool {
        self.static_rules.iter().chain(&self.variable_rules).any(|r| r.motif.is_some())
    }

    /// Is a variable modification limited to the site library?
    pub fn has_library_mods(&self) -> bool {
        self.variable_rules.iter().any(|r| r.library.is_some())
    }

    /// Number of distinct sites (e.g. `M`, `^`) with variable modifications
    fn variable_specificities(&self) -> usize {
        // `variable_mods` is sorted by specificity
        self.variable_mods
            .iter()
            .map(|(specificity, _, _)| specificity)
            .dedup()
            .count()
    }

    /// Serializable config entry for a registered modification
    fn entry(
        &self,
        specificity: ModificationSpecificity,
        id: ModId,
        max_count: Option<usize>,
        rule: Option<&crate::site_rules::SiteRule>,
        search_mode: SearchMode,
    ) -> ModEntry {
        use ModificationSpecificity::*;
        let modification = self.registry.get(id);
        let residue = match specificity {
            Residue(r) => Some(r),
            PeptideN(r) | PeptideC(r) | ProteinN(r) | ProteinC(r) => r,
        };
        let neutral_losses: Vec<f32> = self.registry.neutral_losses(id, residue.unwrap_or(0)).collect();
        let rule = rule.cloned().unwrap_or_default();
        match (&modification.name, max_count) {
            (None, None) if rule.is_any() && neutral_losses.is_empty() && search_mode.is_indexed() => {
                ModEntry::Mass(modification.mass)
            }
            (name, max_count) => ModEntry::Detailed(DetailedMod {
                mass: modification.mass,
                max_count,
                name: name.clone(),
                position: rule.position,
                site_mode: match rule.library {
                    Some(_) => SiteMode::Library,
                    None => SiteMode::Exhaustive,
                },
                motif: rule.motif,
                neutral_losses,
                search_mode,
            }),
        }
    }

    /// Register the configured modifications. Ids depend only on the
    /// configuration, so every database built from the same parameters (e.g.
    /// each prefilter chunk) agrees on them.
    pub fn new(
        static_mods: &HashMap<ModificationSpecificity, ModEntry>,
        variable_mods: &HashMap<ModificationSpecificity, Vec<ModEntry>>,
    ) -> Result<Self, ModificationError> {
        let mut registry = Modifications::default();

        if let Some((specificity, _)) = static_mods
            .iter()
            .find(|(_, entry)| entry.max_count().is_some())
        {
            return Err(ModificationError::StaticMaxCount(*specificity));
        }
        for (specificity, entry) in static_mods {
            if !entry.search_mode().is_indexed() || !entry.site_mode().is_exhaustive() {
                return Err(ModificationError::Rule(format!(
                    "static modification on {specificity}: search_mode and site_mode are for variable modifications"
                )));
            }
        }

        let mut static_entries = static_mods.iter().collect::<Vec<_>>();
        static_entries.sort_unstable_by_key(|(specificity, _)| **specificity);
        let mut static_rules = Vec::with_capacity(static_entries.len());
        let static_mods = static_entries
            .into_iter()
            .map(|(specificity, entry)| {
                let id = registry.insert(entry.mass(), entry.name())?;
                static_rules.push(register_rule(&mut registry, *specificity, id, entry)?);
                Ok((*specificity, id))
            })
            .collect::<Result<Vec<_>, _>>()?;

        // Flatten variable modifications into a stable order. This matters when
        // `max_combinations` truncates variants: equivalent configurations must
        // retain the same variants regardless of randomized `HashMap` iteration.
        let mut variable_entries = variable_mods
            .iter()
            .flat_map(|(specificity, entries)| {
                entries
                    .iter()
                    .enumerate()
                    .map(|(entry_order, entry)| (*specificity, entry_order, entry))
            })
            .collect::<Vec<_>>();
        variable_entries.sort_unstable_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
        let mut variable_rules = Vec::with_capacity(variable_entries.len());
        let mut offset_mods = Vec::new();
        let mut variable = Vec::with_capacity(variable_entries.len());
        for (specificity, _, entry) in variable_entries {
            let id = registry.insert(entry.mass(), entry.name())?;
            let rule = register_rule(&mut registry, specificity, id, entry)?;
            match entry.search_mode() {
                SearchMode::Indexed => {
                    variable_rules.push(rule);
                    variable.push((specificity, id, entry.max_count()));
                }
                SearchMode::MassOffset => {
                    let why = if entry.mass() == 0.0 || !entry.mass().is_finite() {
                        Some("needs a non-zero mass")
                    } else if rule.motif.is_some() {
                        Some("cannot have a motif (offsets are placed on the peptide alone)")
                    } else if rule.library.is_some() {
                        Some("cannot use site_mode library")
                    } else {
                        None
                    };
                    if let Some(why) = why {
                        return Err(ModificationError::Rule(format!("mass offset on {specificity} {why}")));
                    }
                    offset_mods.push((specificity, id, rule));
                }
            }
        }
        let variable_mods = variable;
        if let Some((specificity, id, _)) =
            offset_mods.iter().find(|(_, id, _)| variable_mods.iter().any(|(_, v, _)| v == id))
        {
            return Err(ModificationError::Rule(format!(
                "modification {} on {specificity} is both indexed and a mass offset; search it one way",
                registry.get(*id)
            )));
        }
        if offset_mods.iter().map(|(_, id, _)| id).unique().count() > MAX_MASS_OFFSETS {
            return Err(ModificationError::Rule(format!("at most {MAX_MASS_OFFSETS} mass offsets")));
        }

        Ok(ResolvedModifications {
            registry,
            static_mods,
            variable_mods,
            static_rules,
            variable_rules,
            offset_mods,
        })
    }
}

impl Serialize for ResolvedModifications {
    /// Same shape as the `static_mods` and `variable_mods` config fields
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        use serde::ser::SerializeStruct;

        let static_mods = self
            .static_mods
            .iter()
            .enumerate()
            .map(|(i, &(specificity, id))| {
                (specificity, self.entry(specificity, id, None, self.static_rules.get(i), SearchMode::Indexed))
            })
            .collect::<BTreeMap<_, _>>();
        let mut variable_mods = BTreeMap::<_, Vec<_>>::new();
        for (i, &(specificity, id, max_count)) in self.variable_mods.iter().enumerate() {
            variable_mods.entry(specificity).or_default().push(self.entry(
                specificity,
                id,
                max_count,
                self.variable_rules.get(i),
                SearchMode::Indexed,
            ));
        }
        for (specificity, id, rule) in &self.offset_mods {
            variable_mods.entry(*specificity).or_default().push(self.entry(
                *specificity,
                *id,
                None,
                Some(rule),
                SearchMode::MassOffset,
            ));
        }

        let mut state = serializer.serialize_struct("ResolvedModifications", 2)?;
        state.serialize_field("static_mods", &static_mods)?;
        state.serialize_field("variable_mods", &variable_mods)?;
        state.end()
    }
}

#[derive(Hash, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize)]
#[repr(transparent)]
pub struct PeptideIx(pub u32);

// This is unsafe for use outside of this crate
impl Default for PeptideIx {
    fn default() -> Self {
        Self(u32::MAX)
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Serialize)]
pub struct Theoretical {
    pub peptide_index: PeptideIx,
    pub fragment_mz: f32,
}

/// The size of a search's database, from [`Parameters::estimate`]
#[derive(Debug, Clone, Default, PartialEq, Serialize)]
pub struct DatabaseEstimate {
    pub proteins: usize,
    /// Indexed peptide forms, targets and decoys
    pub peptides: usize,
    pub decoys: usize,
    /// Theoretical fragments in the index
    pub fragments: usize,
    /// Memory of the peptides and the fragment index (bytes)
    pub bytes: u64,
}

#[derive(Default)]
pub struct IndexedDatabase {
    pub peptides: Vec<Peptide>,
    pub fragments: Vec<Theoretical>,
    pub ion_kinds: Vec<Kind>,
    pub min_value: Vec<f32>,
    pub modifications: Modifications,
    pub bucket_size: usize,
    pub generate_decoys: bool,
    pub decoy_tag: String,
    /// SILAC labels, when searched with them
    pub silac: Option<Arc<crate::silac::Labels>>,
    /// SILAC channel of each peptide (by index), [`NO_CHANNEL`] for one
    /// without labelled residues; empty without SILAC
    pub silac_channels: Vec<u8>,
    /// Uniformly 15N / 13C labelled proteins, when searched with them
    pub protein_labels: Option<Arc<crate::labelled::ProteinLabels>>,
    /// Two crosslink databases: of each peptide (by index), [`DATABASE_1`]
    /// and/or [`DATABASE_2`]; empty with one database
    pub database_groups: Vec<u8>,
    /// Variable modifications searched as mass offsets
    pub mass_offsets: Vec<MassOffset>,
    /// Peptides with a mass offset placed, given an index after the search
    /// by [`IndexedDatabase::materialize_mass_offsets`] (indices from
    /// `peptides.len()` on)
    pub offset_peptides: Vec<Peptide>,
}

/// At most this many mass offsets (a candidate keeps its offset in a byte,
/// 0 for none)
pub const MAX_MASS_OFFSETS: usize = u8::MAX as usize - 1;

/// A variable modification searched as a mass offset (`search_mode:
/// mass_offset`, after Sage Plus, MIT). It is not expanded into the fragment
/// index: every spectrum is also searched at its precursor mass less the
/// offset, matching fragments both as indexed and shifted by the offset, and
/// each candidate is then scored with the offset on every site it may take.
/// At most one offset goes on a peptide, on a residue (or terminus) without
/// another modification, besides the indexed variable modifications.
#[derive(Clone, Debug)]
pub struct MassOffset {
    pub id: ModId,
    pub mass: f32,
    /// The sites it may take, with their position rules
    pub sites: Vec<(ModificationSpecificity, SiteRule)>,
}

/// A mass offset placed on a candidate: [`IndexedDatabase::mass_offsets`]
/// index, and site
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct OffsetPlacement {
    pub offset: u8,
    pub site: Site,
}

/// SILAC channel of a peptide without labelled residues
pub const NO_CHANNEL: u8 = u8::MAX;

impl IndexedDatabase {
    /// The uniform isotope label of a peptide (from a labelled protein)
    pub fn isotope_label(&self, ix: PeptideIx) -> Option<crate::labelled::Isotope> {
        self.protein_labels.as_ref()?.isotope(&self[ix])
    }

    /// Modifications that are only heavy isotopes (SILAC labels, uniform
    /// labels, and SILAC labels on uniformly labelled residues): a residue
    /// carrying one is still free to react, e.g. to crosslink
    pub fn isotope_mods(&self) -> Vec<ModId> {
        let uniform: Vec<ModId> = self
            .protein_labels
            .as_ref()
            .map(|l| l.label_ids().into_iter().map(|(_, id)| id).collect())
            .unwrap_or_default();
        let mut ids = uniform.clone();
        if let Some(silac) = &self.silac {
            ids.extend(silac.label_ids());
            ids.extend(silac.combined_with(&uniform));
        }
        ids.sort();
        ids.dedup();
        ids
    }

    /// SILAC channel of a peptide: `Some(0)` light, `None` without SILAC or
    /// labelled residues
    pub fn silac_channel(&self, ix: PeptideIx) -> Option<u8> {
        self.silac_channels.get(ix.0 as usize).copied().filter(|&c| c != NO_CHANNEL)
    }

    /// Can these peptides be crosslinked into one species? Without SILAC, or
    /// with mixed crosslinks allowed, always; otherwise they share a channel
    /// (peptides without labelled residues fit any)
    pub fn same_silac_channel(&self, peptides: &[PeptideIx]) -> bool {
        if self.silac_channels.is_empty() || self.silac.as_ref().is_some_and(|l| l.mixed_crosslinks) {
            return true;
        }
        let mut channel = None;
        for &p in peptides {
            if let Some(c) = self.silac_channel(p) {
                match channel {
                    Some(x) if x != c => return false,
                    _ => channel = Some(c),
                }
            }
        }
        true
    }

    /// Searched with two crosslink databases?
    pub fn two_databases(&self) -> bool {
        !self.database_groups.is_empty()
    }

    /// Database(s) of a peptide's proteins, bits [`DATABASE_1`] and
    /// [`DATABASE_2`]; `None` with one database
    pub fn database_group(&self, ix: PeptideIx) -> Option<u8> {
        self.database_groups.get(ix.0 as usize).copied()
    }

    /// Can a link join these two peptides as far as the databases go? With
    /// two databases, one peptide must come from Database 1 and the other from
    /// Database 2 (a peptide of proteins in both counts as either)
    pub fn cross_database(&self, a: PeptideIx, b: PeptideIx) -> bool {
        match (self.database_group(a), self.database_group(b)) {
            (Some(x), Some(y)) => {
                (x & DATABASE_1 != 0 && y & DATABASE_2 != 0) || (x & DATABASE_2 != 0 && y & DATABASE_1 != 0)
            }
            _ => true,
        }
    }

    /// May these two peptides be crosslinked? SILAC channels and the two
    /// crosslink databases
    pub fn may_crosslink(&self, a: PeptideIx, b: PeptideIx) -> bool {
        self.cross_database(a, b) && self.same_silac_channel(&[a, b])
    }

    /// Create a new [`IndexedQuery`] for a specific [`ProcessedSpectrum`]
    ///
    /// All matches returned by the query will be within the specified tolerance
    /// parameters
    pub fn query(
        &self,
        precursor_mass: f32,
        precursor_tol: Tolerance,
        fragment_tol: Tolerance,
    ) -> IndexedQuery<'_> {
        let (precursor_lo, precursor_hi) = precursor_tol.bounds(precursor_mass);

        let (pre_idx_lo, pre_idx_hi) = binary_search_slice(
            &self.peptides,
            |p, bounds| p.monoisotopic.total_cmp(bounds),
            precursor_lo,
            precursor_hi,
        );

        IndexedQuery {
            db: self,
            precursor_mass,
            precursor_tol,
            fragment_tol,
            pre_idx_lo,
            pre_idx_hi,
        }
    }

    pub fn size(&self) -> usize {
        self.fragments.len()
    }

    pub fn buckets(&self) -> &[f32] {
        &self.min_value
    }

    pub fn serialize(&self) {
        use std::io::Write;
        let mut wtr = std::io::BufWriter::new(std::fs::File::create("fragments.bin").unwrap());
        for fragment in &self.fragments {
            let _ = wtr.write(&fragment.fragment_mz.to_le_bytes()).unwrap();
            let _ = wtr.write(&fragment.peptide_index.0.to_le_bytes()).unwrap();
        }
        wtr.flush().unwrap();

        let mut wtr = std::io::BufWriter::new(std::fs::File::create("peptides.csv").unwrap());
        writeln!(wtr, "peptide,proteins,monoisotopic,decoy").unwrap();
        for fragment in &self.peptides {
            writeln!(
                wtr,
                "{},{},{},{}",
                fragment.display(&self.modifications),
                fragment.proteins(&self.decoy_tag, self.generate_decoys),
                fragment.monoisotopic,
                fragment.decoy
            )
            .unwrap();
        }
        wtr.flush().unwrap();
    }
}

impl IndexedDatabase {
    /// Sites of `peptide` where mass offset `offset` may go
    pub fn offset_sites(&self, peptide: &Peptide, offset: usize) -> Vec<Site> {
        let o = &self.mass_offsets[offset];
        peptide.offset_sites(&o.sites, o.id)
    }

    /// `peptide` with a mass offset placed
    pub fn place(&self, peptide: &Peptide, placement: OffsetPlacement) -> Peptide {
        let o = &self.mass_offsets[placement.offset as usize];
        peptide.with_offset(placement.site, o.id, o.mass)
    }

    /// The peptide a PSM scored, with its mass offset (before
    /// [`IndexedDatabase::materialize_mass_offsets`])
    pub fn resolve(&self, feature: &Feature) -> Cow<'_, Peptide> {
        match feature.mass_offset {
            Some(placement) => Cow::Owned(self.place(&self[feature.peptide_idx], placement)),
            None => Cow::Borrowed(&self[feature.peptide_idx]),
        }
    }

    /// Give every PSM with a mass offset the index of its peptide form, so
    /// that FDR, quantification and the outputs see ordinary peptides: a form
    /// that is also indexed keeps that index, others are added once to
    /// `offset_peptides`. Returns the number of PSMs with an offset.
    pub fn materialize_mass_offsets(&mut self, features: &mut [Feature]) -> usize {
        type Key = (bool, Arc<[u8]>, Box<[Option<ModId>]>, Option<ModId>, Option<ModId>);
        let key = |p: &Peptide| -> Key { (p.decoy, p.sequence.clone(), p.modifications.clone(), p.nterm, p.cterm) };
        let base = self.peptides.len();
        let mut known: HashMap<Key, PeptideIx> = self
            .offset_peptides
            .iter()
            .enumerate()
            .map(|(i, p)| (key(p), PeptideIx((base + i) as u32)))
            .collect();
        let mut placed = 0;
        for feature in features.iter_mut() {
            let Some(placement) = feature.mass_offset.take() else { continue };
            placed += 1;
            let peptide = self.place(&self[feature.peptide_idx], placement);
            let lo = self.peptides.partition_point(|p| p.monoisotopic < peptide.monoisotopic - 0.01);
            let indexed = self.peptides[lo..]
                .iter()
                .take_while(|p| p.monoisotopic <= peptide.monoisotopic + 0.01)
                .position(|p| p.same_form(&peptide));
            feature.peptide_idx = match indexed {
                Some(i) => PeptideIx((lo + i) as u32),
                None => *known.entry(key(&peptide)).or_insert_with(|| {
                    self.offset_peptides.push(peptide);
                    PeptideIx((base + self.offset_peptides.len() - 1) as u32)
                }),
            };
        }
        placed
    }
}

impl std::ops::Index<PeptideIx> for IndexedDatabase {
    type Output = Peptide;

    fn index(&self, index: PeptideIx) -> &Self::Output {
        let index = index.0 as usize;
        match self.peptides.get(index) {
            Some(peptide) => peptide,
            None => &self.offset_peptides[index - self.peptides.len()],
        }
    }
}

pub struct IndexedQuery<'d> {
    db: &'d IndexedDatabase,
    precursor_mass: f32,
    precursor_tol: Tolerance,
    fragment_tol: Tolerance,
    pub pre_idx_lo: usize,
    pub pre_idx_hi: usize,
}

impl IndexedQuery<'_> {
    /// Search for a specified `fragment_mz` within the database
    pub fn page_search(&self, mass: f32) -> impl Iterator<Item = &Theoretical> {
        let (fragment_lo, fragment_hi) = self.fragment_tol.bounds(mass);
        let (precursor_lo, precursor_hi) = self.precursor_tol.bounds(self.precursor_mass);

        // Locate the left and right page indices that contain matching fragments
        // Note that we need to multiply by `bucket_size` to transform these into
        // indices that can be used with `self.db.fragments`
        let (left_idx, right_idx) = binary_search_slice(
            &self.db.min_value,
            |min, bounds| min.total_cmp(bounds),
            fragment_lo,
            fragment_hi,
        );

        // It is absolutely critical that we do not cross page boundaries!
        // If we do, we can no longer rely on total ordering of peptide_index (precursor m/z)
        (left_idx..right_idx).flat_map(move |page| {
            let left_idx = page * self.db.bucket_size;
            // Last chunk not guaranted to be modulo bucket size, make sure we don't
            // accidentally go out of bounds!
            let right_idx = ((page + 1) * self.db.bucket_size).min(self.db.fragments.len());

            // Narrow down into our region of interest, then perform another binary
            // search to further refine down to the slice of matching precursor mzs
            let slice = &&self.db.fragments[left_idx..right_idx];

            let (inner_left, inner_right) = binary_search_slice(
                slice,
                |frag, bounds| (frag.peptide_index.0 as usize).cmp(bounds),
                self.pre_idx_lo,
                self.pre_idx_hi,
            );

            // Finally, filter down our slice into exact matches only
            slice[inner_left..inner_right].iter().filter(move |frag| {
                // This looks somewhat complicated, but it's a consequence of
                // how the `binary_search_slice` function works - it will return
                // the set of indices that maximally cover the desired range - the exact
                // `left` and `right` indices may be valid, or just outside of the range.
                // Anything interior of `left` and `right` is guaranteed to be within the
                // precursor tolerance, so we just need to check the edge cases
                //
                // Previously, a direct lookup to check the mass of the current fragment was
                // performed, but the pointer indirection + float comparison can slow down
                // open searches by as much as 2x!!
                // e.g. used to be `self.db[frag.peptide_index].monoisotopic >= precursor_lo`
                (frag.peptide_index.0 > self.pre_idx_lo as u32
                    || (frag.peptide_index.0 == self.pre_idx_lo as u32
                        && self.db[frag.peptide_index].monoisotopic >= precursor_lo))
                    && (frag.peptide_index.0 < self.pre_idx_hi as u32
                        || (frag.peptide_index.0 == self.pre_idx_hi as u32
                            && self.db[frag.peptide_index].monoisotopic <= precursor_hi))
                    && frag.fragment_mz >= fragment_lo
                    && frag.fragment_mz <= fragment_hi
            })
        })
    }
}

/// Return the widest `left` and `right` indices into a `slice` (sorted by the
/// function `key`) such that all values between `low` and `high` are
/// contained in `slice[left..right]`
///
/// # Invariants
///
/// * `slice[left] <= low || left == 0`
/// * `slice[right] > high || right == slice.len()`
/// * `0 <= left <= right <= slice.len()`
#[inline]
pub fn binary_search_slice<T, F, S>(slice: &[T], key: F, low: S, high: S) -> (usize, usize)
where
    F: Fn(&T, &S) -> Ordering,
{
    let left_idx = slice
        .partition_point(|a| key(a, &low) == Ordering::Less)
        .saturating_sub(1);

    let right_idx =
        slice[left_idx..].partition_point(|a| key(a, &high) != Ordering::Greater) + left_idx;

    (left_idx, right_idx)
}

#[cfg(test)]
mod test {
    use std::sync::Arc;

    use super::*;

    #[test]
    fn binary_search_slice_smoke() {
        // Make sure that our query returns the maximal set of indices
        let data = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0];
        let bounds = binary_search_slice(&data, |a: &f64, b| a.total_cmp(b), 1.75, 3.5);
        assert_eq!(bounds, (1, 6));
        assert!(data[bounds.0] <= 1.75);
        assert_eq!(&data[bounds.0..bounds.1], &[1.5, 2.0, 2.5, 3.0, 3.5]);

        let bounds = binary_search_slice(&data, |a: &f64, b| a.total_cmp(b), 0.0, 5.0);
        assert_eq!(bounds, (0, data.len()));
    }

    #[test]
    fn binary_search_slice_run() {
        // Make sure that our query returns the maximal set of indices
        let data = [1.0, 1.5, 1.5, 1.5, 1.5, 2.0, 2.5, 3.0, 3.0, 3.5, 4.0];
        let (left, right) = binary_search_slice(&data, |a: &f64, b| a.total_cmp(b), 1.5, 3.25);
        assert!(data[left] <= 1.5);
        assert!(data[right] > 3.25);
        assert_eq!(
            &data[left..right],
            &[1.0, 1.5, 1.5, 1.5, 1.5, 2.0, 2.5, 3.0, 3.0]
        );
    }

    #[test]
    fn structured_variable_mod_config_round_trips() {
        let builder: Builder = serde_json::from_value(serde_json::json!({
            "fasta": "none",
            "variable_mods": {
                "M": [15.9949],
                "K": [
                    {"mass": 42.0106, "max_count": 1},
                    {"mass": 14.0157}
                ]
            },
            "max_variable_mods": 2,
            "max_combinations": 0
        }))
        .unwrap();

        let params = builder.make_parameters().unwrap();
        assert_eq!(params.max_variable_mods, 2);
        assert_eq!(params.max_combinations, Some(1));

        let resolved = &params.modifications;
        let mods = &resolved.variable_mods;
        let mass = |id| resolved.registry.mass(Some(id));
        assert_eq!(mods.len(), 3);
        assert_eq!(mods[0].0, ModificationSpecificity::Residue(b'K'));
        assert!((mass(mods[0].1) - 42.0106).abs() < 1e-4);
        assert_eq!(mods[0].2, Some(1));
        assert_eq!(mods[1].0, ModificationSpecificity::Residue(b'K'));
        assert!((mass(mods[1].1) - 14.0157).abs() < 1e-4);
        assert_eq!(mods[1].2, None);
        assert_eq!(mods[2].0, ModificationSpecificity::Residue(b'M'));
        assert!((mass(mods[2].1) - 15.9949).abs() < 1e-4);
        assert_eq!(mods[2].2, None);

        let serialized = serde_json::to_value(params).unwrap();
        let k_entries = &serialized["variable_mods"]["K"];
        assert!(k_entries[0].is_object());
        assert_eq!(k_entries[0]["max_count"], 1);
        // `{"mass": 14.0157}` carries nothing beyond the mass, so it is written back
        // in the bare-mass form
        assert_eq!(k_entries[1], serde_json::json!(14.0157f32));
        assert!(serialized["variable_mods"]["M"][0].is_number());
    }

    /// With `clip_n_term_met`, a protein starting M + small residue is also
    /// searched without its Met: the clipped peptides are protein N-terminal
    /// (they take protein N-terminal acetylation), the unclipped ones stay,
    /// and each has its decoy; a protein M + large residue is not clipped
    #[test]
    fn initiator_methionine_clipping() {
        let digest = |clip: bool| -> (Vec<Peptide>, Vec<String>) {
            let builder: Builder = serde_json::from_value(serde_json::json!({
                "fasta": "none",
                "enzyme": {"missed_cleavages": 1, "min_len": 4, "max_len": 30},
                "variable_mods": {"[": [42.010565]},
                "clip_n_term_met": clip
            }))
            .unwrap();
            let params = builder.make_parameters().unwrap();
            assert_eq!(params.clip_n_term_met, clip);
            let fasta = Fasta::parse(">sp|P1|ONE
MSDEREVAEAKLLR
>sp|P2|TWO
MLDEREVAEAKLLR
".into(), "rev_", true);
            let peptides = params.digest(&fasta);
            let registry = params.modifications.registry.clone();
            let mut names: Vec<String> =
                peptides.iter().filter(|p| !p.decoy).map(|p| p.display(&registry).to_string()).collect();
            names.sort();
            (peptides, names)
        };
        let (off, off_names) = digest(false);
        let (on, on_names) = digest(true);
        assert!(!off_names.iter().any(|p| p == "SDER" || p.ends_with("-SDER")), "{off_names:?}");
        for p in ["SDER", "SDEREVAEAK", "MSDER"] {
            assert!(on_names.contains(&p.to_string()), "{p} in {on_names:?}");
        }
        // the clipped N-termini take protein N-terminal acetylation
        assert!(on_names.iter().any(|p| p.starts_with('[') && p.ends_with("SDER")), "{on_names:?}");
        assert!(on_names.iter().any(|p| p.starts_with('[') && p.ends_with("-SDEREVAEAK")), "{on_names:?}");
        // M + L is not clipped
        assert!(!on_names.iter().any(|p| p.trim_start_matches(|c| c != 'L').starts_with("LDER") && !p.contains('M')), "{on_names:?}");
        let sder = on.iter().find(|p| !p.decoy && p.sequence.as_ref() == b"SDER").unwrap();
        assert_eq!(sder.position, crate::enzyme::Position::Nterm);
        // targets and decoys stay balanced; 2 clipped peptides, each also acetylated
        let count = |peptides: &[Peptide], decoy: bool| peptides.iter().filter(|p| p.decoy == decoy).count();
        assert_eq!(count(&on, false), count(&on, true));
        assert_eq!(count(&on, false), count(&off, false) + 4);
    }

    /// The estimate counts what the search would index
    #[test]
    fn estimate_counts_the_database() {
        let builder: Builder = serde_json::from_value(serde_json::json!({
            "fasta": "none",
            "enzyme": {"missed_cleavages": 1, "min_len": 5, "max_len": 30},
            "static_mods": {"C": 57.021464},
            "variable_mods": {"M": [15.994915], "[": [42.010565]},
            "peptide_min_mass": 300.0
        }))
        .unwrap();
        let params = builder.make_parameters().unwrap();
        let fasta = Fasta::parse(
            ">sp|P1|ONE\nMSDEREVAEAKLLRCGMKPEPTIDEK\n>sp|P2|TWO\nMLDEREVAEAKLLRQQQMMK\n".into(),
            "rev_",
            true,
        );
        let estimate = params.estimate(&fasta);
        let built = params.clone().build(fasta.clone());
        assert_eq!(estimate.proteins, 2);
        // the search merges duplicate forms, so the estimate is an upper bound, and close
        assert!(estimate.peptides >= built.peptides.len(), "{} vs {}", estimate.peptides, built.peptides.len());
        assert!(estimate.peptides <= built.peptides.len() + 2, "{} vs {}", estimate.peptides, built.peptides.len());
        assert!(estimate.fragments >= built.fragments.len());
        assert!(estimate.decoys * 2 >= estimate.peptides - 2);
        assert!(estimate.bytes > 8 * estimate.fragments as u64);
        assert_eq!(params.estimate(&Fasta::parse(String::new(), "rev_", true)).peptides, 0);
    }

    #[test]
    fn digestion() {
        let fasta = r#"
        >sp|AAAAA
        MEWKLEQSMREQALLKAQLTQLK
        >sp|BBBBB
        RMEWKLEQSMREQALLKAQLTQLK
        "#;

        let fasta = Fasta::parse(fasta.into(), "rev_", false);

        // Make sure that FASTA parsed OK
        assert_eq!(
            fasta.targets,
            vec![
                (
                    Arc::from("sp|AAAAA".to_string()),
                    "MEWKLEQSMREQALLKAQLTQLK".into()
                ),
                (
                    Arc::from("sp|BBBBB".to_string()),
                    "RMEWKLEQSMREQALLKAQLTQLK".into()
                ),
            ]
        );

        let variable_mods = [(
            ModificationSpecificity::ProteinN(None),
            vec![ModEntry::Mass(42.0)],
        )]
        .into_iter()
        .collect();
        let params = Parameters {
            bucket_size: 128,
            enzyme: EnzymeBuilder {
                missed_cleavages: Some(1),
                min_len: Some(6),
                max_len: Some(10),
                ..Default::default()
            },
            peptide_min_mass: 150.0,
            peptide_max_mass: 5000.0,
            ion_kinds: vec![Kind::B, Kind::Y],
            min_ion_index: 2,
            modifications: ResolvedModifications::new(&HashMap::default(), &variable_mods).unwrap(),
            max_variable_mods: 2,
            max_combinations: None,
            decoy_tag: "rev_".into(),
            generate_decoys: false,
            fasta: "none".into(),
            prefilter: false,
            prefilter_chunk_size: 0,
            prefilter_low_memory: true,
            clip_n_term_met: false,
            silac: None,
            labels: None,
            labelled_proteins: None,
            protein_labels: None,
            crosslink_databases: None,
            database_groups: None,
            site_library: None,
            library: None,
        };

        let peptides = params.digest(&fasta);

        let expected = [
            "EQALLK",
            "LEQSMR",
            "AQLTQLK",
            "MEWKLEQSMR",
            "[+42]-MEWKLEQSMR",
        ]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();

        let mods = params.modifications.registry.clone();
        let sequences = peptides
            .iter()
            .map(|p| p.display(&mods).to_string())
            .collect::<Vec<_>>();
        assert_eq!(expected, sequences);

        // All peptides are shared except for the protein N-term mod
        for peptide in &peptides[..4] {
            assert_eq!(peptide.proteins.len(), 2, "{:?}", peptide);
        }
        // Ensure that this mod is uniquely called as the first protein
        assert_eq!(
            peptides.last().unwrap().proteins,
            vec!["sp|AAAAA".to_string().into()]
        );
    }

    #[test]
    fn named_modifications() {
        let builder: Builder = serde_json::from_value(serde_json::json!({
            "fasta": "none",
            "static_mods": {
                "C": {"mass": 57.0215, "name": "UNIMOD:4"},
                "^": 229.1629
            },
            "variable_mods": {
                "M": [{"mass": 15.9949, "name": "UNIMOD:35"}],
                "[": [{"mass": 42.0106, "name": "Acetyl"}]
            },
            "generate_decoys": false,
            "peptide_min_mass": 0.0
        }))
        .unwrap();
        let params = builder.make_parameters().unwrap();

        // Registry ids are independent of HashMap iteration order
        let resolved = &params.modifications;
        assert_eq!(resolved.registry.len(), 4);

        // Serializes back to the config it was built from
        let serialized = serde_json::to_value(&params).unwrap();
        assert_eq!(
            serialized["static_mods"],
            serde_json::json!({"C": {"mass": 57.0215f32, "name": "UNIMOD:4"}, "^": 229.1629f32})
        );
        assert_eq!(
            serialized["variable_mods"],
            serde_json::json!({
                "M": [{"mass": 15.9949f32, "name": "UNIMOD:35"}],
                "[": [{"mass": 42.0106f32, "name": "Acetyl"}]
            })
        );
        let reparsed: Builder = serde_json::from_value(serde_json::json!({
            "fasta": "none",
            "static_mods": serialized["static_mods"],
            "variable_mods": serialized["variable_mods"],
            "generate_decoys": false,
            "peptide_min_mass": 0.0
        }))
        .unwrap();
        let reparsed = reparsed.make_parameters().unwrap().modifications;
        assert_eq!(reparsed.registry, resolved.registry);
        assert_eq!(reparsed.static_mods, resolved.static_mods);
        assert_eq!(reparsed.variable_mods, resolved.variable_mods);

        let fasta = Fasta::parse(">sp|AAAAA\nMCPEPK".into(), "rev_", false);
        let db = params.build(fasta);
        let mut sequences = db
            .peptides
            .iter()
            .map(|p| p.display(&db.modifications).to_string())
            .collect::<Vec<_>>();
        sequences.sort();
        assert_eq!(
            sequences,
            vec![
                "[+229.1629]-MC[UNIMOD:4]PEPK",
                "[+229.1629]-M[UNIMOD:35]C[UNIMOD:4]PEPK",
                "[Acetyl]-MC[UNIMOD:4]PEPK",
                "[Acetyl]-M[UNIMOD:35]C[UNIMOD:4]PEPK",
            ]
        );
    }

    #[test]
    fn invalid_modification_config() {
        let parse = |mods: serde_json::Value| {
            let mut config = serde_json::json!({"fasta": "none"});
            config
                .as_object_mut()
                .unwrap()
                .extend(mods.as_object().unwrap().clone());
            serde_json::from_value::<Builder>(config)
        };

        // Rejected while parsing
        assert!(parse(
            serde_json::json!({"variable_mods": {"M": [{"mass": 15.9949, "name": "Ox]"}]}})
        )
        .is_err());
        assert!(parse(serde_json::json!({"variable_mods": {"Z": [15.9949]}})).is_err());
        assert!(parse(serde_json::json!({"static_mods": {"^Z": 42.0}})).is_err());

        // Rejected when building parameters
        let builder =
            parse(serde_json::json!({"static_mods": {"C": {"mass": 57.0215, "max_count": 1}}}));
        assert_eq!(
            builder.unwrap().make_parameters().err(),
            Some(ModificationError::StaticMaxCount(
                ModificationSpecificity::Residue(b'C')
            ))
        );

        let variable_mods = (0..256)
            .map(|i| (i as f32).into())
            .collect::<Vec<serde_json::Value>>();
        let builder = parse(serde_json::json!({"variable_mods": {"M": variable_mods}}));
        assert_eq!(
            builder.unwrap().make_parameters().err(),
            Some(ModificationError::TooMany)
        );
    }

    #[test]
    fn two_crosslink_databases_are_merged_and_tagged() {
        let builder: Builder = serde_json::from_value(serde_json::json!({
            "crosslink_databases": {"database_1": "complex1.fasta", "database_2": "complex2.fasta"},
            "enzyme": {"min_len": 5, "max_len": 30},
            "peptide_min_mass": 400.0
        }))
        .unwrap();
        let mut params = builder.make_parameters().unwrap();
        assert_eq!(params.fasta_files(), vec!["complex1.fasta", "complex2.fasta"]);
        // a protein in both files (same accession and sequence) is read once
        let db1 = Fasta::parse(">sp|P1|ONE\nMSDLANSEKGLEKVSTR\n>sp|P3|SHARED\nAAWLDEKTTPYK\n".into(), "rev_", true);
        let db2 = Fasta::parse(">sp|P2|TWO\nTDYHFPNEKLIR\n>sp|P3|SHARED\nAAWLDEKTTPYK\n".into(), "rev_", true);
        let fasta = params.merge_databases(db1, db2);
        assert_eq!(fasta.targets.len(), 3);
        let db = params.build(fasta);
        assert!(db.two_databases());
        let group = |seq: &str, decoy: bool| {
            let ix = db.peptides.iter().position(|p| p.decoy == decoy && p.sequence.as_ref() == seq.as_bytes());
            db.database_group(PeptideIx(ix.unwrap_or_else(|| panic!("{seq}")) as u32)).unwrap()
        };
        assert_eq!(group("GLEKVSTR", false), DATABASE_1);
        assert_eq!(group("TDYHFPNEK", false), DATABASE_2);
        assert_eq!(group("AAWLDEK", false), DATABASE_1 | DATABASE_2);
        // decoys come from their target's database
        assert_eq!(group("GTSVKELR", true), DATABASE_1);
        assert_eq!(group("TENPFHYDK", true), DATABASE_2);
        let ix = |seq: &str| {
            PeptideIx(db.peptides.iter().position(|p| !p.decoy && p.sequence.as_ref() == seq.as_bytes()).unwrap() as u32)
        };
        assert!(db.cross_database(ix("GLEKVSTR"), ix("TDYHFPNEK")));
        assert!(!db.cross_database(ix("GLEKVSTR"), ix("GLEKVSTR")));
        assert!(db.cross_database(ix("AAWLDEK"), ix("AAWLDEK")));
        assert!(db.cross_database(ix("AAWLDEK"), ix("GLEKVSTR")));
    }

    /// Named modifications: motif and position rules decide the sites (a
    /// motif may reach past the peptide into its protein); decoys keep their
    /// target's sites; rules and neutral losses survive serialization
    #[test]
    fn named_modification_sites() {
        let builder: Builder = serde_json::from_value(serde_json::json!({
            "fasta": "none",
            "enzyme": {"missed_cleavages": 1, "min_len": 5},
            "variable_mods": {
                "motif:N-{P}-[ST]": [{"mass": 0.984016, "name": "UNIMOD:7"}],
                "internal_residue:K": [{"mass": 114.042927, "name": "UNIMOD:121"}],
                "S": [{"mass": 79.966331, "name": "UNIMOD:21", "neutral_losses": [97.976896]}],
                "first_residue:Q": [{"mass": -17.026549, "name": "UNIMOD:28"}]
            }
        }))
        .unwrap();
        let params = builder.make_parameters().unwrap();
        let fasta = crate::fasta::Fasta::parse(
            ">p\nMAGLLNASGKLLNPTEGKLLLLNRSAAAAAKGGKLLLAAKQLLLLGR\n".into(),
            "rev_",
            false,
        );
        let peptides = params.digest(&fasta);
        let show = |p: &Peptide| p.display(&params.modifications.registry).to_string();
        let targets: Vec<String> = peptides.iter().filter(|p| !p.decoy).map(show).collect();
        let has = |name: &str| targets.iter().any(|t| t == name);
        // the sequon N-A-S; not N-P-T
        assert!(has("MAGLLN[UNIMOD:7]ASGK"), "{targets:?}");
        assert!(!targets.iter().any(|t| t.contains("N[UNIMOD:7]PT")), "{targets:?}");
        // N-R|S: the S/T is the next protein residue, past the peptide
        assert!(has("LLLLN[UNIMOD:7]R"), "{targets:?}");
        // GlyGly on an internal lysine only
        assert!(has("GGK[UNIMOD:121]LLLAAK"), "{targets:?}");
        assert!(!targets.iter().any(|t| t.ends_with("K[UNIMOD:121]")), "{targets:?}");
        // named terminal site: pyro-glu on a peptide-initial Q
        assert!(has("Q[UNIMOD:28]LLLLGR"), "{targets:?}");
        // decoys are reversed targets: they keep the sites
        assert!(peptides.iter().any(|p| p.decoy && show(p).contains("N[UNIMOD:7]")));

        // phospho's loss is registered on S
        let phospho = params.modifications.variable_mods.iter().find(|(s, _, _)| *s == ModificationSpecificity::Residue(b'S')).unwrap().1;
        let losses: Vec<f32> = params.modifications.registry.neutral_losses(phospho, b'S').collect();
        assert!((losses[0] - 97.976896).abs() < 1e-4);
        assert_eq!(params.modifications.registry.neutral_losses(phospho, b'T').count(), 0);

        // the configuration written back keeps rules and losses, and reads again
        let written = serde_json::to_value(&params.modifications).unwrap();
        assert_eq!(written["variable_mods"]["N"][0]["motif"], "N-{P}-[ST]");
        assert_eq!(written["variable_mods"]["K"][0]["position"], "internal");
        assert!((written["variable_mods"]["S"][0]["neutral_losses"][0].as_f64().unwrap() - 97.976896).abs() < 1e-4);
        let again: Builder = serde_json::from_value(serde_json::json!({
            "fasta": "none",
            "enzyme": {"missed_cleavages": 1, "min_len": 5},
            "variable_mods": written["variable_mods"].clone()
        }))
        .unwrap();
        let again = again.make_parameters().unwrap();
        assert_eq!(again.modifications.variable_rules, params.modifications.variable_rules);
        let mut a: Vec<String> = again.digest(&fasta).iter().map(|p| p.display(&again.modifications.registry).to_string()).collect();
        let mut b: Vec<String> = peptides.iter().map(show).collect();
        a.sort();
        b.sort();
        assert_eq!(a, b);

        // rules that cannot apply
        for bad in [
            serde_json::json!({"^": [{"mass": 1.0, "position": "internal"}]}),
            serde_json::json!({"K": [{"mass": 1.0, "motif": "N-x-[ST]"}]}),
            serde_json::json!({"S": [{"mass": 1.0, "neutral_losses": [-3.0]}]}),
        ] {
            let builder: Builder =
                serde_json::from_value(serde_json::json!({"fasta": "none", "variable_mods": bad})).unwrap();
            assert!(builder.make_parameters().is_err(), "{bad}");
        }
        let unknown = serde_json::from_value::<Builder>(serde_json::json!({"fasta": "none", "variable_mods": {"elsewhere:K": [1.0]}}));
        assert!(unknown.is_err());

        // a static modification with a rule: on internal lysines only
        let builder: Builder = serde_json::from_value(serde_json::json!({
            "fasta": "none",
            "enzyme": {"missed_cleavages": 1, "min_len": 5},
            "static_mods": {"internal_residue:K": 42.010565}
        }))
        .unwrap();
        let params = builder.make_parameters().unwrap();
        let shown: Vec<String> = params.digest(&fasta).iter().filter(|p| !p.decoy).map(|p| p.display(&params.modifications.registry).to_string()).collect();
        assert!(shown.iter().any(|t| t.starts_with("GGK[+42.01") && t.ends_with("]LLLAAK")), "{shown:?}");
        assert!(!shown.iter().any(|t| t.ends_with(']')), "no C-terminal lysine carries it: {shown:?}");
    }
}

#[cfg(test)]
mod site_library_tests {
    use super::*;

    fn parameters(variable: serde_json::Value, library: bool) -> Result<Parameters, ModificationError> {
        let mut config = serde_json::json!({
            "fasta": "none",
            "enzyme": {"missed_cleavages": 0, "min_len": 5},
            "variable_mods": variable,
        });
        if library {
            config["site_library"] = serde_json::json!({"path": "sites.tsv"});
        }
        let builder: Builder = serde_json::from_value(config).unwrap();
        builder.make_parameters()
    }

    fn phospho() -> serde_json::Value {
        serde_json::json!({
            "S": [{"mass": 79.966331, "name": "UNIMOD:21", "site_mode": "library"}],
            "T": [{"mass": 79.966331, "name": "UNIMOD:21", "site_mode": "library"}],
            "M": [15.9949]
        })
    }

    /// Phospho only at the library's sites, on any protein of a peptide;
    /// decoys keep their target's site
    #[test]
    fn library_sites_limit_a_modification() {
        let mut params = parameters(phospho(), true).unwrap();
        // P1: MSSTPLGEKSASTPAAR; P2 shares SASTPAAR at 1-8
        let library = params
            .load_site_library(
                "protein\tposition\tresidue\tmodification\n\
                 P1\t3\tS\tUNIMOD:21\n\
                 P2\t4\tT\t+79.96633\n\
                 P1\t2\tS\tUNIMOD:35\n\
                 P1\t5\tP\tUNIMOD:21\n",
            )
            .unwrap();
        assert_eq!((library.rows, library.used), (4, 2), "an unknown modification and a P are not used");
        let json = serde_json::to_value(&params).unwrap();
        assert_eq!(json["variable_mods"]["S"][0]["site_mode"], "library");
        assert_eq!(json["site_library"]["path"], "sites.tsv");

        let fasta = crate::fasta::Fasta::parse(">P1\nMSSTPLGEKSASTPAAR\n>P2\nSASTPAARWEEK\n".into(), "rev_", true);
        let db = params.clone().build(fasta);
        let phospho = params.modifications.variable_mods.iter().find(|m| m.0 == ModificationSpecificity::Residue(b'S')).unwrap().1;
        let sites = |sequence: &[u8], decoy: bool| -> Vec<Vec<usize>> {
            let mut out: Vec<Vec<usize>> = db
                .peptides
                .iter()
                .filter(|p| p.decoy == decoy && p.sequence.as_ref() == sequence)
                .map(|p| (0..p.sequence.len()).filter(|&i| p.modifications[i] == Some(phospho)).collect())
                .collect();
            out.sort();
            out
        };
        // MSSTPLGEK: only S3 (index 2), with or without oxidised M
        assert_eq!(sites(b"MSSTPLGEK", false), vec![vec![], vec![], vec![2], vec![2]]);
        // SASTPAAR: T4 from P2's row (index 3)
        assert_eq!(sites(b"SASTPAAR", false), vec![vec![], vec![3]]);
        // its decoy SAAPTSAR: the same residue, reversed (index 4)
        assert_eq!(sites(b"SAAPTSAR", true), vec![vec![], vec![4]]);
    }

    #[test]
    fn library_needs_a_library_and_residues() {
        assert!(parameters(phospho(), false).is_err(), "no site_library");
        let terminal = serde_json::json!({"^": [{"mass": 42.0106, "site_mode": "library"}]});
        assert!(parameters(terminal, true).is_err());
        let builder: Builder = serde_json::from_value(serde_json::json!({
            "fasta": "none",
            "static_mods": {"C": {"mass": 57.021465, "site_mode": "library"}},
        }))
        .unwrap();
        assert!(builder.make_parameters().is_err(), "static");
    }
}
