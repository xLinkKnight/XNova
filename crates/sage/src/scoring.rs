use crate::database::{IndexedDatabase, OffsetPlacement, PeptideIx, Theoretical};
use crate::heap::bounded_min_heapify;
use crate::ion_series::{IonSeries, Kind};
use crate::mass::{Tolerance, NEUTRON, PROTON};
use crate::peptide::Peptide;
use crate::spectrum::{Precursor, ProcessedSpectrum};
use serde::{Deserialize, Serialize};
use std::ops::AddAssign;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub enum ScoreType {
    SageHyperScore,
    OpenMSHyperScore,
}

/// Structure to hold temporary scores
#[derive(Copy, Clone, Default, Debug, PartialEq)]
struct Score {
    peptide: PeptideIx,
    matched_b: u16,
    matched_y: u16,
    summed_b: f32,
    summed_y: f32,
    longest_b: usize,
    longest_y: usize,
    hyperscore: f64,
    ppm_difference: f32,
    precursor_charge: u8,
    isotope_error: i8,
}

impl Eq for Score {}

impl PartialOrd for Score {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Score {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.hyperscore.total_cmp(&other.hyperscore)
    }
}

/// Preliminary score - # of matched peaks for each candidate peptide
#[derive(Copy, Clone, Default, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct PreScore {
    matched: u16,
    peptide: PeptideIx,
    precursor_charge: u8,
    isotope_error: i8,
    /// Mass offset tried (1-based into [`IndexedDatabase::mass_offsets`]), 0
    /// for none
    offset: u8,
}

/// A fully scored candidate: score, matched fragments, and the mass offset
/// placed on the peptide (if any)
type Scored = (Score, Option<Fragments>, Option<OffsetPlacement>);

/// Store preliminary scores & stats for first pass search for a query spectrum
#[derive(Clone, Default)]
struct InitialHits {
    matched_peaks: usize,
    // Number of peptide candidates with > 0 matched peaks
    scored_candidates: usize,
    preliminary: Vec<PreScore>,
}

impl AddAssign<InitialHits> for InitialHits {
    fn add_assign(&mut self, rhs: InitialHits) {
        self.matched_peaks += rhs.matched_peaks;
        self.scored_candidates += rhs.scored_candidates;

        self.preliminary.extend(rhs.preliminary);
    }
}

#[derive(Serialize, Clone, Debug, Default)]
/// Features of a candidate peptide spectrum match
pub struct Feature {
    #[serde(skip_serializing)]
    pub peptide_idx: PeptideIx,
    // psm_id help to match with matched fragments table.
    pub psm_id: usize,
    pub peptide_len: usize,
    /// Spectrum id
    pub spec_id: String,
    /// File identifier
    pub file_id: usize,
    /// PSM rank
    pub rank: u32,
    /// Target/Decoy label, -1 is decoy, 1 is target
    pub label: i32,
    /// Experimental mass
    pub expmass: f32,
    /// Calculated mass
    pub calcmass: f32,
    /// Reported precursor charge
    pub charge: u8,
    /// Retention time
    pub rt: f32,
    /// Globally aligned retention time
    pub aligned_rt: f32,
    /// Predicted RT, if enabled
    pub predicted_rt: f32,
    /// Difference between predicted & observed RT
    pub delta_rt_model: f32,
    /// Ion mobility
    pub ims: f32,
    /// Predicted ion mobility, if enabled
    pub predicted_ims: f32,
    /// Difference between predicted & observed ion mobility
    pub delta_ims_model: f32,
    /// Difference between expmass and calcmass
    pub delta_mass: f32,
    /// C13 isotope error
    pub isotope_error: f32,
    /// Average ppm delta mass for matched fragments
    pub average_ppm: f32,
    /// X!Tandem hyperscore
    pub hyperscore: f64,
    /// Difference between hyperscore of this candidate, and the next best candidate
    pub delta_next: f64,
    /// Difference between hyperscore of this candidate, and the best candidate
    pub delta_best: f64,
    /// Number of matched theoretical fragment ions
    pub matched_peaks: u32,
    /// Longest b-ion series
    pub longest_b: u32,
    /// Longest y-ion series
    pub longest_y: u32,
    /// Longest y-ion series, divided by peptide length
    pub longest_y_pct: f32,
    /// Number of missed cleavages
    pub missed_cleavages: u8,
    /// Fraction of matched MS2 intensity
    pub matched_intensity_pct: f32,
    /// Number of scored candidates for this spectrum
    pub scored_candidates: u32,
    /// Probability of matching exactly N peaks across all candidates Pr(x=k)
    pub poisson: f64,
    /// Combined score from linear discriminant analysis, used for FDR calc
    pub discriminant_score: f32,
    /// Posterior error probability for this PSM / local FDR
    pub posterior_error: f32,
    /// Assigned q_value
    pub spectrum_q: f32,
    pub peptide_q: f32,
    pub protein_q: f32,
    pub protein_group_q: f32,

    pub ms2_intensity: f32,

    pub protein_groups: Option<String>,
    pub num_protein_groups: u32,

    pub fragments: Option<Fragments>,

    /// Site localization of the variable modifications (rank 1, with
    /// `localization` on)
    #[serde(skip_serializing)]
    pub localization: Option<crate::localization::Localization>,

    /// The mass offset placed on `peptide_idx` (until
    /// [`IndexedDatabase::materialize_mass_offsets`] gives the placed peptide
    /// an index of its own)
    #[serde(skip_serializing)]
    pub mass_offset: Option<OffsetPlacement>,
}

/// Matching Fragment details
#[derive(Serialize, Default, Clone, Debug)]
pub struct Fragments {
    /// Observed fragment charge state.
    #[serde(skip_serializing)]
    pub charges: Vec<i32>,
    pub kinds: Vec<Kind>,
    pub fragment_ordinals: Vec<i32>,
    pub intensities: Vec<f32>,
    pub mz_calculated: Vec<f32>,
    pub mz_experimental: Vec<f32>,
}

static PSM_COUNTER: AtomicUsize = AtomicUsize::new(1);

fn increment_psm_counter() -> usize {
    PSM_COUNTER.fetch_add(1, Ordering::Relaxed)
}

/// Stirling's approximation for log factorial
fn lnfact(n: u16) -> f64 {
    if n == 0 {
        1.0
    } else {
        let n = n as f64;
        n * n.ln() - n + 0.5 * n.ln() + 0.5 * (std::f64::consts::PI * 2.0 * n).ln()
    }
}

impl ScoreType {
    pub fn score(&self, matched_b: u16, matched_y: u16, summed_b: f32, summed_y: f32) -> f64 {
        let score = match self {
            // Calculate the X!Tandem hyperscore
            Self::SageHyperScore => {
                let i = (summed_b + 1.0) as f64 * (summed_y + 1.0) as f64;

                i.ln() + lnfact(matched_b) + lnfact(matched_y)
            }
            // Calculate the OpenMS flavour hyperscore
            Self::OpenMSHyperScore => {
                let summed_intensity = summed_b + summed_y;

                summed_intensity.ln_1p() as f64 + lnfact(matched_b) + lnfact(matched_y)
            }
        };
        if score.is_finite() {
            score
        } else {
            255.0
        }
    }
}

impl Score {
    /// Calculate the hyperscore for a given PSM choosing between implementations based on `score_type`
    fn hyperscore(&self, score_type: ScoreType) -> f64 {
        score_type.score(self.matched_b, self.matched_y, self.summed_b, self.summed_y)
    }
}

pub struct Scorer<'db> {
    pub db: &'db IndexedDatabase,
    pub precursor_tol: Tolerance,
    pub fragment_tol: Tolerance,
    /// What is the minimum number of matched b and y ion peaks to report PSMs for?
    pub min_matched_peaks: u16,
    /// Precursor isotope error lower bounds (e.g. -1)
    pub min_isotope_err: i8,
    /// Precursor isotope error upper bounds (e.g. 3)
    pub max_isotope_err: i8,
    pub min_precursor_charge: u8,
    pub max_precursor_charge: u8,
    pub override_precursor_charge: bool,
    pub max_fragment_charge: Option<u8>,
    pub chimera: bool,
    pub report_psms: usize,

    // Rather than use a fixed precursor tolerance, dynamically alter
    // the precursor tolerance window based on MS2 isolation window and charge
    pub wide_window: bool,
    pub annotate_matches: bool,
    pub score_type: ScoreType,
}

#[inline(always)]
/// Calculate upper bound (excluded) of the charge state range to use for
/// searching fragment ions (1..N)
/// If user has configured max_fragment_charge, potentially override precursor
/// charge
fn max_fragment_charge(max_fragment_charge: Option<u8>, precursor_charge: u8) -> u8 {
    precursor_charge
        .min(
            max_fragment_charge
                .map(|c| c + 1)
                .unwrap_or(precursor_charge),
        )
        .max(2)
}

impl<'db> Scorer<'db> {
    /// Perform a quick first-pass scoring, where we consider a peptide "identified"
    /// if it meets the following criterion:
    ///  * prefilter_low_memory = true: in the top `report_psms` hits for a spectrum
    ///  * prefilter_low_memory = false: has at least `min_matched_peaks` fragment ion matches
    /// * `keep`: A vector of atomic bools is used to maintain an identification list across scans
    pub fn quick_score(
        &self,
        query: &ProcessedSpectrum,
        prefilter_low_memory: bool,
        keep: &[AtomicBool],
    ) {
        assert_eq!(
            query.level, 2,
            "internal bug, trying to score a non-MS2 scan!"
        );
        let precursor = query.precursors.first().unwrap_or_else(|| {
            panic!("missing MS1 precursor for {}", query.id);
        });
        let hits = self.initial_hits(query, precursor);

        if prefilter_low_memory {
            let mut score_vector = hits
                .preliminary
                .iter()
                .filter_map(|pre| {
                    if pre.peptide == PeptideIx::default() {
                        return None;
                    }
                    let score = self.score_hypotheses(query, pre, false).into_iter().map(|s| s.0).max()?;
                    if (score.matched_b + score.matched_y) < self.min_matched_peaks {
                        return None;
                    }
                    Some(score)
                })
                .collect::<Vec<_>>();

            let k = self.report_psms.min(score_vector.len());
            bounded_min_heapify(&mut score_vector, k);
            for score in &score_vector[..k] {
                keep[score.peptide.0 as usize].store(true, Ordering::Relaxed);
            }
        } else {
            for pre in &hits.preliminary {
                if pre.peptide != PeptideIx::default() {
                    keep[pre.peptide.0 as usize].store(true, Ordering::Relaxed);
                }
            }
        }
    }

    pub fn score(&self, query: &ProcessedSpectrum) -> Vec<Feature> {
        assert_eq!(
            query.level, 2,
            "internal bug, trying to score a non-MS2 scan!"
        );
        match self.chimera {
            true => self.score_chimera_fast(query),
            false => self.score_standard(query),
        }
    }

    /// Perform a k-select and truncation of an [`InitialHits`] list.
    ///
    /// Determine how many candidates to actually calculate hyperscore for.
    /// Hyperscore is relatively computationally expensive, so we don't want
    /// to calculate it for every possible candidate (100s - 10,000s depending on search)
    /// when we are only going to report a few PSMs. But we also want to calculate
    /// it for enough candidates that we don't accidentally miss the best hit!
    ///
    /// Given that hyperscore is dominated by the number of matched peaks, it seems
    /// reasonable to assume that the highest hyperscore will belong to one of the
    /// top 50 candidates sorted by # of matched peaks.
    fn trim_hits(&self, hits: &mut InitialHits) {
        let k = 50.clamp(
            (self.report_psms * 2).min(hits.preliminary.len()),
            hits.preliminary.len(),
        );
        bounded_min_heapify(&mut hits.preliminary, k);
        hits.preliminary.truncate(k);
    }

    /// [`Scorer::trim_hits`] for each mass offset tried on its own: an
    /// offset's candidates match fragments both shifted and not, so their
    /// counts do not compare with those of the plain search
    fn trim_hits_by_offset(&self, hits: &mut InitialHits) {
        if self.db.mass_offsets.is_empty() {
            return self.trim_hits(hits);
        }
        let mut preliminary = std::mem::take(&mut hits.preliminary);
        preliminary.sort_unstable_by_key(|p| p.offset);
        let mut trimmed = Vec::with_capacity(preliminary.len());
        for group in preliminary.chunk_by_mut(|a, b| a.offset == b.offset) {
            let k = 50.clamp((self.report_psms * 2).min(group.len()), group.len());
            bounded_min_heapify(group, k);
            trimmed.extend_from_slice(&group[..k]);
        }
        hits.preliminary = trimmed;
    }

    /// Preliminary Score, return # of matched peaks per candidate
    /// Returned hits are guaranteed to be the top-K hits (see above comment)
    /// from among all potential candidates, but the returned vector is not
    /// in sorted order.
    fn matched_peaks_with_isotope(
        &self,
        query: &ProcessedSpectrum,
        precursor_mass: f32,
        precursor_charge: u8,
        precursor_tol: Tolerance,
        isotope_error: i8,
        offset: u8,
    ) -> InitialHits {
        let observed = precursor_mass - isotope_error as f32 * NEUTRON;
        // a mass offset: candidates at the precursor mass less the offset (the
        // window as wide as at the observed mass), and fragments carrying it
        // shifted by it
        let (mass, tolerance, shift) = match offset {
            0 => (observed, precursor_tol, None),
            o => {
                let o = &self.db.mass_offsets[o as usize - 1];
                let (lo, hi) = precursor_tol.bounds(observed);
                (observed - o.mass, Tolerance::Da(lo - observed, hi - observed), Some(o.mass))
            }
        };
        let candidates = self.db.query(mass, tolerance, self.fragment_tol);

        let max_fragment_charge = max_fragment_charge(self.max_fragment_charge, precursor_charge);
        // Allocate space for all potential candidates - many potential candidates
        let potential = candidates.pre_idx_hi - candidates.pre_idx_lo + 1;
        let mut hits = InitialHits {
            matched_peaks: 0,
            scored_candidates: 0,
            preliminary: vec![PreScore::default(); potential],
        };

        let mut count = |frag: &Theoretical| {
            let idx = frag.peptide_index.0 as usize - candidates.pre_idx_lo;
            let sc = &mut hits.preliminary[idx];
            if sc.matched == 0 {
                hits.scored_candidates += 1;
                sc.precursor_charge = precursor_charge;
                sc.peptide = frag.peptide_index;
                sc.isotope_error = isotope_error;
                sc.offset = offset;
            }

            sc.matched += 1;
            hits.matched_peaks += 1;
        };
        for peak_mass in query.masses.iter() {
            for charge in 1..max_fragment_charge {
                let mass = peak_mass * charge as f32;
                candidates.page_search(mass).for_each(&mut count);
                if let Some(shift) = shift {
                    candidates.page_search(mass - shift).for_each(&mut count);
                }
            }
        }
        if hits.matched_peaks == 0 {
            return hits;
        }

        self.trim_hits(&mut hits);
        hits
    }

    /// Preliminary scores without a mass offset, then with each
    fn matched_peaks(
        &self,
        query: &ProcessedSpectrum,
        precursor_mass: f32,
        precursor_charge: u8,
        precursor_tol: Tolerance,
    ) -> InitialHits {
        let mut hits = InitialHits::default();
        for offset in 0..=self.db.mass_offsets.len() as u8 {
            hits += self.matched_peaks_with_offset(query, precursor_mass, precursor_charge, precursor_tol, offset);
        }
        hits
    }

    fn matched_peaks_with_offset(
        &self,
        query: &ProcessedSpectrum,
        precursor_mass: f32,
        precursor_charge: u8,
        precursor_tol: Tolerance,
        offset: u8,
    ) -> InitialHits {
        if self.min_isotope_err != self.max_isotope_err {
            let mut hits = (self.min_isotope_err..=self.max_isotope_err).fold(
                InitialHits::default(),
                |mut hits, isotope| {
                    hits += self.matched_peaks_with_isotope(
                        query,
                        precursor_mass,
                        precursor_charge,
                        precursor_tol,
                        isotope,
                        offset,
                    );
                    hits
                },
            );
            self.trim_hits(&mut hits);
            hits
        } else {
            self.matched_peaks_with_isotope(
                query,
                precursor_mass,
                precursor_charge,
                precursor_tol,
                0,
                offset,
            )
        }
    }

    fn initial_hits(&self, query: &ProcessedSpectrum, precursor: &Precursor) -> InitialHits {
        // Sage operates on masses without protons; [M] instead of [MH+]
        let mz = precursor.mz - PROTON;

        // Search in wide-window/DIA mode
        let mut hits = if self.wide_window {
            (self.min_precursor_charge..=self.max_precursor_charge).fold(
                InitialHits::default(),
                |mut hits, precursor_charge| {
                    let precursor_mass = mz * precursor_charge as f32;
                    let precursor_tol = precursor
                        .isolation_window
                        .unwrap_or(Tolerance::Da(-2.4, 2.4))
                        * precursor_charge as f32;
                    hits +=
                        self.matched_peaks(query, precursor_mass, precursor_charge, precursor_tol);
                    hits
                },
            )
        } else if precursor.charge.is_some() && !self.override_precursor_charge {
            let charge = precursor.charge.unwrap();
            // Charge state is already annotated for this precusor, only search once
            let precursor_mass = mz * charge as f32;
            self.matched_peaks(query, precursor_mass, charge, self.precursor_tol)
        } else {
            // Not all selected ion precursors have charge states annotated (or user has set
            // `override_precursor_charge`)
            // assume it could be z=2, z=3, z=4 and search all three
            (self.min_precursor_charge..=self.max_precursor_charge).fold(
                InitialHits::default(),
                |mut hits, precursor_charge| {
                    let precursor_mass = mz * precursor_charge as f32;
                    hits += self.matched_peaks(
                        query,
                        precursor_mass,
                        precursor_charge,
                        self.precursor_tol,
                    );
                    hits
                },
            )
        };
        self.trim_hits_by_offset(&mut hits);
        hits
    }

    /// Score a single [`ProcessedSpectrum`] against the database
    pub fn score_standard(&self, query: &ProcessedSpectrum) -> Vec<Feature> {
        let precursor = query.precursors.first().unwrap_or_else(|| {
            panic!("missing MS1 precursor for {}", query.id);
        });

        let hits = self.initial_hits(query, precursor);
        let mut features = Vec::with_capacity(self.report_psms);
        self.build_features(query, precursor, &hits, self.report_psms, &mut features);
        features
    }

    /// Given a set of [`InitialHits`] against a query spectrum, prepare N=`report_psms`
    /// best PSMs ([`Feature`])
    fn build_features(
        &self,
        query: &ProcessedSpectrum,
        precursor: &Precursor,
        hits: &InitialHits,
        report_psms: usize,
        features: &mut Vec<Feature>,
    ) {
        let mut score_vector = hits
            .preliminary
            .iter()
            .filter(|score| score.peptide != PeptideIx::default())
            .flat_map(|pre| self.score_hypotheses(query, pre, self.annotate_matches))
            .filter(|s| (s.0.matched_b + s.0.matched_y) >= self.min_matched_peaks)
            .collect::<Vec<_>>();

        // Hyperscore is our primary score function for PSMs
        score_vector.sort_by(|a, b| b.0.hyperscore.total_cmp(&a.0.hyperscore));

        // Expected value for poisson distribution
        // (average # of matches peaks/peptide candidate)
        let lambda = hits.matched_peaks as f64 / hits.scored_candidates as f64;

        // Sage operates on masses without protons; [M] instead of [MH+]
        let mz = precursor.mz - PROTON;

        for idx in 0..report_psms.min(score_vector.len()) {
            let score = score_vector[idx].0;
            let fragments: Option<Fragments> = score_vector[idx].1.take();
            let mass_offset = score_vector[idx].2;
            let psm_id = increment_psm_counter();

            let placed = mass_offset.map(|p| self.db.place(&self.db[score.peptide], p));
            let peptide = placed.as_ref().unwrap_or(&self.db[score.peptide]);
            let precursor_mass = mz * score.precursor_charge as f32;

            let next = score_vector
                .get(idx + 1)
                .map(|score| score.0.hyperscore)
                .unwrap_or_default();

            let best = score_vector
                .first()
                .map(|score| score.0.hyperscore)
                .expect("we know that index 0 is valid");

            // Poisson distribution log10 probability mass function
            // Computed directly in log space to avoid overflow from lambda.powi(k)
            // log10(PMF) = (k*ln(lambda) - lambda - lnfact(k)) / ln(10)
            let k = score.matched_b + score.matched_y;
            let log10_poisson =
                (k as f64 * lambda.ln() - lambda - lnfact(k)) / std::f64::consts::LN_10;

            let isotope_error = score.isotope_error as f32 * NEUTRON;
            let delta_mass = (precursor_mass - peptide.monoisotopic - isotope_error) * 2E6
                / (precursor_mass - isotope_error + peptide.monoisotopic);

            // let (num_proteins, proteins) = self.db.assign_proteins(peptide);

            features.push(Feature {
                // Identifiers
                psm_id,
                peptide_idx: score.peptide,
                spec_id: query.id.clone(),
                file_id: query.file_id,
                rank: idx as u32 + 1,
                label: peptide.label(),
                expmass: precursor_mass,
                calcmass: peptide.monoisotopic,
                // Features
                charge: score.precursor_charge,
                rt: query.scan_start_time,
                ims: query
                    .precursors
                    .first()
                    .unwrap()
                    .inverse_ion_mobility
                    .unwrap_or(0.0),
                delta_mass,
                isotope_error,
                average_ppm: score.ppm_difference,
                hyperscore: score.hyperscore,
                delta_next: score.hyperscore - next,
                delta_best: best - score.hyperscore,
                matched_peaks: k as u32,
                matched_intensity_pct: 100.0 * (score.summed_b + score.summed_y)
                    / query.total_ion_current,
                poisson: if log10_poisson.is_finite() {
                    log10_poisson
                } else {
                    f64::NEG_INFINITY
                },
                longest_b: score.longest_b as u32,
                longest_y: score.longest_y as u32,
                longest_y_pct: score.longest_y as f32 / (peptide.sequence.len() as f32),
                peptide_len: peptide.sequence.len(),
                scored_candidates: hits.scored_candidates as u32,
                missed_cleavages: peptide.missed_cleavages,

                // Outputs
                discriminant_score: 0.0,
                posterior_error: 1.0,
                spectrum_q: 1.0,
                protein_q: 1.0,
                peptide_q: 1.0,
                predicted_rt: 0.0,
                predicted_ims: 0.0,
                aligned_rt: query.scan_start_time,
                delta_rt_model: 0.999,
                delta_ims_model: 0.999,
                ms2_intensity: score.summed_b + score.summed_y,

                //Fragments
                protein_groups: None,
                num_protein_groups: 0,
                fragments,
                protein_group_q: 1.0,
                localization: None,
                mass_offset,
            })
        }
    }

    /// Remove peaks matching a PSM from a query spectrum
    fn remove_matched_peaks(&self, query: &mut ProcessedSpectrum, psm: &Feature) {
        let peptide = self.db.resolve(psm);
        let fragments = self
            .db
            .ion_kinds
            .iter()
            .flat_map(|kind| IonSeries::new(&peptide, *kind, &self.db.modifications));

        let max_fragment_charge = max_fragment_charge(self.max_fragment_charge, psm.charge);

        // Remove MS2 peaks matched by previous match
        let mut to_remove = Vec::new();
        for frag in fragments {
            for charge in 1..max_fragment_charge {
                // Experimental peaks are multipled by charge, therefore theoretical are divided
                if let Some(peak_idx) = crate::spectrum::select_most_intense_peak(
                    &query.masses,
                    &query.intensities,
                    frag.monoisotopic_mass / charge as f32,
                    self.fragment_tol,
                    None,
                ) {
                    to_remove.push((
                        query.masses[peak_idx],
                        query.intensities[peak_idx],
                        query.charges[peak_idx],
                    ));
                }
            }
        }

        let mut masses = Vec::with_capacity(query.masses.len());
        let mut intensities = Vec::with_capacity(query.intensities.len());
        let mut charges = Vec::with_capacity(query.charges.len());
        let mut mobilities = Vec::with_capacity(query.mobilities.len());

        for idx in 0..query.masses.len() {
            let peak = (
                query.masses[idx],
                query.intensities[idx],
                query.charges[idx],
            );
            if !to_remove.contains(&peak) {
                masses.push(query.masses[idx]);
                intensities.push(query.intensities[idx]);
                charges.push(query.charges[idx]);
                if !query.mobilities.is_empty() {
                    mobilities.push(query.mobilities[idx]);
                }
            }
        }

        query.masses = masses;
        query.intensities = intensities;
        query.charges = charges;
        query.mobilities = mobilities;
        query.total_ion_current = query.intensities.iter().sum::<f32>();
    }

    /// Return multiple PSMs for each spectra - first is the best match, second PSM is the best match
    /// after all theoretical peaks assigned to the best match are removed, etc
    pub fn score_chimera_fast(&self, query: &ProcessedSpectrum) -> Vec<Feature> {
        let precursor = query.precursors.first().unwrap_or_else(|| {
            panic!("missing MS1 precursor for {}", query.id);
        });

        let mut query = query.clone();
        let hits = self.initial_hits(&query, precursor);

        let mut candidates: Vec<Feature> = Vec::with_capacity(self.report_psms);

        let mut prev = 0;
        while candidates.len() < self.report_psms {
            self.build_features(&query, precursor, &hits, 1, &mut candidates);
            if candidates.len() > prev {
                if let Some(feat) = candidates.get_mut(prev) {
                    self.remove_matched_peaks(&mut query, feat);
                    feat.rank = prev as u32 + 1;
                }
                prev = candidates.len()
            } else {
                break;
            }
        }
        candidates
    }

    /// Score a preliminary candidate: as it is, or, for a mass offset, with
    /// the offset on each site it may take (each a candidate of its own, as
    /// an indexed modification's forms would be)
    fn score_hypotheses(&self, query: &ProcessedSpectrum, pre: &PreScore, annotate: bool) -> Vec<Scored> {
        if pre.offset == 0 {
            let (score, fragments) = self.score_peptide(query, pre, &self.db[pre.peptide], annotate);
            return vec![(score, fragments, None)];
        }
        let offset = pre.offset - 1;
        let base = &self.db[pre.peptide];
        self.db
            .offset_sites(base, offset as usize)
            .into_iter()
            .map(|site| {
                let placement = OffsetPlacement { offset, site };
                let (score, fragments) = self.score_peptide(query, pre, &self.db.place(base, placement), annotate);
                (score, fragments, Some(placement))
            })
            .collect()
    }

    /// Calculate full hyperscore for a given PSM
    #[allow(dead_code)]
    fn score_candidate(
        &self,
        query: &ProcessedSpectrum,
        pre_score: &PreScore,
    ) -> (Score, Option<Fragments>) {
        self.score_peptide(query, pre_score, &self.db[pre_score.peptide], self.annotate_matches)
    }

    /// Hyperscore of any `peptide` (e.g. a modification moved to another
    /// residue) against `query`, as the search scores its candidates
    pub fn rescore(&self, query: &ProcessedSpectrum, peptide: &Peptide, precursor_charge: u8) -> f64 {
        let pre = PreScore { precursor_charge, ..Default::default() };
        self.score_peptide(query, &pre, peptide, false).0.hyperscore
    }

    fn score_peptide(
        &self,
        query: &ProcessedSpectrum,
        pre_score: &PreScore,
        peptide: &Peptide,
        annotate: bool,
    ) -> (Score, Option<Fragments>) {
        let mut score = Score {
            peptide: pre_score.peptide,
            precursor_charge: pre_score.precursor_charge,
            isotope_error: pre_score.isotope_error,
            ..Default::default()
        };
        let max_fragment_charge =
            max_fragment_charge(self.max_fragment_charge, score.precursor_charge);

        // Regenerate theoretical ions - initial database search might be
        // using only a subset of all possible ions (e.g. no b1/b2/y1/y2)
        // so we need to completely re-score this candidate
        let fragments =
            self.db.ion_kinds.iter().flat_map(|kind| {
                IonSeries::new(peptide, *kind, &self.db.modifications).enumerate()
            });

        let mut b_run = Run::default();
        let mut y_run = Run::default();

        let mut fragments_details = Fragments::default();

        // modifications whose fragments lose a neutral mass (phospho S/T)
        let lossy = match self.db.modifications.has_neutral_losses() {
            true => crate::ion_series::lossy_sites(peptide, &self.db.modifications),
            false => Vec::new(),
        };
        let mut matched: Vec<usize> = Vec::new();

        for (idx, frag) in fragments {
            for charge in 1..max_fragment_charge {
                // Experimental peaks are multipled by charge, therefore theoretical are divided
                let mz = frag.monoisotopic_mass / charge as f32;

                if let Some(peak_idx) = crate::spectrum::select_most_intense_peak(
                    &query.masses,
                    &query.intensities,
                    mz,
                    self.fragment_tol,
                    None,
                ) {
                    let peak_mass = query.masses[peak_idx];
                    let peak_intensity = query.intensities[peak_idx];
                    let fragment_charge = query.charges[peak_idx].max(charge);
                    if !lossy.is_empty() {
                        matched.push(peak_idx);
                    }

                    score.ppm_difference +=
                        peak_intensity * (mz - peak_mass).abs() * 2E6 / (mz + peak_mass);

                    let exp_mz = query.peak_mz(peak_idx);
                    let calc_mz = frag.monoisotopic_mass / fragment_charge as f32 + PROTON;

                    match frag.kind {
                        Kind::A | Kind::B | Kind::C => {
                            score.matched_b += 1;
                            score.summed_b += peak_intensity;
                            b_run.matched(idx);
                        }
                        Kind::X | Kind::Y | Kind::Z | Kind::ZRadical => {
                            score.matched_y += 1;
                            score.summed_y += peak_intensity;
                            y_run.matched(idx);
                        }
                    }

                    if annotate {
                        let idx = match frag.kind {
                            Kind::A | Kind::B | Kind::C => idx as i32 + 1,
                            Kind::X | Kind::Y | Kind::Z | Kind::ZRadical => {
                                peptide.sequence.len().saturating_sub(1) as i32 - idx as i32
                            }
                        };
                        fragments_details.kinds.push(frag.kind);
                        fragments_details.charges.push(fragment_charge as i32);
                        fragments_details.mz_experimental.push(exp_mz);
                        fragments_details.mz_calculated.push(calc_mz);
                        fragments_details.fragment_ordinals.push(idx);
                        fragments_details.intensities.push(peak_intensity);
                    }
                }
            }
        }

        // fragments containing a lossy site, less the loss: matched peaks
        // and intensity (not ion ladders, nor annotations); a peak counts once
        if !lossy.is_empty() {
            let fragments = self.db.ion_kinds.iter().flat_map(|kind| {
                IonSeries::new(peptide, *kind, &self.db.modifications).enumerate()
            });
            for (idx, frag) in fragments {
                for loss in crate::ion_series::fragment_losses(frag.kind, idx, &lossy) {
                    for charge in 1..max_fragment_charge {
                        let mz = (frag.monoisotopic_mass - loss) / charge as f32;
                        let Some(peak_idx) = crate::spectrum::select_most_intense_peak(
                            &query.masses,
                            &query.intensities,
                            mz,
                            self.fragment_tol,
                            None,
                        ) else {
                            continue;
                        };
                        if matched.contains(&peak_idx) {
                            continue;
                        }
                        matched.push(peak_idx);
                        let peak_mass = query.masses[peak_idx];
                        let peak_intensity = query.intensities[peak_idx];
                        score.ppm_difference += peak_intensity * (mz - peak_mass).abs() * 2E6 / (mz + peak_mass);
                        match frag.kind {
                            Kind::A | Kind::B | Kind::C => {
                                score.matched_b += 1;
                                score.summed_b += peak_intensity;
                            }
                            Kind::X | Kind::Y | Kind::Z | Kind::ZRadical => {
                                score.matched_y += 1;
                                score.summed_y += peak_intensity;
                            }
                        }
                    }
                }
            }
        }

        score.hyperscore = score.hyperscore(self.score_type);
        score.longest_b = b_run.longest;
        score.longest_y = y_run.longest;
        score.ppm_difference /= score.summed_b + score.summed_y;

        if annotate {
            (score, Some(fragments_details))
        } else {
            // drop(fragments_details);
            (score, None)
        }
    }
}

/// Maintain information about the longest continous ion ladder for a series
#[derive(Default)]
struct Run {
    start: usize,
    length: usize,
    last: usize,
    pub longest: usize,
}

impl Run {
    pub fn matched(&mut self, index: usize) {
        if self.last == index {
            return;
        } else if self.start + self.length == index {
            self.length += 1;
            self.longest = self.longest.max(self.length);
        } else {
            self.start = index;
            self.length = 1;
            self.longest = self.longest.max(self.length);
        }
        self.last = index;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn longest_series() {
        let mut run = Run::default();

        run.matched(1);
        run.matched(2);
        run.matched(3);
        run.matched(3);
        run.matched(3);

        assert_eq!(run.length, 3);
        assert_eq!(run.longest, 3);

        run.matched(5);
        run.matched(5);
        assert_eq!(run.length, 1);
        assert_eq!(run.longest, 3);
        run.matched(6);
        assert_eq!(run.length, 2);
    }

    #[test]
    fn test_max_fragment_charge() {
        assert_eq!(max_fragment_charge(None, 1), 2);
        assert_eq!(max_fragment_charge(None, 2), 2);
        assert_eq!(max_fragment_charge(None, 3), 3);
        assert_eq!(max_fragment_charge(None, 4), 4);
        assert_eq!(max_fragment_charge(Some(1), 2), 2);
        assert_eq!(max_fragment_charge(Some(1), 3), 2);
        assert_eq!(max_fragment_charge(Some(2), 4), 3);
        assert_eq!(max_fragment_charge(Some(4), 1), 2);
    }

    #[test]
    fn score_ordering_issue_236() {
        let mut scores = [(1, 100.0), (2, 50.0), (3, 1.0)]
            .into_iter()
            .map(|(index, hyperscore)| Score {
                peptide: PeptideIx(index),
                hyperscore,
                ..Default::default()
            })
            .collect::<Vec<_>>();
        bounded_min_heapify(&mut scores, 2);
        let mut retained = scores[..2]
            .iter()
            .map(|score| score.hyperscore)
            .collect::<Vec<_>>();
        retained.sort_by(f64::total_cmp);
        assert_eq!(retained, vec![50.0, 100.0]);
    }
}

#[cfg(test)]
mod mass_offset_tests {
    use super::*;
    use crate::database::{Builder, OffsetPlacement, Parameters};
    use crate::peptide::Site;

    fn parameters(variable: serde_json::Value) -> Result<Parameters, crate::modification::ModificationError> {
        let builder: Builder = serde_json::from_value(serde_json::json!({
            "fasta": "none",
            "enzyme": {"missed_cleavages": 1, "min_len": 5},
            "static_mods": {"C": 57.021465},
            "variable_mods": variable,
        }))
        .unwrap();
        builder.make_parameters()
    }

    fn dsso() -> serde_json::Value {
        serde_json::json!({
            "M": [15.9949],
            "K": [{"mass": 176.01433, "name": "DSSO-OH", "search_mode": "mass_offset", "position": "not_c_term"}],
            "[": [{"mass": 176.01433, "name": "DSSO-OH", "search_mode": "mass_offset"}]
        })
    }

    fn scorer(db: &IndexedDatabase) -> Scorer<'_> {
        Scorer {
            db,
            precursor_tol: Tolerance::Ppm(-10.0, 10.0),
            fragment_tol: Tolerance::Ppm(-10.0, 10.0),
            min_matched_peaks: 4,
            min_isotope_err: 0,
            max_isotope_err: 0,
            min_precursor_charge: 2,
            max_precursor_charge: 3,
            override_precursor_charge: false,
            max_fragment_charge: Some(1),
            chimera: false,
            report_psms: 2,
            wide_window: false,
            annotate_matches: false,
            score_type: ScoreType::SageHyperScore,
        }
    }

    /// A doubly charged spectrum of `peptide`'s b and y ions, with noise
    fn spectrum(db: &IndexedDatabase, peptide: &Peptide) -> ProcessedSpectrum {
        let mut peaks: Vec<(f32, f32)> = [Kind::B, Kind::Y]
            .iter()
            .flat_map(|kind| IonSeries::new(peptide, *kind, &db.modifications))
            .map(|ion| (ion.monoisotopic_mass, 100.0))
            .collect();
        peaks.extend((0..40).map(|k| (150.31 + k as f32 * 23.137, 10.0)));
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
        ProcessedSpectrum {
            level: 2,
            id: "scan=1".into(),
            precursors: vec![Precursor {
                mz: peptide.monoisotopic / 2.0 + PROTON,
                charge: Some(2),
                ..Default::default()
            }],
            total_ion_current: peaks.iter().map(|p| p.1).sum(),
            masses: peaks.iter().map(|p| p.0).collect(),
            intensities: peaks.iter().map(|p| p.1).collect(),
            charges: vec![1; peaks.len()],
            ..Default::default()
        }
    }

    #[test]
    fn configuration() {
        let params = parameters(dsso()).unwrap();
        let mods = &params.modifications;
        assert_eq!(mods.variable_mods.len(), 1, "only M is indexed");
        assert_eq!(mods.offset_mods.len(), 2);
        // written back as configured
        let json = serde_json::to_value(&params).unwrap();
        assert_eq!(json["variable_mods"]["K"][0]["search_mode"], "mass_offset");
        assert_eq!(json["variable_mods"]["K"][0]["position"], "not_c_term");
        assert_eq!(json["variable_mods"]["["][0]["search_mode"], "mass_offset");
        let builder: Builder = serde_json::from_value(serde_json::json!({
            "fasta": "none",
            "static_mods": json["static_mods"],
            "variable_mods": json["variable_mods"],
        }))
        .unwrap();
        assert_eq!(builder.make_parameters().unwrap().modifications.offset_mods.len(), 2);

        for bad in [
            serde_json::json!({"N": [{"mass": 1.0, "search_mode": "mass_offset", "motif": "N-x-[ST]"}]}),
            serde_json::json!({"K": [{"mass": 0.0, "search_mode": "mass_offset"}]}),
            serde_json::json!({"K": [{"mass": 42.0106, "search_mode": "mass_offset"}], "^": [42.0106]}),
            serde_json::json!({"K": [{"mass": 42.0106, "search_mode": "mass_offset", "site_mode": "library"}]}),
        ] {
            assert!(parameters(bad.clone()).is_err(), "{bad}");
        }
    }

    /// An offset is not indexed; the search finds it on the right lysine,
    /// and materializing gives the placed peptide an index of its own
    #[test]
    fn offsets_are_searched_and_materialized() {
        let params = parameters(dsso()).unwrap();
        let fasta = crate::fasta::Fasta::parse(">p\nMAGEKVLDSGAPIKAYRWNNPEDLLGFR\n".into(), "rev_", true);
        let mut db = params.build(fasta);
        assert_eq!(db.mass_offsets.len(), 1, "one offset with two kinds of site");
        assert_eq!(db.mass_offsets[0].sites.len(), 2);
        let offset = db.mass_offsets[0].id;
        assert!(db.peptides.iter().all(|p| !p.modifications.contains(&Some(offset)) && p.nterm != Some(offset)));

        let base = db.peptides.iter().position(|p| !p.decoy && p.sequence.as_ref() == b"VLDSGAPIKAYR").unwrap();
        let base = PeptideIx(base as u32);
        // K8, not the C-terminal residue (there is none here), nor the protein N-terminus
        assert_eq!(db.offset_sites(&db[base], 0), vec![Site::Sequence(8)]);
        // MAGEK: the protein N-terminus only (its K is C-terminal)
        let first = db.peptides.iter().find(|p| !p.decoy && p.sequence.as_ref() == b"MAGEK").unwrap();
        assert_eq!(db.offset_sites(first, 0), vec![Site::Nterm]);
        // a decoy has as many sites as its target
        let reversed = db[base].reverse();
        assert_eq!(reversed.sequence.as_ref(), b"VYAKIPAGSDLR");
        let decoy = db.peptides.iter().find(|p| p.decoy && p.sequence == reversed.sequence).expect("decoy");
        assert_eq!(db.offset_sites(decoy, 0).len(), 1);

        let placement = OffsetPlacement { offset: 0, site: Site::Sequence(8) };
        let placed = db.place(&db[base], placement);
        assert!((placed.monoisotopic - db[base].monoisotopic - 176.01433).abs() < 1e-3);
        let query = spectrum(&db, &placed);
        let mut features = {
            let scorer = scorer(&db);
            let features = scorer.score(&query);
            assert!(!features.is_empty());
            let best = &features[0];
            assert_eq!(best.mass_offset, Some(placement), "{features:?}");
            assert!(db.resolve(best).same_form(&placed));
            assert!(best.delta_mass.abs() < 1.0, "{}", best.delta_mass);
            assert!((best.calcmass - placed.monoisotopic).abs() < 1e-3);
            // chimeric search removes the placed peptide's peaks
            let chimeric = Scorer { chimera: true, ..scorer }.score(&query);
            assert_eq!(chimeric[0].mass_offset, Some(placement));
            features
        };
        // without the offset, nothing explains the precursor
        let plain = parameters(serde_json::json!({"M": [15.9949]})).unwrap();
        let plain_db = plain.build(crate::fasta::Fasta::parse(
            ">p\nMAGEKVLDSGAPIKAYRWNNPEDLLGFR\n".into(),
            "rev_",
            true,
        ));
        assert!(scorer(&plain_db).score(&query).iter().all(|f| f.delta_mass.abs() > 10.0));

        features.push(features[0].clone());
        let last = features.len() - 1;
        let n = db.peptides.len();
        let with_offsets = features.iter().filter(|f| f.mass_offset.is_some()).count();
        assert_eq!(db.materialize_mass_offsets(&mut features), with_offsets);
        assert!(features[0].peptide_idx.0 as usize >= n);
        assert_eq!(features[0].peptide_idx, features[last].peptide_idx, "one index per placed peptide");
        assert!(features[0].mass_offset.is_none());
        assert!(db[features[0].peptide_idx].same_form(&placed));
        assert_eq!(db.offset_peptides.len(), with_offsets - 1);
        assert_eq!(
            db[features[0].peptide_idx].display(&db.modifications).to_string(),
            "VLDSGAPIK[DSSO-OH]AYR"
        );
    }
}
