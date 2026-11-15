use flate2::{Decompress, FlushDecompress, Status};
use quick_xml::events::Event;
use quick_xml::Reader;
use rayon::prelude::*;
use sage_core::spectrum::{Precursor, Representation};
use sage_core::{mass::Tolerance, spectrum::RawSpectrum};
use tokio::io::AsyncBufRead;

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
/// Which tag are we inside?
enum State {
    Spectrum,
    Scan,
    BinaryDataArray,
    Binary,
    Precursor,
    SelectedIon,
}

#[derive(Copy, Clone, Debug)]
enum BinaryKind {
    Intensity,
    Mz,
    Noise,
}

#[derive(Copy, Clone, Debug)]
enum Dtype {
    F32,
    F64,
}

/// Base64 held before the queued spectra are decoded in parallel. A 1.2 GB
/// Thermo mzML is about 600 MB of base64, so a batch is a small part of a file.
/// It bounds the queue: each worker thread also keeps a base64 and an inflate
/// buffer as large as the largest array it has seen.
const DECODE_BATCH_BYTES: usize = 32 * 1024 * 1024;

/// One binary data array, still base64-encoded
struct PendingArray {
    kind: BinaryKind,
    dtype: Dtype,
    compression: bool,
    base64: Vec<u8>,
}

/// A parsed spectrum whose arrays are not decoded yet
struct PendingSpectrum {
    spectrum: RawSpectrum,
    arrays: Vec<PendingArray>,
}

/// Little-endian floats as `f32`
fn decode_float_array(bytes: &[u8], dtype: Dtype) -> Vec<f32> {
    match dtype {
        Dtype::F32 => bytes
            .chunks_exact(4)
            .map(|chunk| f32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]))
            .collect(),
        Dtype::F64 => bytes
            .chunks_exact(8)
            .map(|chunk| {
                let mut buf = [0u8; 8];
                buf.copy_from_slice(chunk);
                f64::from_le_bytes(buf) as f32
            })
            .collect(),
    }
}

/// Inflate one zlib-compressed binary array into `output`, reusing the
/// decompressor and the buffer across arrays (rather than a decoder, a future
/// and an allocation per array, of which a large mzML holds ~300,000)
fn inflate_into(inflater: &mut Decompress, input: &[u8], output: &mut Vec<u8>) -> Result<(), MzMLError> {
    fn invalid(error: impl std::fmt::Display) -> MzMLError {
        MzMLError::IOError(std::io::Error::new(std::io::ErrorKind::InvalidData, error.to_string()))
    }
    inflater.reset(true);
    output.clear();
    // arrays inflate to about 1.7x; ask for twice that, grow when a stream beats it
    output.reserve(input.len().saturating_mul(2).max(64));
    loop {
        let consumed = inflater.total_in() as usize;
        let before = inflater.total_out();
        let status = inflater
            .decompress_vec(&input[consumed..], output, FlushDecompress::Finish)
            .map_err(invalid)?;
        match status {
            Status::StreamEnd => return Ok(()),
            // `decompress_vec` writes only into spare capacity: no progress
            // with room left means the stream is cut short
            Status::Ok | Status::BufError => {
                if inflater.total_out() == before && inflater.total_in() as usize == consumed {
                    if output.len() < output.capacity() {
                        return Err(invalid("truncated zlib stream in a binary array"));
                    }
                    output.reserve(output.capacity().max(64));
                }
            }
        }
    }
}

/// Decode a batch of spectra in parallel (base64, inflate, floats) and append
/// them in file order. The XML has to be read in order; decoding, most of the
/// cost of reading an mzML, need not. With S/N on, intensities at that MS
/// level are divided by the spectrum's noise array. The first bad array in
/// file order is the error reported, as when decoding in order.
fn decode_batch(batch: Vec<PendingSpectrum>, signal_to_noise: Option<u8>, spectra: &mut Vec<RawSpectrum>) -> Result<(), MzMLError> {
    if batch.is_empty() {
        return Ok(());
    }
    let decoded: Vec<Result<RawSpectrum, MzMLError>> = batch
        .into_par_iter()
        .map_init(
            || (Vec::<u8>::new(), Vec::<u8>::new(), Decompress::new(true)),
            |(base64_buffer, inflate_buffer, inflater), pending| {
                let PendingSpectrum { mut spectrum, arrays } = pending;
                let mut noise = Vec::new();
                for array in &arrays {
                    base64_buffer.clear();
                    base64::decode_config_buf(&array.base64, base64::STANDARD, base64_buffer)?;
                    let bytes = match array.compression {
                        false => &base64_buffer[..],
                        true => {
                            inflate_into(inflater, base64_buffer, inflate_buffer)?;
                            &inflate_buffer[..]
                        }
                    };
                    let values = decode_float_array(bytes, array.dtype);
                    match array.kind {
                        BinaryKind::Intensity => spectrum.intensity = values,
                        BinaryKind::Mz => spectrum.mz = values,
                        BinaryKind::Noise => noise = values,
                    }
                }
                if signal_to_noise == Some(spectrum.ms_level) && !noise.is_empty() {
                    spectrum
                        .intensity
                        .iter_mut()
                        .zip(noise.iter())
                        .for_each(|(int, noise)| *int /= noise);
                }
                Ok(spectrum)
            },
        )
        .collect();
    for spectrum in decoded {
        spectra.push(spectrum?);
    }
    Ok(())
}

// MUST supply only one of the following
const ZLIB_COMPRESSION: &[u8] = b"MS:1000574";
const NO_COMPRESSION: &[u8] = b"MS:1000576";

// Numpress encodings, alone or followed by zlib. There is no numpress decoder:
// such an array is refused by name, as its bytes are not IEEE floats and read
// as floats would give wrong numbers without a word.
const NUMPRESS: [(&[u8], &str); 6] = [
    (b"MS:1002312", "MS-Numpress linear prediction compression"),
    (b"MS:1002313", "MS-Numpress positive integer compression"),
    (b"MS:1002314", "MS-Numpress short logged float compression"),
    (b"MS:1002746", "MS-Numpress linear prediction compression followed by zlib"),
    (b"MS:1002747", "MS-Numpress positive integer compression followed by zlib"),
    (b"MS:1002748", "MS-Numpress short logged float compression followed by zlib"),
];

// MUST supply only one of the following
const INTENSITY_ARRAY: &[u8] = b"MS:1000515";
const MZ_ARRAY: &[u8] = b"MS:1000514";
const NOISE_ARRAY: &[u8] = b"MS:1002744";

// MUST supply only one of the following
const FLOAT_64: &[u8] = b"MS:1000523";
const FLOAT_32: &[u8] = b"MS:1000521";

const MS_LEVEL: &[u8] = b"MS:1000511";
const PROFILE: &[u8] = b"MS:1000128";
const CENTROID: &[u8] = b"MS:1000127";
const TOTAL_ION_CURRENT: &[u8] = b"MS:1000285";

const SCAN_START_TIME: &[u8] = b"MS:1000016";
const UNIT_SECONDS: &[u8] = b"UO:0000010";
const UNIT_MINUTES: &[u8] = b"UO:0000031";
const ION_INJECTION_TIME: &[u8] = b"MS:1000927";

const SELECTED_ION_MZ: &[u8] = b"MS:1000744";
const SELECTED_ION_INT: &[u8] = b"MS:1000042";
const SELECTED_ION_CHARGE: &[u8] = b"MS:1000041";

const ISO_WINDOW_TARGET: &[u8] = b"MS:1000827";
const ISO_WINDOW_LOWER: &[u8] = b"MS:1000828";
const ISO_WINDOW_UPPER: &[u8] = b"MS:1000829";

const INVERSE_ION_MOBILITY: &[u8] = b"MS:1002815";

pub struct MzMLReader {
    ms_level: Option<u8>,
    // If set to Some(level) and noise intensities are present in the MzML file,
    // divide intensities at this MS-level by noise to calculate S/N
    signal_to_noise: Option<u8>,

    file_id: usize,
    /// Base64 per decode batch: [`DECODE_BATCH_BYTES`], lowered only by tests
    decode_batch_bytes: usize,
}

impl MzMLReader {
    /// Create a new [`MzMlReader`] with a minimum MS level filter
    ///
    /// # Example
    ///
    /// A minimum level of 2 will not parse or return MS1 scans
    pub fn with_file_id_and_level_filter(file_id: usize, ms_level: u8) -> Self {
        Self {
            ms_level: Some(ms_level),
            file_id,
            signal_to_noise: None,
            decode_batch_bytes: DECODE_BATCH_BYTES,
        }
    }

    pub fn with_file_id(file_id: usize) -> Self {
        Self {
            ms_level: None,
            signal_to_noise: None,
            file_id,
            decode_batch_bytes: DECODE_BATCH_BYTES,
        }
    }

    /// Decode a batch once it holds this much base64 (tests: to cross batches
    /// without a 32 MB file)
    #[cfg(test)]
    fn set_decode_batch_bytes(&mut self, bytes: usize) -> &mut Self {
        self.decode_batch_bytes = bytes.max(1);
        self
    }

    pub fn set_file_id(&mut self, file_id: usize) -> &mut Self {
        self.file_id = file_id;
        self
    }

    pub fn set_signal_to_noise(&mut self, sn: Option<u8>) -> &mut Self {
        self.signal_to_noise = sn;
        self
    }

    /// Here be dragons -
    /// Seriously, this kinda sucks because it's a giant imperative, stateful loop.
    /// But I also don't want to spend any more time working on an mzML parser...
    pub async fn parse<B: AsyncBufRead + Unpin>(
        &self,
        b: B,
    ) -> Result<Vec<RawSpectrum>, MzMLError> {
        let mut reader = Reader::from_reader(b);
        let mut buf = Vec::new();

        let mut state = None;
        let mut compression = false;
        let mut binary_dtype = Dtype::F64;
        let mut binary_array = None;
        // the arrays of the spectrum being read, and the spectra read, still
        // encoded: decoded in parallel every `decode_batch_bytes` of base64
        let mut pending_arrays: Vec<PendingArray> = Vec::new();
        let mut pending: Vec<PendingSpectrum> = Vec::new();
        let mut pending_bytes = 0usize;

        let mut spectrum = RawSpectrum::default_with_file_id(self.file_id);
        let mut precursor = Precursor::default();
        let mut iso_window_lo: Option<f32> = None;
        let mut iso_window_hi: Option<f32> = None;
        let mut spectra = Vec::new();

        macro_rules! extract {
            ($ev:expr, $key:expr) => {
                $ev.try_get_attribute($key)?
                    .ok_or(MzMLError::Malformed)?
                    .value
            };
        }

        macro_rules! extract_value {
            ($ev:expr) => {{
                let s = $ev
                    .try_get_attribute(b"value")?
                    .ok_or(MzMLError::Malformed)?
                    .value;
                std::str::from_utf8(&s)?.parse()?
            }};
        }

        loop {
            match reader.read_event_into_async(&mut buf).await {
                Ok(Event::Start(ref ev)) => {
                    // State transition into child tag
                    state = match (ev.name().into_inner(), state) {
                        (b"spectrum", _) => Some(State::Spectrum),
                        (b"scan", Some(State::Spectrum)) => Some(State::Scan),
                        (b"binaryDataArray", Some(State::Spectrum)) => Some(State::BinaryDataArray),
                        (b"binary", Some(State::BinaryDataArray)) => Some(State::Binary),
                        (b"precursor", Some(State::Spectrum)) => Some(State::Precursor),
                        (b"selectedIon", Some(State::Precursor)) => Some(State::SelectedIon),
                        _ => state,
                    };
                    match ev.name().into_inner() {
                        b"spectrum" => {
                            let id = extract!(ev, b"id");
                            let id = std::str::from_utf8(&id)?;
                            spectrum.id = id.to_string();
                        }
                        b"precursor" => {
                            // Not all precursor fields have a spectrumRef
                            if let Some(scan) = ev.try_get_attribute(b"spectrumRef")? {
                                let scan = std::str::from_utf8(&scan.value)?;
                                precursor.spectrum_ref = Some(scan.to_string())
                            }
                        }
                        _ => {}
                    }
                }
                Ok(Event::Empty(ref ev)) => match (state, ev.name().into_inner()) {
                    (Some(State::BinaryDataArray), b"cvParam") => {
                        let accession = extract!(ev, b"accession");
                        match accession.as_ref() {
                            ZLIB_COMPRESSION => compression = true,
                            NO_COMPRESSION => compression = false,
                            FLOAT_64 => binary_dtype = Dtype::F64,
                            FLOAT_32 => binary_dtype = Dtype::F32,
                            INTENSITY_ARRAY => binary_array = Some(BinaryKind::Intensity),
                            MZ_ARRAY => binary_array = Some(BinaryKind::Mz),
                            NOISE_ARRAY => binary_array = Some(BinaryKind::Noise),
                            accession if NUMPRESS.iter().any(|(a, _)| *a == accession) => {
                                let name = NUMPRESS.iter().find(|(a, _)| *a == accession).map_or("", |(_, n)| n);
                                return Err(MzMLError::UnsupportedCV(format!(
                                    "{} ({name}): zlib-compressed or uncompressed binary arrays only",
                                    String::from_utf8_lossy(accession)
                                )));
                            }
                            _ => {
                                // Unknown CV - perhaps noise
                                binary_array = None;
                            }
                        }
                    }
                    (Some(State::Spectrum), b"cvParam") => {
                        let accession = extract!(ev, b"accession");
                        match accession.as_ref() {
                            MS_LEVEL => {
                                let level = extract_value!(ev);
                                if let Some(filter) = self.ms_level {
                                    if level != filter {
                                        spectrum = RawSpectrum::default_with_file_id(self.file_id);
                                        state = None;
                                    }
                                }
                                spectrum.ms_level = level;
                            }
                            PROFILE => spectrum.representation = Representation::Profile,
                            CENTROID => spectrum.representation = Representation::Centroid,
                            TOTAL_ION_CURRENT => {
                                let value = extract_value!(ev);
                                if value == 0.0 {
                                    // No ion current, break out of current state
                                    spectrum = RawSpectrum::default_with_file_id(self.file_id);
                                    state = None;
                                } else {
                                    spectrum.total_ion_current = value;
                                }
                            }
                            _ => {}
                        }
                    }
                    (Some(State::Precursor), b"cvParam") => {
                        let accession = extract!(ev, b"accession");
                        match accession.as_ref() {
                            ISO_WINDOW_TARGET => {
                                // use isolation window target for precursor m/z, e.g. to handle
                                // DIA setups where the mzML conversion software doesn't write
                                // a selection ion tag
                                if precursor.mz == 0.0 {
                                    precursor.mz = extract_value!(ev)
                                }
                            }
                            ISO_WINDOW_LOWER => iso_window_lo = Some(extract_value!(ev)),
                            ISO_WINDOW_UPPER => iso_window_hi = Some(extract_value!(ev)),
                            _ => {}
                        }
                    }
                    (Some(State::SelectedIon), b"cvParam") => {
                        let accession = extract!(ev, b"accession");
                        match accession.as_ref() {
                            SELECTED_ION_CHARGE => {
                                precursor.charge = Some(extract_value!(ev));
                            }
                            SELECTED_ION_MZ => {
                                let val = extract_value!(ev);
                                if val != 0.0 {
                                    precursor.mz = val;
                                }
                            }
                            SELECTED_ION_INT => {
                                precursor.intensity = Some(extract_value!(ev));
                            }
                            INVERSE_ION_MOBILITY => {
                                precursor.inverse_ion_mobility = Some(extract_value!(ev));
                            }
                            _ => {}
                        }
                    }
                    (Some(State::Scan), b"cvParam") => {
                        let accession = extract!(ev, b"accession");
                        match accession.as_ref() {
                            SCAN_START_TIME => {
                                let scan_start_time = extract_value!(ev);
                                let unit = extract!(ev, b"unitAccession");

                                spectrum.scan_start_time = match unit.as_ref() {
                                    UNIT_SECONDS => scan_start_time / 60.0,
                                    UNIT_MINUTES => scan_start_time,
                                    _ => return Err(MzMLError::Malformed),
                                };
                            }
                            ION_INJECTION_TIME => {
                                spectrum.ion_injection_time = extract_value!(ev);
                            }
                            INVERSE_ION_MOBILITY => {
                                precursor.inverse_ion_mobility = Some(extract_value!(ev));
                            }
                            _ => {}
                        }
                    }

                    _ => {}
                },
                Ok(Event::Text(text)) => {
                    if let Some(State::Binary) = state {
                        if let Some(filter) = self.ms_level {
                            if spectrum.ms_level != filter {
                                continue;
                            }
                        }
                        let raw = text.unescape()?;
                        // There are occasionally empty binary data arrays, or unknown CVs
                        if raw.is_empty() || binary_array.is_none() {
                            continue;
                        }
                        // kept encoded, decoded with the batch (`decode_batch`)
                        if let Some(kind) = binary_array.take() {
                            pending_bytes += raw.len();
                            pending_arrays.push(PendingArray {
                                kind,
                                dtype: binary_dtype,
                                compression,
                                base64: raw.into_owned().into_bytes(),
                            });
                        }
                    }
                }
                Ok(Event::End(ev)) => {
                    state = match (state, ev.name().into_inner()) {
                        (Some(State::Binary), b"binary") => Some(State::BinaryDataArray),
                        (Some(State::BinaryDataArray), b"binaryDataArray") => Some(State::Spectrum),
                        (Some(State::SelectedIon), b"selectedIon") => Some(State::Precursor),
                        (Some(State::Precursor), b"precursor") => {
                            if precursor.mz != 0.0 {
                                precursor.isolation_window = match (iso_window_lo, iso_window_hi) {
                                    (Some(lo), Some(hi)) => Some(Tolerance::Da(-lo, hi)),
                                    _ => None,
                                };
                                spectrum.precursors.push(precursor);
                                precursor = Precursor::default();
                            }
                            Some(State::Spectrum)
                        }
                        (Some(State::Scan), b"scan") => Some(State::Spectrum),
                        (_, b"spectrum") => {
                            let allow = self
                                .ms_level
                                .as_ref()
                                .map(|&level| level == spectrum.ms_level)
                                .unwrap_or(true);

                            if allow {
                                pending.push(PendingSpectrum { spectrum, arrays: std::mem::take(&mut pending_arrays) });
                            } else {
                                pending_arrays.clear();
                            }
                            if pending_bytes >= self.decode_batch_bytes {
                                decode_batch(std::mem::take(&mut pending), self.signal_to_noise, &mut spectra)?;
                                pending_bytes = 0;
                            }
                            spectrum = RawSpectrum::default_with_file_id(self.file_id);
                            None
                        }
                        _ => state,
                    };
                }
                Ok(Event::Eof) => break,
                Ok(_) => {}
                Err(err) => {
                    log::error!("unhandled XML error while parsing mzML: {}", err)
                }
            }
            buf.clear();
        }
        decode_batch(pending, self.signal_to_noise, &mut spectra)?;
        Ok(spectra)
    }
}

#[derive(thiserror::Error, Debug)]
pub enum MzMLError {
    #[error("malformed MzML")]
    Malformed,
    #[error("unsupported cvParam {0}")]
    UnsupportedCV(String),
    #[error("XML parsing error: {0}")]
    XMLError(#[from] quick_xml::Error),
    #[error("io error: {0}")]
    IOError(#[from] std::io::Error),
    #[error("utf8 error: {0}")]
    Utf8Error(#[from] std::str::Utf8Error),
    #[error("error parsing float: {0}")]
    FloatError(#[from] std::num::ParseFloatError),
    #[error("error parsing int: {0}")]
    IntError(#[from] std::num::ParseIntError),
    #[error("error decoding base64: {0}")]
    Base64Error(#[from] base64::DecodeError),
}

#[cfg(test)]
mod test {
    use sage_core::{mass::Tolerance, spectrum::Representation};

    use super::{MzMLError, MzMLReader};

    #[tokio::test]
    async fn parse_spectrum_issue_78() -> Result<(), MzMLError> {
        let s = r#"
        <spectrum id="spectrum=2442" index="286" defaultArrayLength="102" dataProcessingRef="dp_sp_1">
            <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" />
            <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2" />
            <cvParam cvRef="MS" accession="MS:1000294" name="mass spectrum" />
            <cvParam cvRef="MS" accession="MS:1000130" name="positive scan" />
            <cvParam cvRef="MS" accession="MS:1000504" name="base peak m/z" value="638.352905273437955"/>
            <cvParam cvRef="MS" accession="MS:1000505" name="base peak intensity" value="113.885513305664006"/>
            <cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="793.395202636718977"/>
            <cvParam cvRef="MS" accession="MS:1000528" name="lowest observed m/z" value="147.290603637695"/>
            <cvParam cvRef="MS" accession="MS:1000527" name="highest observed m/z" value="769.255798339843977"/>
            <userParam name="filter string" type="xsd:string" value="ITMS + c NSI d w Full ms2 457.72@cid35.00 [115.00-930.00]"/>
            <userParam name="preset scan configuration" type="xsd:string" value="2"/>
            <scanList count="1">
                <cvParam cvRef="MS" accession="MS:1000795" name="no combination" />
                <scan >
                    <cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="1503.96166992188" unitAccession="UO:0000010" unitName="second" unitCvRef="UO" />
                    <userParam name="[Thermo Trailer Extra]Monoisotopic M/Z:" type="xsd:double" value="457.723968505858977"/>
                    <scanWindowList count="1">
                        <scanWindow>
                            <cvParam cvRef="MS" accession="MS:1000501" name="scan window lower limit" value="115" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                            <cvParam cvRef="MS" accession="MS:1000500" name="scan window upper limit" value="930" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        </scanWindow>
                    </scanWindowList>
                </scan>
            </scanList>
            <precursorList count="1">
                <precursor>
                    <isolationWindow>
                        <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="457.723968505859" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="1.5" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="0.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                    </isolationWindow>
                    <selectedIonList count="1">
                        <selectedIon>
                            <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="457.723968505859" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                            <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="2" />
                            <cvParam cvRef="MS" accession="MS:1002815" name="inverse reduced ion mobility" value="1.078628" unitAccession="MS:1002814" unitName="volt-second per square centimeter"/>
                        </selectedIon>
                    </selectedIonList>
                    <activation>
                        <cvParam cvRef="MS" accession="MS:1000133" name="collision-induced dissociation" />
                        <cvParam cvRef="MS" accession="MS:1000045" name="collision energy" value="35.0"/>
                    </activation>
                </precursor>
            </precursorList>
            <binaryDataArrayList count="2">
                <binaryDataArray encodedLength="1088">
                    <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                    <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" />
                    <cvParam cvRef="MS" accession="MS:1000576" name="no compression" />
                    <binary>AAAAoExpYkAAAACA3MpkQAAAAACph2VAAAAAAE4wZkAAAACAlMdmQAAAAECZAmdAAAAAwP9jaEAAAADgj4ZoQAAAAGC7HWlAAAAAAOXFaUAAAADg+4dqQAAAAMC1pmpAAAAA4IGFa0AAAACAaUZsQAAAACBzYW1AAAAAANCjbUAAAACAQ6duQAAAAIDsxG5AAAAAQKIlb0AAAACA5z9vQAAAAIDuw29AAAAAAJQicEAAAAAg9UZwQAAAAKCeVHBAAAAAIInEcEAAAACAcs5wQAAAAOA6BHFAAAAAADoOcUAAAAAgfcRxQAAAAOA68nFAAAAAoPExckAAAADATKVyQAAAAMC10nJAAAAAwBJHc0AAAAAA7FNzQAAAAIAYkXNAAAAAgJzRc0AAAABgE2R0QAAAAMCrc3RAAAAAgE+zdEAAAAAAhMR0QAAAAIC64XRAAAAA4Cf/dEAAAADgy3B1QAAAAMCVgnVAAAAAoDugdUAAAACAX/Z1QAAAAAAAB3ZAAAAAgO4XdkAAAABAqEJ2QAAAAIDp8nZAAAAAIAgRd0AAAACggzR3QAAAAODwT3dAAAAAIHJsd0AAAAAA4YJ3QAAAAGC91ndAAAAAAL3id0AAAADg0xZ4QAAAAOA5NXhAAAAAYDaPeEAAAACgK7p4QAAAACCm0XhAAAAA4GHkeEAAAADgyPJ4QAAAAOB5/3hAAAAAoFtNeUAAAADA8H15QAAAAGAHtXlAAAAAoD7HeUAAAAAAEtR5QAAAAGCx5XlAAAAA4NEJekAAAAAgtVN6QAAAACDCX3pAAAAAIAqmekAAAACg4OR6QAAAAGDymnxAAAAAICV/fUAAAAAgd6Z9QAAAAKDYA4BAAAAAoCoVgEAAAACA/kOAQAAAAKCpYoBAAAAA4MycgEAAAADA3DyBQAAAAKCbrIFAAAAAoPC6gUAAAADgV22CQAAAACABY4NAAAAAQE+qg0AAAADA0vKDQAAAAEDz+oNAAAAAoIxrhEAAAADg6euEQAAAAIAuDIVAAAAAoOwjhUAAAACgZUuFQAAAAADdm4VAAAAAoCzrh0AAAABgYvWHQAAAAOALCohA</binary>
                </binaryDataArray>
                <binaryDataArray encodedLength="544">
                    <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" unitAccession="MS:1000131" unitName="number of detector counts" unitCvRef="MS"/>
                    <cvParam cvRef="MS" accession="MS:1000521" name="32-bit float" />
                    <cvParam cvRef="MS" accession="MS:1000576" name="no compression" />
                    <binary>3FlbQDg/ZUB8w3FAV2fMQMiOnkCXfP4/T2I2QC6qskAnhOZA/NU2QCc2QEAI1UhAQcAbQRrziUBmHq5AXutSQWZDbkAZGWdAzt6lQYNptUDSFDNBoY4IQAYaQEDeT7Q/16HGP9GtXUCITrQ/Rxu0Pzhc6j9mpjZAX1X8P7tPQ0AqxS5BZTzZPye+m0B7Sa5AfPsPQRr/W0CYwBRBwDh3QMAmtD/nq6E/bJHGPxJ9UUDsy/dAoCYMQRM2a0BkAR9Boo5pQMV0VEArYu5A4kaMQAyTI0BQPRJAML3TQCKVCED85+tArObGP1BVP0EtJuVAdyKAQFjctkFQa2NBixMTQXyyjUFX8eo/IHelQTdFcEFo1zZAhagsQAO53EBIugRB0M+gQfhBgkH0MsJAbGlIQZXg+EHe6CZBsbA2QHMHOECtW6BAjE2oQUpZckBasZ1AtKl3QEZYIUHkip1AQX7TQPqF60GNuaE/USk2QGLF40Im65ZAmXqlQBGuSUC70KBAAneMQeK3aEB87MVA5NigQE/Wb0BO475A</binary>
                </binaryDataArray>
            </binaryDataArrayList>
        </spectrum>
        "#;
        let mut spectra = MzMLReader::with_file_id(0).parse(s.as_bytes()).await?;

        assert_eq!(spectra.len(), 1);
        let s = spectra.pop().unwrap();

        assert_eq!(s.id, "spectrum=2442");
        assert_eq!(s.ms_level, 2);
        assert_eq!(s.representation, Representation::Centroid);
        assert_eq!(s.precursors.len(), 1);
        assert_eq!(s.precursors[0].charge, Some(2));
        assert!((s.precursors[0].mz - 457.723968) < 0.0001);
        assert!(match s.precursors[0].inverse_ion_mobility {
            Some(x) => (x - 1.0786) < 0.0001,
            None => false,
        });
        assert_eq!(
            s.precursors[0].isolation_window,
            Some(Tolerance::Da(-1.5, 0.75))
        );
        assert!((s.scan_start_time - 25.066).abs() < 0.0001);
        assert_eq!(s.ion_injection_time, 0.0);
        assert_eq!(s.intensity.len(), s.mz.len());
        Ok(())
    }

    #[tokio::test]
    async fn parse_spectrum_issue_117() -> Result<(), MzMLError> {
        // The issue was that some converters write the ion mobility as part of the selected ion (as in the last test)
        // and some write it as part of the scan, as in this test. This test checks that it can be read
        // fom the scan section.
        let s = r#"
        <spectrum id="spectrum=8678309" index="8678309" defaultArrayLength="102" dataProcessingRef="dp_sp_1">
            <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" />
            <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2" />
            <cvParam cvRef="MS" accession="MS:1000294" name="mass spectrum" />
            <cvParam cvRef="MS" accession="MS:1000130" name="positive scan" />
            <cvParam cvRef="MS" accession="MS:1000504" name="base peak m/z" value="638.352905273437955"/>
            <cvParam cvRef="MS" accession="MS:1000505" name="base peak intensity" value="113.885513305664006"/>
            <cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="793.395202636718977"/>
            <userParam name="filter string" type="xsd:string" value="ITMS + c NSI d w Full ms2 457.72@cid35.00 [115.00-930.00]"/>
            <scanList count="1">
                <cvParam cvRef="MS" accession="MS:1000795" name="no combination" />
                <scan >
                    <cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="1503.96166992188" unitAccession="UO:0000010" unitName="second" unitCvRef="UO" />
                    <cvParam cvRef="MS" accession="MS:1002815" name="inverse reduced ion mobility" value="1.078628" unitAccession="MS:1002814" unitName="volt-second per square centimeter"/>
                    <userParam name="[Thermo Trailer Extra]Monoisotopic M/Z:" type="xsd:double" value="457.723968505858977"/>
                    <scanWindowList count="1">
                        <scanWindow>
                            <cvParam cvRef="MS" accession="MS:1000501" name="scan window lower limit" value="115" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                            <cvParam cvRef="MS" accession="MS:1000500" name="scan window upper limit" value="930" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        </scanWindow>
                    </scanWindowList>
                </scan>
            </scanList>
            <precursorList count="1">
                <precursor>
                    <isolationWindow>
                        <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="457.723968505859" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="1.5" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="0.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                    </isolationWindow>
                    <selectedIonList count="1">
                        <selectedIon>
                            <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="457.723968505859" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                            <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="2" />
                        </selectedIon>
                    </selectedIonList>
                    <activation>
                        <cvParam cvRef="MS" accession="MS:1000133" name="collision-induced dissociation" />
                        <cvParam cvRef="MS" accession="MS:1000045" name="collision energy" value="35.0"/>
                    </activation>
                </precursor>
            </precursorList>
            <binaryDataArrayList count="2">
                <binaryDataArray encodedLength="1088">
                    <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                    <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" />
                    <cvParam cvRef="MS" accession="MS:1000576" name="no compression" />
                    <binary>AAAAoExpYkAAAACA3MpkQAAAAACph2VAAAAAAE4wZkAAAACAlMdmQAAAAECZAmdAAAAAwP9jaEAAAADgj4ZoQAAAAGC7HWlAAAAAAOXFaUAAAADg+4dqQAAAAMC1pmpAAAAA4IGFa0AAAACAaUZsQAAAACBzYW1AAAAAANCjbUAAAACAQ6duQAAAAIDsxG5AAAAAQKIlb0AAAACA5z9vQAAAAIDuw29AAAAAAJQicEAAAAAg9UZwQAAAAKCeVHBAAAAAIInEcEAAAACAcs5wQAAAAOA6BHFAAAAAADoOcUAAAAAgfcRxQAAAAOA68nFAAAAAoPExckAAAADATKVyQAAAAMC10nJAAAAAwBJHc0AAAAAA7FNzQAAAAIAYkXNAAAAAgJzRc0AAAABgE2R0QAAAAMCrc3RAAAAAgE+zdEAAAAAAhMR0QAAAAIC64XRAAAAA4Cf/dEAAAADgy3B1QAAAAMCVgnVAAAAAoDugdUAAAACAX/Z1QAAAAAAAB3ZAAAAAgO4XdkAAAABAqEJ2QAAAAIDp8nZAAAAAIAgRd0AAAACggzR3QAAAAODwT3dAAAAAIHJsd0AAAAAA4YJ3QAAAAGC91ndAAAAAAL3id0AAAADg0xZ4QAAAAOA5NXhAAAAAYDaPeEAAAACgK7p4QAAAACCm0XhAAAAA4GHkeEAAAADgyPJ4QAAAAOB5/3hAAAAAoFtNeUAAAADA8H15QAAAAGAHtXlAAAAAoD7HeUAAAAAAEtR5QAAAAGCx5XlAAAAA4NEJekAAAAAgtVN6QAAAACDCX3pAAAAAIAqmekAAAACg4OR6QAAAAGDymnxAAAAAICV/fUAAAAAgd6Z9QAAAAKDYA4BAAAAAoCoVgEAAAACA/kOAQAAAAKCpYoBAAAAA4MycgEAAAADA3DyBQAAAAKCbrIFAAAAAoPC6gUAAAADgV22CQAAAACABY4NAAAAAQE+qg0AAAADA0vKDQAAAAEDz+oNAAAAAoIxrhEAAAADg6euEQAAAAIAuDIVAAAAAoOwjhUAAAACgZUuFQAAAAADdm4VAAAAAoCzrh0AAAABgYvWHQAAAAOALCohA</binary>
                </binaryDataArray>
                <binaryDataArray encodedLength="544">
                    <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" unitAccession="MS:1000131" unitName="number of detector counts" unitCvRef="MS"/>
                    <cvParam cvRef="MS" accession="MS:1000521" name="32-bit float" />
                    <cvParam cvRef="MS" accession="MS:1000576" name="no compression" />
                    <binary>3FlbQDg/ZUB8w3FAV2fMQMiOnkCXfP4/T2I2QC6qskAnhOZA/NU2QCc2QEAI1UhAQcAbQRrziUBmHq5AXutSQWZDbkAZGWdAzt6lQYNptUDSFDNBoY4IQAYaQEDeT7Q/16HGP9GtXUCITrQ/Rxu0Pzhc6j9mpjZAX1X8P7tPQ0AqxS5BZTzZPye+m0B7Sa5AfPsPQRr/W0CYwBRBwDh3QMAmtD/nq6E/bJHGPxJ9UUDsy/dAoCYMQRM2a0BkAR9Boo5pQMV0VEArYu5A4kaMQAyTI0BQPRJAML3TQCKVCED85+tArObGP1BVP0EtJuVAdyKAQFjctkFQa2NBixMTQXyyjUFX8eo/IHelQTdFcEFo1zZAhagsQAO53EBIugRB0M+gQfhBgkH0MsJAbGlIQZXg+EHe6CZBsbA2QHMHOECtW6BAjE2oQUpZckBasZ1AtKl3QEZYIUHkip1AQX7TQPqF60GNuaE/USk2QGLF40Im65ZAmXqlQBGuSUC70KBAAneMQeK3aEB87MVA5NigQE/Wb0BO475A</binary>
                </binaryDataArray>
            </binaryDataArrayList>
        </spectrum>
        "#;
        let mut spectra = MzMLReader::with_file_id(0).parse(s.as_bytes()).await?;

        assert_eq!(spectra.len(), 1);
        let s = spectra.pop().unwrap();
        assert!(match s.precursors[0].inverse_ion_mobility {
            Some(x) => (x - 1.0786) < 0.0001,
            None => false,
        });

        // The rest of these assertions just make sure the integrity of the spectrum is maintained
        assert_eq!(s.id, "spectrum=8678309");
        assert_eq!(s.ms_level, 2);
        assert_eq!(s.representation, Representation::Centroid);
        assert_eq!(s.precursors.len(), 1);
        assert_eq!(s.precursors[0].charge, Some(2));
        assert!((s.precursors[0].mz - 457.723968) < 0.0001);
        assert_eq!(
            s.precursors[0].isolation_window,
            Some(Tolerance::Da(-1.5, 0.75))
        );
        assert!((s.scan_start_time - 25.066).abs() < 0.0001);
        assert_eq!(s.ion_injection_time, 0.0);
        assert_eq!(s.intensity.len(), s.mz.len());
        Ok(())
    }

    #[tokio::test]
    async fn parse_spectrum_issue_210() -> Result<(), MzMLError> {
        // Handle cases where both isolation window target m/z is set and different than selected ion m/z
        let s = r#"
        <spectrum id="spectrum=8678309" index="8678309" defaultArrayLength="102" dataProcessingRef="dp_sp_1">
            <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" />
            <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2" />
            <precursorList count="1">
                <precursor>
                    <isolationWindow>
                        <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="457.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="1.5" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="0.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                    </isolationWindow>
                    <selectedIonList count="1">
                        <selectedIon>
                            <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="457.723968505859" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                            <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="2" />
                        </selectedIon>
                    </selectedIonList>
                </precursor>
            </precursorList>
        </spectrum>
        "#;
        let mut spectra = MzMLReader::with_file_id(0).parse(s.as_bytes()).await?;

        assert_eq!(spectra.len(), 1);
        let s = spectra.pop().unwrap();
        assert!((s.precursors[0].mz - 457.723968) < 0.0001);
        assert_eq!(
            s.precursors[0].isolation_window,
            Some(Tolerance::Da(-1.5, 0.75))
        );

        // Check different ordering of fields in mzML
        let s = r#"
        <spectrum id="spectrum=8678309" index="8678309" defaultArrayLength="102" dataProcessingRef="dp_sp_1">
            <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" />
            <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2" />
            <precursorList count="1">
                <precursor>
                    <selectedIonList count="1">
                        <selectedIon>
                            <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="457.723968505859" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                            <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="2" />
                        </selectedIon>
                    </selectedIonList>
                    <isolationWindow>
                        <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="457.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="1.5" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="0.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                    </isolationWindow>
                </precursor>
            </precursorList>
        </spectrum>
        "#;
        let mut spectra = MzMLReader::with_file_id(0).parse(s.as_bytes()).await?;

        assert_eq!(spectra.len(), 1);
        let s = spectra.pop().unwrap();
        assert!((s.precursors[0].mz - 457.723968) < 0.0001);
        assert_eq!(
            s.precursors[0].isolation_window,
            Some(Tolerance::Da(-1.5, 0.75))
        );

        // Check fallback keeping iso window m/z of fields in mzML
        let s = r#"
        <spectrum id="spectrum=8678309" index="8678309" defaultArrayLength="102" dataProcessingRef="dp_sp_1">
            <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" />
            <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2" />
            <precursorList count="1">
                <precursor>
                    <isolationWindow>
                        <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="457.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="1.5" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="0.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                    </isolationWindow>
                    <selectedIonList count="1">
                        <selectedIon>
                            <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="2" />
                        </selectedIon>
                    </selectedIonList>
                </precursor>
            </precursorList>
        </spectrum>
        "#;
        let mut spectra = MzMLReader::with_file_id(0).parse(s.as_bytes()).await?;

        assert_eq!(spectra.len(), 1);
        let s = spectra.pop().unwrap();
        assert!((s.precursors[0].mz - 457.75) < 0.0001);
        assert_eq!(
            s.precursors[0].isolation_window,
            Some(Tolerance::Da(-1.5, 0.75))
        );
        Ok(())
    }

    /// Little-endian bytes of `values`, as 64-bit or 32-bit floats
    fn floats(values: &[f64], f64s: bool) -> Vec<u8> {
        match f64s {
            true => values.iter().flat_map(|v| v.to_le_bytes()).collect(),
            false => values.iter().flat_map(|v| (*v as f32).to_le_bytes()).collect(),
        }
    }

    fn zlib(bytes: &[u8]) -> Vec<u8> {
        use std::io::Write;
        let mut encoder = flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(bytes).unwrap();
        encoder.finish().unwrap()
    }

    /// One binary data array: its kind, float size and compression terms, and
    /// its payload (base64 of the bytes, zlib-compressed first with `zlib`)
    fn array(kind: &str, f64s: bool, zlib_bytes: bool, bytes: &[u8]) -> String {
        let (dtype, compression) = (if f64s { "MS:1000523" } else { "MS:1000521" }, if zlib_bytes { "MS:1000574" } else { "MS:1000576" });
        let payload = if zlib_bytes { zlib(bytes) } else { bytes.to_vec() };
        array_with(kind, &format!(r#"<cvParam accession="{dtype}"/><cvParam accession="{compression}"/>"#), &base64::encode(payload))
    }

    fn array_with(kind: &str, params: &str, blob: &str) -> String {
        format!(r#"<binaryDataArray>{params}<cvParam accession="{kind}"/><binary>{blob}</binary></binaryDataArray>"#)
    }

    /// An MS2 spectrum holding `arrays`
    fn spectrum(id: &str, arrays: &[String]) -> String {
        format!(
            r#"<spectrum id="{id}"><cvParam accession="MS:1000511" value="2"/><cvParam accession="MS:1000127"/>
               <scanList><scan><cvParam accession="MS:1000016" value="90" unitAccession="UO:0000010"/></scan></scanList>
               <binaryDataArrayList>{}</binaryDataArrayList></spectrum>"#,
            arrays.concat()
        )
    }

    fn mzml(spectra: &[String]) -> String {
        format!("<mzML><run><spectrumList>{}</spectrumList></run></mzML>", spectra.concat())
    }

    const MZ: &str = "MS:1000514";
    const INTENSITY: &str = "MS:1000515";
    const NOISE: &str = "MS:1002744";

    /// Arrays decoded in parallel batches come out exact and in file order,
    /// whatever the batch size: zlib and plain, 64- and 32-bit
    #[tokio::test]
    async fn arrays_decode_the_same_in_any_batch() -> Result<(), MzMLError> {
        let expected: Vec<(Vec<f64>, Vec<f64>)> = (0..7)
            .map(|n| {
                let mz = (0..20 + n).map(|i| 150.0 + i as f64 * 13.000_1 + n as f64 * 0.25).collect();
                let intensity = (0..20 + n).map(|i| 1000.0 + (i * 37 + n) as f64).collect();
                (mz, intensity)
            })
            .collect();
        let input = mzml(
            &expected
                .iter()
                .enumerate()
                .map(|(n, (mz, intensity))| {
                    let (zlib_bytes, f64s) = (n % 2 == 0, n % 3 != 0);
                    spectrum(
                        &format!("scan={n}"),
                        &[array(MZ, f64s, zlib_bytes, &floats(mz, f64s)), array(INTENSITY, false, zlib_bytes, &floats(intensity, false))],
                    )
                })
                .collect::<Vec<_>>(),
        );
        for batch in [usize::MAX, 1, 200] {
            let mut reader = MzMLReader::with_file_id(0);
            if batch != usize::MAX {
                reader.set_decode_batch_bytes(batch);
            }
            let spectra = reader.parse(input.as_bytes()).await?;
            assert_eq!(spectra.len(), expected.len(), "batch {batch}");
            for (n, (s, (mz, intensity))) in spectra.iter().zip(&expected).enumerate() {
                assert_eq!(s.id, format!("scan={n}"), "batch {batch}: file order");
                assert_eq!(s.mz, mz.iter().map(|&v| v as f32).collect::<Vec<_>>(), "batch {batch}, scan {n}");
                assert_eq!(s.intensity, intensity.iter().map(|&v| v as f32).collect::<Vec<_>>());
            }
        }
        Ok(())
    }

    /// With S/N on, intensities at that MS level are divided by the spectrum's
    /// own noise array
    #[tokio::test]
    async fn signal_to_noise_divides_by_the_noise_array() -> Result<(), MzMLError> {
        let input = mzml(&[spectrum(
            "scan=1",
            &[
                array(MZ, true, true, &floats(&[400.0, 500.0], true)),
                array(INTENSITY, false, true, &floats(&[100.0, 300.0], false)),
                array(NOISE, false, true, &floats(&[10.0, 20.0], false)),
            ],
        )]);
        let mut reader = MzMLReader::with_file_id(0);
        reader.set_signal_to_noise(Some(2));
        let spectra = reader.parse(input.as_bytes()).await?;
        assert_eq!(spectra[0].intensity, [10.0, 15.0]);
        Ok(())
    }

    /// Numpress arrays are refused by name: their bytes are not floats, and
    /// read as floats they would give wrong m/z values without a word
    #[tokio::test]
    async fn numpress_arrays_are_refused_by_name() {
        for accession in ["MS:1002312", "MS:1002313", "MS:1002314", "MS:1002746", "MS:1002747", "MS:1002748"] {
            let blob = base64::encode(floats(&[400.0, 500.0], true));
            let input = mzml(&[spectrum("scan=1", &[array_with(MZ, &format!(r#"<cvParam accession="{accession}"/>"#), &blob)])]);
            match MzMLReader::with_file_id(0).parse(input.as_bytes()).await {
                Err(MzMLError::UnsupportedCV(message)) => assert!(message.contains(accession) && message.contains("Numpress"), "{message}"),
                other => panic!("{accession} was accepted: {other:?}"),
            }
        }
    }

    #[tokio::test]
    async fn truncated_or_not_zlib_arrays_are_invalid_data() {
        let complete = zlib(&floats(&[400.0, 500.0, 600.0], true));
        for payload in [complete[..complete.len() - 4].to_vec(), b"these bytes are not a zlib stream".to_vec()] {
            let params = r#"<cvParam accession="MS:1000523"/><cvParam accession="MS:1000574"/>"#;
            let input = mzml(&[spectrum("scan=1", &[array_with(MZ, params, &base64::encode(payload))])]);
            let result = MzMLReader::with_file_id(0).parse(input.as_bytes()).await;
            assert!(
                matches!(&result, Err(MzMLError::IOError(e)) if e.kind() == std::io::ErrorKind::InvalidData),
                "{result:?}"
            );
        }
    }

    /// Two bad arrays in one batch, failing in different ways: the error is
    /// the earlier spectrum's, not whichever thread finished first
    #[tokio::test]
    async fn the_first_bad_array_in_file_order_is_reported() {
        let input = mzml(&[
            spectrum("scan=1", &[array_with(MZ, r#"<cvParam accession="MS:1000576"/>"#, "@@@@not base64")]),
            spectrum(
                "scan=2",
                &[array_with(MZ, r#"<cvParam accession="MS:1000574"/>"#, &base64::encode(b"these bytes are not a zlib stream"))],
            ),
        ]);
        let result = MzMLReader::with_file_id(0).parse(input.as_bytes()).await;
        assert!(matches!(result, Err(MzMLError::Base64Error(_))), "{result:?}");
    }
}
