use clap::{value_parser, Arg, Command, ValueHint};
use rayon::ThreadPoolBuilder;
use sage_cli::input::Input;
use sage_cli::memory::{gib, MemoryGuard};
use sage_cli::runner::Runner;

fn main() -> anyhow::Result<()> {
    env_logger::Builder::default()
        .filter_level(log::LevelFilter::Error)
        .parse_env(env_logger::Env::default().filter_or("SAGE_LOG", "error,sage=info"))
        .init();

    let matches = Command::new("sage")
        .version(clap::crate_version!())
        .author("Michael Lazear <michaellazear92@gmail.com>")
        .about("\u{1F52E} Sage \u{1F9D9} - Proteomics searching so fast it feels like magic!")
        .arg(
            Arg::new("parameters")
                .required_unless_present("gpu-info")
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help("Path to configuration parameters (JSON file)")
                .value_hint(ValueHint::FilePath),
        )
        .arg(
            Arg::new("mzml_paths")
                .num_args(1..)
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help(
                    "Paths to mzML files to process. Overrides mzML files listed in the \
                     configuration file.",
                )
                .value_hint(ValueHint::FilePath),
        )
        .arg(
            Arg::new("fasta")
                .short('f')
                .long("fasta")
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help(
                    "Path to FASTA database. Overrides the FASTA file \
                     specified in the configuration file.",
                )
                .value_hint(ValueHint::FilePath),
        )
        .arg(
            Arg::new("output_directory")
                .short('o')
                .long("output_directory")
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help(
                    "Path where search and quant results will be written. \
                     Overrides the directory specified in the configuration file.",
                )
                .value_hint(ValueHint::DirPath),
        )
        .arg(
            Arg::new("batch-size")
                .long("batch-size")
                .value_parser(value_parser!(u16).range(1..))
                .help("Number of files to load and search in parallel (default = # of CPUs/2)")
                .value_hint(ValueHint::Other),
        )
        .arg(
            Arg::new("parquet")
                .long("parquet")
                .action(clap::ArgAction::SetTrue)
                .help("Write search output in parquet format instead of tsv"),
        )
        .arg(
            Arg::new("annotate-matches")
                .long("annotate-matches")
                .action(clap::ArgAction::SetTrue)
                .help("Write matched fragments output file."),
        )
        .arg(
            Arg::new("write-pin")
                .long("write-pin")
                .action(clap::ArgAction::SetTrue)
                .help("Write percolator-compatible `.pin` output files"),
        )
        .arg(
            Arg::new("write-report")
                .long("write-report")
                .action(clap::ArgAction::SetTrue)
                .help("Write `.html` report file"),
        )
        .arg(
            Arg::new("write-annotator")
                .long("write-annotator")
                .action(clap::ArgAction::SetTrue)
                .help(
                    "Write `results.annotator.csv` / `crosslinks.annotator.csv` to view annotated \
                     spectra in Annotator (github.com/snijderlab/annotator)",
                ),
        )
        .arg(
            Arg::new("disable-telemetry")
                .long("disable-telemetry-i-dont-want-to-improve-sage")
                .action(clap::ArgAction::SetFalse)
                .help("Disable sending telemetry data"),
        )
        .arg(
            Arg::new("threads")
                .long("threads")
                .value_parser(value_parser!(u16).range(1..))
                .help("Maximum number of CPU threads to use (default = all logical CPUs)")
                .value_hint(ValueHint::Other),
        )
        .arg(
            Arg::new("gpu-info")
                .long("gpu-info")
                .action(clap::ArgAction::SetTrue)
                .help("Report whether this build can use a GPU, list the adapters it sees, and exit"),
        )
        .arg(
            Arg::new("max-memory")
                .long("max-memory")
                .value_parser(value_parser!(f64))
                .help("Stop the search when sage.exe uses more memory than this, in GB (also max_memory_gb)")
                .value_hint(ValueHint::Other),
        )
        .arg(
            Arg::new("estimate")
                .long("estimate")
                .action(clap::ArgAction::SetTrue)
                .help("Print the size of the database the search would build (peptides, fragments, memory) and exit"),
        )
        .arg(
            Arg::new("stack-size")
                .long("stack-size")
                .value_parser(value_parser!(u32).range(1..))
                .help("Set Rayon worker thread stack size in MiB (default: 2 MiB)")
                .value_hint(ValueHint::Other),
        )
        .help_template(
            "{usage-heading} {usage}\n\n\
             {about-with-newline}\n\
             Written by {author-with-newline}Version {version}\n\n\
             {all-args}{after-help}",
        )
        .get_matches();

    if matches.get_flag("gpu-info") {
        gpu_info();
        return Ok(());
    }

    let stack_size_mib = matches.get_one::<u32>("stack-size").copied().unwrap_or(2);
    let stack_size_bytes = stack_size_mib as usize * 1024 * 1024;
    log::trace!(
        "setting Rayon worker thread stack size to {} MiB",
        stack_size_mib
    );
    let threads = matches
        .get_one::<u16>("threads")
        .map(|&n| n as usize)
        .unwrap_or_else(num_cpus::get);
    log::info!("using {} CPU thread(s)", threads);
    ThreadPoolBuilder::new()
        .stack_size(stack_size_bytes)
        .num_threads(threads)
        .build_global()
        .expect("configure Rayon pool");

    let parallel = matches
        .get_one::<u16>("batch-size")
        .map(|&n| n as usize)
        .unwrap_or_else(|| (threads / 2).max(1));

    let parquet = matches.get_one::<bool>("parquet").copied().unwrap_or(false);
    let send_telemetry = matches
        .get_one::<bool>("disable-telemetry")
        .copied()
        .unwrap_or(true);

    let estimate = matches.get_flag("estimate");
    let max_memory = matches.get_one::<f64>("max-memory").copied();

    let input = Input::from_arguments(matches)?;
    let mut parameters = input.build()?;

    if estimate {
        return print_estimate(&mut parameters, parallel);
    }

    // the limit given on the command line, else the configuration's
    let limit = max_memory.filter(|gb| *gb > 0.0).or(parameters.max_memory_gb);
    if let Some(gb) = limit {
        log::info!("memory limit: {:.2} GB", gb);
    }
    let guard = MemoryGuard::start(limit);

    let runner = Runner::new(parameters, parallel)?;
    let tel = runner.run(parallel, parquet)?;
    log::info!("peak memory use: {:.2} GB", gib(guard.peak_bytes()));
    drop(guard);

    if send_telemetry {
        tel.send();
    }

    Ok(())
}

/// `--estimate`: the database the search would build, counted without
/// building it, and the spectrum files it would read
fn print_estimate(parameters: &mut sage_cli::input::Search, parallel: usize) -> anyhow::Result<()> {
    let fasta = sage_cli::runner::read_protein_databases(&mut parameters.database)?;
    let e = parameters.database.estimate(&fasta);
    let files = parameters.mzml_paths.len();
    let on_disk = sage_cli::memory::spectra_bytes(&parameters.mzml_paths);
    println!("database estimate (counted without building it; an upper bound, as the search merges duplicate forms):");
    println!("  proteins: {}", e.proteins);
    println!("  peptide forms (targets and decoys): {} ({} decoys)", e.peptides, e.decoys);
    println!("  fragments indexed: {}", e.fragments);
    println!("  database memory: about {:.2} GB", gib(e.bytes));
    println!(
        "  spectrum files: {} ({:.2} GB on disk), read {} at a time (--batch-size)",
        files,
        gib(on_disk),
        parallel.min(files.max(1))
    );
    if parameters.crosslink.is_some() {
        println!("  crosslink search: the candidate pairs are formed per spectrum and are not counted here");
    }
    Ok(())
}

/// `--gpu-info`: one `key: value` line per fact, for people and for sage-gui
fn gpu_info() {
    #[cfg(feature = "gpu")]
    {
        println!("gpu-support: yes");
        for adapter in sage_core::crosslink_gpu::adapters() {
            println!("gpu-adapter: {adapter}");
        }
    }
    #[cfg(not(feature = "gpu"))]
    println!("gpu-support: no (rebuild with `cargo build --release --features gpu`)");
}
