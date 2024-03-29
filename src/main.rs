use anyhow::{Context, Result};
use clap::Parser;
use env_logger::Builder;
use itertools::Itertools;
use log::{info, LevelFilter};
use noodles::fasta;
use noodles::fasta::record::{Definition, Sequence};
use skc::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{stdout, BufReader, BufWriter};
use std::path::{Path, PathBuf};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Target sequence (smallest of the two genomes recommended)
    ///
    /// Can be compressed with gzip, bzip2, xz, or zstd
    #[arg()]
    target: String,
    /// Query sequence
    ///
    /// Can be compressed with gzip, bzip2, xz, or zstd
    #[arg()]
    query: String,
    /// Size of k-mers (max. 32)
    #[arg(short, long, default_value_t = 21, value_parser = clap::value_parser!(u64).range(1..=32))]
    kmer: u64,
    /// Output filepath(s); stdout if not present.
    #[clap(short, long)]
    pub output: Option<PathBuf>,
    /// u: uncompressed; b: Bzip2; g: Gzip; l: Lzma; z: Zstd
    ///
    /// Output compression format is automatically guessed from the filename extension. This option
    /// is used to override that
    #[clap(short = 'O', long, value_name = "u|b|g|l|z", value_parser = parse_compression_format, default_value="u")]
    pub output_type: Option<niffler::compression::Format>,
    /// Compression level to use if compressing output
    #[clap(short = 'l', long, value_parser = parse_level, default_value="6", value_name = "INT")]
    pub compress_level: niffler::Level,
}
fn main() -> Result<()> {
    let args = Args::parse();

    let mut log_builder = Builder::new();
    log_builder
        .filter(None, LevelFilter::Info)
        .format_module_path(false)
        .format_target(false)
        .init();

    let k = args.kmer as usize;

    let (reader, _compression) = niffler::from_path(Path::new(&args.target))
        .with_context(|| format!("Failed to open {}", args.target))?;
    let mut fa_reader = fasta::Reader::new(BufReader::new(reader));

    let mut target_kmers: HashMap<u64, KmerInfo> = HashMap::new();

    for (n_rec, record) in fa_reader.records().enumerate() {
        let record = record.with_context(|| {
            format!("Failed to parse record {} (zero-based) from target", n_rec)
        })?;
        let chrom = record.name();
        let seq = record.sequence().as_ref();
        for i in 0..seq.len() {
            let Some(kmer) = &seq.get(i..i + k) else {continue};
            let h = encode(kmer)[0];
            target_kmers.entry(h).or_default().add_pos(chrom, i);
        }
    }

    info!("{} unique k-mers in target", target_kmers.len());

    let (reader, _compression) = niffler::from_path(Path::new(&args.query))
        .with_context(|| format!("Failed to open {}", args.query))?;
    let mut fa_reader = fasta::Reader::new(BufReader::new(reader));

    let output_handle = match &args.output {
        None => match args.output_type {
            None => Box::new(stdout()),
            Some(fmt) => niffler::basic::get_writer(Box::new(stdout()), fmt, args.compress_level)
                .context("Failed to create writer to stdout")?,
        },
        Some(p) => {
            let out_fd = File::create(p)
                .map(BufWriter::new)
                .context(format!("Failed to create output {:?}", p))?;

            let fmt = match args.output_type {
                None => match args.output {
                    Some(p) => niffler::Format::from_path(&p),
                    None => niffler::Format::No,
                },
                Some(f) => f,
            };
            niffler::get_writer(Box::new(out_fd), fmt, args.compress_level)
                .context("Failed to create writer for output file")?
        }
    };

    let mut fa_writer = fasta::Writer::new(output_handle);

    let mut query_kmers: HashMap<u64, KmerInfo> = HashMap::new();

    for (n_rec, record) in fa_reader.records().enumerate() {
        let record = record.context(format!(
            "Failed to parse record {} (zero-based) from query",
            n_rec
        ))?;
        let chrom = record.name();
        let seq = record.sequence().as_ref();
        for i in 0..seq.len() {
            let Some(kmer) = &seq.get(i..i + k) else { continue };
            let h = encode(kmer)[0];
            if target_kmers.contains_key(&h) {
                query_kmers.entry(h).or_default().add_pos(chrom, i);
            }
        }
    }
    info!(
        "{} shared k-mers between target and query",
        query_kmers.len()
    );

    for (h, query_kmerinfo) in query_kmers {
        let kmer = decode(&[h], k);
        // safe to unwrap as we know the hash is in target
        let target_kmerinfo = target_kmers.get(&h).unwrap();
        let mut description = format!(
            "tcount={} qcount={} ",
            target_kmerinfo.count(),
            query_kmerinfo.count()
        );
        let target_positions = target_kmerinfo.positions.iter().join(",");
        let query_positions = query_kmerinfo.positions.iter().join(",");
        let pos_descr = format!("tpos={} qpos={}", target_positions, query_positions);
        description.push_str(&pos_descr);
        let definition = Definition::new(h.to_string(), Some(description));
        let seq = Sequence::from(kmer);
        let record = fasta::Record::new(definition, seq);

        fa_writer
            .write_record(&record)
            .context("Failed to write record to output")?;
    }

    Ok(())
}
