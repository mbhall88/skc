use anyhow::Result;
use clap::Parser;
use noodles::fasta;
use noodles::fasta::record::{Definition, Sequence};
use shared::*;
use std::collections::HashSet;
use std::fs::File;
use std::io::{stdout, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

/// Simple program to greet a person
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Target sequence
    #[arg()]
    target: String,
    /// Query sequence
    #[arg()]
    query: String,
    /// Size of k-mers to use
    #[arg(short, long, default_value_t = 21, value_parser = clap::value_parser!(u64).range(1..=32))]
    kmer: u64,
    /// Output filepath(s); stdout if not present.
    #[clap(short, long)]
    pub output: Option<PathBuf>,
    /// u: uncompressed; b: Bzip2; g: Gzip; l: Lzma; z: Zstd
    ///
    /// Will attempt to infer the output compression format automatically from the filename
    /// extension. This option is used to override that. If writing to stdout, the default is
    /// uncompressed
    #[clap(short = 'O', long, value_name = "u|b|g|l|z", value_parser = parse_compression_format, default_value="u")]
    pub output_type: Option<niffler::compression::Format>,
    /// Compression level to use if compressing output
    #[clap(short = 'l', long, value_parser = parse_level, default_value="6", value_name = "INT")]
    pub compress_level: niffler::Level,
}
fn main() -> Result<()> {
    let args = Args::parse();
    let k = args.kmer as usize;

    let (reader, _compression) = niffler::from_path(Path::new(&args.target))?;
    let mut fa_reader = fasta::Reader::new(BufReader::new(reader));

    // todo: add count = switch from HashSet to BTreeMap
    let mut target_kmers = HashSet::new();

    for record in fa_reader.records() {
        let record = record?;
        let seq = record.sequence().as_ref();
        for i in 0..seq.len() {
            let Some(kmer) = &seq.get(i..i + k) else {continue};
            let h = encode(kmer)[0];
            target_kmers.insert(h);
        }
    }

    println!("{} unique k-mers in target", target_kmers.len());

    let (reader, _compression) = niffler::from_path(Path::new(&args.query))?;
    let mut fa_reader = fasta::Reader::new(BufReader::new(reader));
    let mut query_kmers = HashSet::new();

    let mut output_handle = match &args.output {
        None => match args.output_type {
            None => Box::new(stdout()),
            Some(fmt) => niffler::basic::get_writer(Box::new(stdout()), fmt, args.compress_level)?,
        },
        Some(p) => {
            let out_fd = File::create(p).map(BufWriter::new)?;

            let fmt = match args.output_type {
                None => match args.output {
                    Some(p) => niffler::Format::from_path(&p),
                    None => niffler::Format::No,
                },
                Some(f) => f,
            };
            niffler::get_writer(Box::new(out_fd), fmt, args.compress_level)?
        }
    };

    let mut fa_writer = fasta::Writer::new(output_handle);

    for record in fa_reader.records() {
        let record = record?;
        let seq = record.sequence().as_ref();
        for i in 0..seq.len() {
            let Some(kmer) = &seq.get(i..i + k) else { continue };
            let h = encode(kmer)[0];
            query_kmers.insert(h);
        }

        let shared_kmers: Vec<_> = target_kmers.intersection(&query_kmers).collect();

        println!(
            "{} shared k-mers between target and query",
            shared_kmers.len()
        );

        for h in shared_kmers {
            let kmer = decode(&[*h], k);
            let definition = Definition::new(h.to_string(), Some(format!("k={}", k)));
            let seq = Sequence::from(kmer);
            let record = fasta::Record::new(definition, seq);

            fa_writer.write_record(&record)?;
        }
    }

    Ok(())
}
