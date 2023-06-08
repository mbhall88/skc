use std::collections::HashSet;
use std::fs::File;
use anyhow::Result;
use clap::Parser;
use noodles::fasta;
use shared::n_to_bits_lut;
use std::io::{BufReader, Write};
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
    query: Option<String>,
    /// Size of k-mers to use
    #[arg(short, long, default_value_t = 21, value_parser = clap::value_parser!(u64).range(1..=32))]
    kmer: u64,
    /// Write target k-mer index to file
    ///
    /// Useful if you are going to use this target for multiple queries
    #[arg(short = 'x', long = "index")]
    write_index: Option<PathBuf>
}
fn main() -> Result<()> {
    let args = Args::parse();
    let k = args.kmer as usize;

    let (reader, _compression) = niffler::from_path(Path::new(&args.target))?;
    let mut fa_reader = fasta::Reader::new(BufReader::new(reader));

    let mut kmers = HashSet::new();

    for record in fa_reader.records() {
        let record = record?;
        let seq = record.sequence().as_ref();
        for i in 0..seq.len() {
            let Some(kmer) = &seq.get(i..i + k) else {continue};
            let h = n_to_bits_lut(kmer)[0] as usize;
            kmers.insert(h);
        }
    }

    println!("{} distinct k-mers in target", kmers.len());

    if let Some(index) = args.write_index {
        let mut fd = File::create(index)?;
        for i in &kmers {
            writeln!(fd, "{}", i)?;
        }
    }

    if let Some(query) = args.query {
        let (reader, _compression) = niffler::from_path(Path::new(&query))?;
        let mut fa_reader = fasta::Reader::new(BufReader::new(reader));
        let mut n_shared = HashSet::new();

        for record in fa_reader.records() {
            let record = record?;
            let seq = record.sequence().as_ref();
            for i in 0..seq.len() {
                let Some(kmer) = &seq.get(i..i + k) else { continue };
                let h = n_to_bits_lut(kmer)[0] as usize;
                if kmers.contains(&h) {
                    n_shared.insert(h);
                }
            }
        }

        println!("{} shared k-mers between target and query", n_shared.len());
    }

    Ok(())
}
