extern crate regex;
extern crate clap;
extern crate bio;
extern crate needletail;


use needletail::{fastx};
use regex::bytes::Regex;

use std::io::{BufWriter, Write};
use std::fs::File;
use clap::{Arg, App};
use bio::alphabets;


fn main() {
    let matches = App::new("fastq_parser")
        .version("0.0.1")
        .author("Oliver P. <oliverpelz@gmail.com>")
        .about("Fast Fastq extractor")
        .arg(Arg::with_name("PATTERN")
            .short("p")
            .long("pattern")
            .value_name("match_pattern")
            .help("PERL style regexp to extract sub sequences")
            .takes_value(true))
        .arg(Arg::with_name("FASTQ")
            .help("fastq/fastq.gz input file to process")
            .short("f")
            .long("fastq-file")
            .value_name("fastq_input_file")
            .required(true))
        .arg(Arg::with_name("REVCOMP")
            .short("c")
            .long("reverse-complement")
            .value_name("reverse_complement")
            .help("set to 'yes' if reverse complement, otherwise (default) set to no")
            .takes_value(true))
        .arg(Arg::with_name("LOGFILE")
            .short("l")
            .long("log-file")
            .value_name("log_file")
            .takes_value(true)
            .help("file name of the log file to output"))

        .get_matches();

    let fastq_file = matches.value_of("FASTQ").expect("cannot get fastq file");
    let fastq_base_name = fastq_file.replace(".fastq", "");
    let fastq_out_file = format!("{}_extracted.fastq", fastq_base_name);

    // define some default arguments for non-required values
    let patt = matches.value_of("PATTERN").unwrap_or(r"ACC(.{20,21})G");
    let is_reverse_str = matches.value_of("REVCOMP").unwrap_or("no");


    let log_file_str = format!("{}_log.txt", fastq_base_name);
    let mut log_out_file =
        BufWriter::new(File::create(matches.value_of("LOGFILE").unwrap_or(&log_file_str)).expect("cannot create out log file"));

    let mut is_reverse = false;

    if is_reverse_str == "yes" {
        is_reverse = true;
    }


    let re = Regex::new(patt).expect("programmer error in accession regex");

    let mut count_total = 0;
    let mut count_extracted = 0;

    let mut out_file =
        BufWriter::new(File::create(fastq_out_file).expect("problem opening output file"));

    fastx::fastx_file(&fastq_file[..], |seq| {
        count_total += 1;
        match re.captures(&seq.seq).as_mut() {
            None => {},
            Some(caps) => {
                match caps.get(1).as_mut() {
                    None => {},
                    Some(mat) => {
                        out_file.write_all(b"@").unwrap();
                        out_file.write_all((&seq.id).as_bytes()).unwrap();
                        out_file.write_all(b"\n").unwrap();
                        if is_reverse {
                            out_file.write_all(&alphabets::dna::revcomp(mat.as_bytes())).expect("reverse complement fails");
                        } else {
                            out_file.write_all(mat.as_bytes()).unwrap();
                        }
                        out_file.write_all(b"\n+\n").unwrap();
                        let qual_ar = seq.qual.expect("cannot extract subseq");
                        out_file.write_all(&qual_ar[mat.start()..mat.end()]).unwrap();
                        out_file.write_all(b"\n").unwrap();
                        count_extracted += 1;
                    }
                }
            }
        }
    });

    log_out_file.write_all(b"Total\tExtracted\n").unwrap();
    log_out_file.write_all(count_total.to_string().as_bytes()).unwrap();
    log_out_file.write_all(b"\t").unwrap();
    log_out_file.write_all(count_extracted.to_string().as_bytes()).unwrap();
    log_out_file.write_all(b"\n").unwrap();

    println!("Fastq data extracted successfully");
    println!("Total Reads in this file:\t{}", count_total);
    println!("Extracted Reads in this file with the matching pattern:\t{}", count_extracted);
    println!("The provided pattern matched:\t{:.2}%", (count_extracted as f32 / count_total as f32 * 100.0) );
}
