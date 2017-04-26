extern crate regex;
extern crate clap;
extern crate bio;

use regex::Regex;

use std::io::{BufReader, BufRead, BufWriter, Write};
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
                               .takes_value(true)
                               .default_value(r"ACC(.{20,21})G")
                          )
                          .arg(Arg::with_name("FASTQ")
                               .help("fastq input file to process")
                              .short("f")
                              .long("fastq-file")
                              .value_name("fastq_input_file")
                              .required(true))
                          .arg(Arg::with_name("REVCOMP")
                               .short("c")
                               .long("reverse-complement")
                              .value_name("reverse_complement")
                              .takes_value(true)
                              .possible_values(&["yes", "no"])
                              .default_value("no")
                              .help("set to 'yes' if reverse complement, otherwise (default) set to no")
                          )
                          .arg(Arg::with_name("LOGFILE")
                              .short("l")
                              .long("log-file")
                              .value_name("log_file")
                              .takes_value(true)
                              .help("file name of the log file to output"))

        .get_matches();

    // Input and output files
    let fastq_file = matches.value_of("FASTQ").expect("cannot get fastq file");
    let fastq_base_name = fastq_file.replace(".fastq", "");
    let file = BufReader::new(File::open(fastq_file).expect("Problem opening fastq file"));
    let fastq_out_file = format!("{}_extracted.fastq", fastq_base_name);
    let log_file_str = format!("{}_log.txt", fastq_base_name);
    let mut log_out_file = BufWriter::new(File::create(
            matches.value_of("LOGFILE").unwrap_or(&log_file_str)
            ).expect("cannot create out log file"));
    let mut out_file = BufWriter::new(File::create(
            fastq_out_file
            ).expect("problem opening output file"));

    // convert non-required args to usable form
    let patt = matches.value_of("PATTERN").expect("PATTERN should have a default value");
    let re = Regex::new(patt).expect("programmer error in accession regex");
    let is_reverse: bool = matches.value_of("REVCOMP") == Some("yes");


    let mut fq_header = String::from("");
    let mut fq_seq = String::with_capacity(200);
    let mut fq_start = 0;
    let mut fq_stop = 0;
    let mut strand = String::from("");

    let mut found_hit = false;

    let mut count_total = 0;
    let mut count_extracted = 0;


    for (line_number, line) in file.lines().enumerate() {
        let l = line.expect("programmer error: no line to unwrap");
        let n = line_number % 4; // offset inside fastq record

        if n == 0 {
            // first line of record: header
            fq_header = l;
            found_hit = false;
            count_total += 1;
        } else if n == 1 {
            match re.captures(&l).as_mut() {
                None => {},
                Some(caps) => {
                    match caps.get(1).as_mut() {
                        None => {},
                        Some(mat) => {
                           fq_seq  = mat.as_str().to_owned();
                           fq_start = mat.start();
                           fq_stop = mat.end();

                           count_extracted += 1;
                           found_hit = true;
                        }
                   }
                }
            }
        } else if n == 2 && found_hit {
            // third line of record: strand
            strand = l
        } else if n == 3 && found_hit {
            // fourth/last line of record: store everything
            out_file.write_all((&fq_header).as_bytes()).unwrap();
            out_file.write_all(b"\n").unwrap();

            //out_file.write_all((&fq_seq[fq_start..fq_stop]).as_bytes()).unwrap();

            let normalised_seq = if is_reverse {
                &alphabets::dna::revcomp(fq_seq.as_bytes())
            } else {
                fq_seq.as_bytes()
            };
            out_file.write_all(normalised_seq);
            out_file.write_all(b"\n").unwrap();

            out_file.write_all((&strand).as_bytes()).unwrap();
            out_file.write_all(b"\n").unwrap();
            out_file.write_all((&l[fq_start..fq_stop]).as_bytes()).unwrap();
            out_file.write_all(b"\n").unwrap();
        }
    }

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
