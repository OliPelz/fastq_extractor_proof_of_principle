//#![feature(alloc_system)]
//:extern crate alloc_system;
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
                               .takes_value(true))
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
    let patt = matches.value_of("PATTERN").unwrap_or("ACC.{20,21}G");
    let is_reverse_str = matches.value_of("REVCOMP").unwrap_or("no");


    let log_file_str = format!("{}_log.txt", fastq_base_name);
    let mut log_out_file =
        BufWriter::new(File::create(matches.value_of("LOGFILE").unwrap_or(&log_file_str)).expect("cannot create out log file"));

    let mut is_reverse = false;

    if is_reverse_str == "yes" {
        is_reverse = true;
    }



    let re = Regex::new(patt).expect("programmer error in accession regex");
    let file = BufReader::new(File::open(fastq_file).expect("Problem opening fastq file"));

    let mut fq_header = String::from("");
    let mut fq_seq = String::with_capacity(200);
    let mut fq_start = 0;
    let mut fq_stop = 0;
    let mut strand = String::from("");

    let mut found_hit = false;

    let mut count_total = 0;
    let mut count_extracted = 0;

    let mut out_file =
        BufWriter::new(File::create(fastq_out_file).expect("problem opening output file"));

    for (line_number, line) in file.lines().enumerate() {
        let l = line.expect("programmer error: no line to unwrap");

        let n = line_number % 4; // offset inside fastq record

        if n == 0 {
            // first line of record: header
            fq_header = l;
            found_hit = false;
            count_total += 1;
        } else if n == 1 {
            // second line of record
            match re.find(&l) {
                None => continue,
                Some(mat) => {
                    found_hit = true;
                    //seq = mat.as_str();
                    fq_start = mat.start() + 3;
                    fq_stop = mat.end() - 1;
                    fq_seq.clone_from(&l);
                    count_extracted += 1;
                }
            };
        } else if n == 2 && found_hit {
            // third line of record: strand
            strand = l
        } else if n == 3 && found_hit {
            // fourth/last line of record: store everything
            out_file.write_all((&fq_header).as_bytes()).unwrap();
            out_file.write_all(b"\n").unwrap();

            //out_file.write_all((&fq_seq[fq_start..fq_stop]).as_bytes()).unwrap();

            if is_reverse  {
                out_file.write_all(&alphabets::dna::revcomp(&fq_seq[fq_start..fq_stop].to_owned().into_bytes())).expect("reverse complement fails");
            }
            else {
                out_file.write_all((&fq_seq[fq_start..fq_stop]).as_bytes()).unwrap();
            }
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
