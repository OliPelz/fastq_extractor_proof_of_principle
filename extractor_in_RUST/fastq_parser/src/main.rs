extern crate regex;

use regex::Regex;

use std::io::{BufReader, BufRead, BufWriter, Write};
use std::fs::File;
use std::env;
use std::process;


fn main() {
    let re = Regex::new("ACC.{20,21}G").expect("programmer error in accession regex");

    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        println!("error, argument missing: a fastq file");
        process::exit(1);
    }
    let fastq_file = &args[1];
    let file = BufReader::new(File::open(fastq_file).expect("Problem opening fastq file"));

    //let mut fq_header: &str = "";
    let mut fq_header = String::from("");
    let mut fq_seq = String::from("");
    let mut fq_start = 0;
    let mut fq_stop = 0;
    let mut strand = String::from("");

    let mut found_hit = false;

    let mut count_total = 0;
    let mut count_extracted = 0;

    let mut out_file =
        BufWriter::new(File::create("/tmp/out").expect("problem opening output file"));

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
                    fq_seq = l.clone();
                    fq_start = mat.start() + 3;
                    fq_stop = mat.end() - 1;
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
            out_file.write_all((&fq_seq[fq_start..fq_stop]).as_bytes()).unwrap();
            out_file.write_all(b"\n").unwrap();
            out_file.write_all((&strand).as_bytes()).unwrap();
            out_file.write_all(b"\n").unwrap();

            out_file.write_all((&l[fq_start..fq_stop]).as_bytes()).unwrap();
            out_file.write_all(b"\n").unwrap();
        }
    }

    println!("Fastq data extracted successfully");
    println!("Total Reads in this file:\t {}", count_total);
    println!("Extracted Reads in this file with the matching pattern:\t {}",
             count_extracted);
    println!("The provided pattern worked in:\t {:.2}%",(count_extracted as f32 / count_total as f32 * 100.0) );
}
