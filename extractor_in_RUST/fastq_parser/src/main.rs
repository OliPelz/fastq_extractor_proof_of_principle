extern crate regex;
use regex::Regex;

use std::io;
use std::io::BufReader;
use std::io::BufRead;
use std::io::BufWriter;
use std::io::Write;
use std::fs::File;
use std::env;
use std::process;

fn main() {

    let re = Regex::new("ACC.{20,21}G").unwrap();

    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
       println!("error, argument missing: a fastq file");
       process::exit(1);
    }
    let fastq_file = &args[1];
    let f = File::open(fastq_file).unwrap();

    let file = BufReader::new(&f);
    let mut writer = BufWriter::new(io::stdout());
    let mut cnt = 1;

    //let mut fq_header: &str = "";
    let mut fq_header = String::from("");
    let mut fq_seq = String::from("");
    let mut fq_start = 0;
    let mut fq_stop = 0;
    let mut strand = String::from("");        

    let mut found_hit = false;

    let mut count_total = 0;
    let mut count_extracted = 0;

    let mut out_file = File::create("/tmp/out").unwrap();
    let mut statistics_file = File::create("/tmp/stats").unwrap();

    for line in file.lines() {
        let l = line.unwrap();
        let n = cnt % 4;
        cnt = cnt + 1;
        if  n == 1 {
            fq_header = l.clone();
            found_hit = false;
            count_total += 1;
        }
        if n == 2 {
            match re.find(&l) {
                None => continue,
                Some(mat) => {
                    found_hit = true;
                    //seq = mat.as_str();
                    fq_seq = l.clone();
                    fq_start = mat.start() + 3;
                    fq_stop = mat.end() - 1;
                    count_extracted += 1;
                },
            };
        }
        if n == 3 && found_hit {
            strand = l.clone()
        }
        if n ==0 && found_hit {
            out_file.write((&fq_header).as_bytes()).unwrap();
            out_file.write(b"\n").unwrap();
            out_file.write((&fq_seq[fq_start..fq_stop]).as_bytes()).unwrap();
            out_file.write(b"\n").unwrap();
            out_file.write((&strand).as_bytes()).unwrap();
            out_file.write(b"\n").unwrap();

            out_file.write((&l[fq_start..fq_stop]).as_bytes()).unwrap();
            out_file.write(b"\n").unwrap();
        }
    }

    // print statistics
    writeln!(statistics_file, "Total\tExtracted\n").unwrap();
    //statistics_file.write("{}\t{}",count_total.to_string(), count_extracted.to_string());
    //statistics_file.write(b"\n").unwrap();

    println!("Fastq data extracted successfully");
    println!("Total Reads in this file:\t {}", count_total);
    println!("Extracted Reads in this file with the matching pattern:\t {}",count_extracted);
    println!("The provided pattern worked in:\t {}", (count_extracted/count_total*100));
}
