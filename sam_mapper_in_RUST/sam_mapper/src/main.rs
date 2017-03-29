extern crate regex;
extern crate argparse;


use regex::Regex;
use std::env;
use std::process;
use std::fs::File;
use argparse::{ArgumentParser, Store};
use std::collections::HashSet;
use std::io::{BufReader, BufRead, BufWriter, Write};

fn main() {


    let mut fasta_file_arg = String::new();
    let mut sam_file_arg = String::new();
    let mut mapping_match_pattern = String::from("M{20,21}$");
    let mut geneid_pattern = String::from("_");
    let mut logfile_out = String::from("./log.out");
    {
        // put the argparsing in its own scope
        let mut cli_parser = ArgumentParser::new();
        cli_parser.set_description("mapper for CRISPRanalyzer");

        cli_parser.refer(&mut fasta_file_arg)
            .add_option(&["-f", "--fasta-file"], Store, "Fasta Library Input File")
            .required();

        cli_parser.refer(&mut sam_file_arg)
            .add_option(&["-s", "--sam-file"], Store, "Sam Input file")
            .required();

        cli_parser.refer(&mut mapping_match_pattern)
            .add_option(&["-m", "--match-pattern"], Store,
            "Mapping match pattern e.g. M{20,21}$");

        cli_parser.refer(&mut geneid_pattern).add_option(&["-g", "--geneid-pattern"],
                                                         Store,
                                                         "GeneId pattern to parse, e.g. '_'");

        cli_parser.refer(&mut logfile_out).add_option(&["-l", "--logfile"],
                                                      Store,
                                                      "Logfile filename");

        cli_parser.parse_args_or_exit();
    }

    // first parse the fasta file

    let mut geneids = HashSet::new();

    let fasta_re = Regex::new(format!("{}{}", String::from("^>(.+)\n$"), geneid_pattern).as_str())
        .expect("programmer error in accession regex");
    let fasta_file =
        BufReader::new(File::open(fasta_file_arg).expect("Problem opening fastq file"));
    let mut ln;
    {
        for line in fasta_file.lines() {
            ln = line.expect("programmer error in reading fasta line by line");

            let caps;
            caps = fasta_re.captures(&ln).unwrap();
            if caps.len() >= 1 {
                geneids.insert(caps.get(1).map_or("", |m| m.as_str()));
            }
        }
    }
    println!("{}", geneids.len());


    // our buffer
    let mut next_line = String::new();

    let mut sam_file =
        BufReader::new(File::open(&mut sam_file_arg).expect("Problem opening fastq file"));

    // fast forward the sam header to the beginning of the
    // alignment section
    let mut alignment = String::from("");

    while sam_file.read_line(&mut next_line).unwrap() > 0 {
        //println!("{}" , next_line);
        if !next_line.starts_with('@') {
            break;
        }
        next_line.clear();
    }


    let mut count_total = 0;

    if !next_line.is_empty() {
        // thats how to construct a do { } while loop in RUST
        loop {
            count_total += 1;
            // ----------the basic algorithm starts here ---
            alignment = next_line.clone();
            //print!("{}", alignment);
            // now split
            let mut al_arr: Vec<&str> = alignment.trim_right().split("\t").collect();

            let mut gene_id = al_arr[2].split("_").nth(0).unwrap();




            // --------- end of basic algorithm ---

            next_line.clear();
            // get next alignment
            if sam_file.read_line(&mut next_line).unwrap() == 0 {
                next_line.clear();
                break;
            }
        }
    }

    println!("Total\tMatched");
    println!("{}\t{}", count_total, count_total);


}
