extern crate regex;
extern crate argparse;


use regex::Regex;
use std::fs::File;
use argparse::{ArgumentParser, Store};
use std::collections::HashSet;
use std::io::{BufReader, BufRead, BufWriter, Write};

fn main() {
    let mut fasta_file_arg = String::new();
    let mut sam_file_arg   = String::new();
    let mut mapping_match_pattern = String::from("M{20,21}$");
    let mut geneid_pattern = String::from("_");
    let mut logfile_out    = String::from("./log.out");
    {
        // put the argparsing in its own scope
        let mut cli_parser = ArgumentParser::new();
        cli_parser.set_description("mapper for CRISPRanalyzer");

        cli_parser.refer(&mut fasta_file_arg)
            .add_option(&["-f", "--fasta-file"], Store,
            "Fasta Library Input File").required();

        cli_parser.refer(&mut sam_file_arg)
            .add_option(&["-s", "--sam-file"], Store,
            "Sam Input file").required();

        cli_parser.refer(&mut mapping_match_pattern)
            .add_option(&["-m", "--match-pattern"], Store,
            "Mapping match pattern e.g. M{20,21}$");

        cli_parser.refer(&mut geneid_pattern)
            .add_option(&["-g", "--geneid-pattern"], Store,
            "GeneId pattern to parse, e.g. '_'");

        cli_parser.refer(&mut logfile_out)
            .add_option(&["-l", "--logfile"], Store,
            "Logfile filename");

        cli_parser.parse_args_or_exit();
    }

// first parse the fasta file

    let mut geneids = HashSet::<String>::new();

    let fasta_re = Regex::new(&format!(r"^>(.+){}", geneid_pattern))
            .expect("programmer error in accession regex");

    let fasta_file = BufReader::new(File::open(fasta_file_arg).expect("Problem opening fastq file"));

    for line in fasta_file.lines() {
        let ln = line.expect("programmer error in reading fasta line by line");

        geneids.extend(
            fasta_re.captures_iter(&ln).map(|captures: regex::Captures| // iterate over all Matches, which may have multiple capture groups each
                captures.get(1) // of this match, take the first capture group
                .expect("fasta regex match should have had first capture group")
                .as_str().to_owned() // make Owned copy of capture-group contents
            )
        );
    }

    // now to the sam parser
    // TODO: implement as XVS stream parser
    // our buffer for the sam parser
    let mut next_line = String::new();
    let mut sam_file = BufReader::new(File::open(&mut sam_file_arg).expect("Problem opening fastq file"));

    // fast forward the sam header to the beginning of the
    // alignment section - skip the header starting with @
    let mut alignment = String::from("");

    while sam_file.read_line(&mut next_line).unwrap() > 0 {
        if !next_line.starts_with('@') {
            break;
        }
        next_line.clear();
    }


    let mut count_total = 0;

    if !next_line.is_empty() {
        let mismatch_in_patt = mapping_match_pattern.contains('x') || mapping_match_pattern.contains('X');

        let sam_mismatch_re = Regex::new(r"MD:Z:([0-9]+)([A-Z]+)([0-9])+" ).expect("programmer error in mismatch regex");
        let match_string_re = Regex::new(r"(([0-9]+)([MID]))+").expect("programmer error in match regex");

        // thats how to construct a do { } while loop in RUST
        loop {
            count_total += 1;
            // ----------the basic algorithm starts here ---
            alignment = next_line.clone();
            
            // now split
            let al_arr: Vec<&str> = alignment.trim_right().split("\t").collect();
            //println!("{}", al_arr[2]);
            let gene_id = al_arr[2].split("_").nth(0).unwrap();

            // the sam file format is so BAD that a certain position of any optional field cannot be
            // predicted for sure, so we need to parse the whole line for the mismatch string
            // at least we know that we have to search from the right end to the left because in the
            // beginning we have mandantory fields (first 11)
            let found_mismatch = sam_mismatch_re.is_match(&alignment);

            // do some prechecks to safe computation time...skip the obvious
            let skip = (!mismatch_in_patt && found_mismatch || mismatch_in_patt && !found_mismatch);
            if !skip {
                // build / expand cigar string, e.g. 20M -> MMMMMMMMMMMMMMMMMMMM, 10M,1I,5D ->
                // MMMMMMMMMMIDDDDD
                let mut match_string = String::new();
                for caps in match_string_re.captures_iter(&al_arr[5]) {
                    println!("{}",&caps[2]);
                    let until: i32 = caps[2].parse().expect("programmer error: cannot convert string to number for iterating");
                    for i in 1..until {
                        //print!("{}", &caps[3]);
                    }
                    println!();
                    //println!("BINGO {}", &caps[1]);
                }
                // now introduce mismatches if needed
                if(found_mismatch){
                    
                }

            }

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
