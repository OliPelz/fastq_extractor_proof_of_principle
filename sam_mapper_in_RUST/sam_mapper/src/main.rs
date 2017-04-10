// #![feature(alloc_system)]
// extern crate alloc_system;
extern crate regex;
extern crate argparse;

use regex::Regex;
use std::fs::File;
use argparse::{ArgumentParser, Store};
use std::collections::HashSet;
use std::collections::BTreeMap;
use std::io::{BufReader, BufRead, BufWriter, Write};

// print $seqgene "Gene\tCount\tdesigns-present\n";
// print $seqgene $target."\t".$reads{$target}{"genematch"}."\t".$reads{$target}{"targetmatched"}."\n";
// open ($stats, ">", $logfile) or die $!; 
// print $stats "Total\tMatched\n";
// print $stats $counttotal . "\t" . $countmatched . "\n";
// close($stats);


fn main() {
    // buffers to hold parsed arguments
    let mut fasta_file_arg = String::new();
    let mut sam_file_arg = String::new();
    let mut mapping_match_pattern = String::from("M{20,21}$");
    let mut geneid_pattern = String::from("_");
    let mut logfile_out = String::from("./log.out");
// TODO: change argparse to clap as suggested by Jules
    parse_args(&mut fasta_file_arg,
               &mut sam_file_arg,
               &mut mapping_match_pattern,
               &mut geneid_pattern,
               &mut logfile_out);

    //let fasta_re = Regex::new(&format!(r"^>(.+){}", geneid_pattern))
    let fasta_re = Regex::new(r"^>(.+)")
            .expect("programmer error in accession regex");

    let mismatch_in_pattern = mapping_match_pattern.contains('x') ||
                              mapping_match_pattern.contains('X');

    let mut gene_matches = BTreeMap::<String, u32>::new();
    let mut ref_libIds = BTreeMap::<String, u32>::new();
    let mut targets_matched = BTreeMap::<String, u32>::new();
    //let mut geneIds = TreeMap::<String, u32>::new();

    // first parse the reference genome from the fasta file
    process_fasta(&fasta_file_arg, &fasta_re, geneid_pattern, &mut gene_matches, &mut ref_libIds);

    // now parse the samfile
    //let (mapped_geneids, count_total) =
    let (count_total, count_matched) =
        process_sam(&sam_file_arg, mismatch_in_pattern, &mapping_match_pattern, &mut gene_matches, &mut ref_libIds);

    let out_base_name = sam_file_arg.replace(".sam", "");

    let mut design_out_file =
        BufWriter::new(File::create(format!("{}-designs.txt", out_base_name)).expect("problem opening output file"));


    design_out_file.write_all(b"sgRNA\tCount\n").unwrap();
//    let mut uniqLibIds
    for (k, v) in &ref_libIds {
        design_out_file.write_all(k.replace("\"", "").as_bytes()).unwrap();
        design_out_file.write_all(b"\t").unwrap();
        design_out_file.write_all(v.to_string().as_bytes()).unwrap();
        design_out_file.write_all(b"\n").unwrap();

        if(*v > 0) {
            let gid = k.split("_").nth(0).unwrap().to_string();
            *targets_matched.entry(gid).or_insert(0) += 1;
        }
        //println!("{}\t{}", k.replace("\"", ""), v);
    }

    let mut genes_out_file =
        BufWriter::new(File::create(format!("{}-genes.txt", out_base_name)).expect("problem opening output file"));
    genes_out_file.write_all(b"Gene\tCount\tdesigns-present\n").unwrap();
    for (k, v) in &gene_matches {
        genes_out_file.write_all(k.as_bytes()).unwrap();
        genes_out_file.write_all(b"\t").unwrap();
        genes_out_file.write_all(v.to_string().as_bytes()).unwrap();
        genes_out_file.write_all(b"\t").unwrap();
        genes_out_file.write_all(targets_matched.get(k).unwrap().to_string().as_bytes()).unwrap();
        genes_out_file.write_all(b"\n").unwrap();

    }


//   foreach $target ( sort keys %reads) {
// print $seqgene $target."\t".$reads{$target}{"genematch"}."\t".$reads{$target}{"targetmatched"}."\n";

    // write log file
    let mut log_file =
        BufWriter::new(File::create(format!("{}_log.txt", fasta_file_arg)).expect("problem opening output file"));
    log_file.write_all(b"Total\tMatched\n").unwrap();
    log_file.write_all(b"\n").unwrap();
    log_file.write_all(count_total.to_string().as_bytes()).unwrap();
    log_file.write_all(b"\t").unwrap();
    log_file.write_all(count_matched.to_string().as_bytes()).unwrap();
    log_file.write_all(b"\n").unwrap();


    // FIXME: two times "count_total"?
    //println!("Total\tMatched");
    //println!("{}\t{}", count_total, count_total);
}


fn parse_args(fasta_file_arg: &mut String,
              sam_file_arg: &mut String,
              mapping_match_pattern: &mut String,
              geneid_pattern: &mut String,
              logfile_out: &mut String) {
    // put the argparsing in its own scope
    let mut cli_parser = ArgumentParser::new();
    cli_parser.set_description("mapper for CRISPRanalyzer");

    cli_parser.refer(fasta_file_arg)
        .add_option(&["-f", "--fasta-file"], Store, "Fasta Library Input File")
        .required();

    cli_parser.refer(sam_file_arg)
        .add_option(&["-s", "--sam-file"], Store, "Sam Input file")
        .required();

    cli_parser.refer(mapping_match_pattern)
        .add_option(&["-m", "--match-pattern"], Store, "Mapping match pattern e.g. M{20,21}$");

    cli_parser.refer(geneid_pattern)
        .add_option(&["-g", "--geneid-pattern"], Store, "GeneId pattern to parse, e.g. '_'");

    cli_parser.refer(logfile_out).add_option(&["-l", "--logfile"], Store, "Logfile filename");

    cli_parser.parse_args_or_exit();
}


fn process_fasta(fasta_file: &str, fasta_re: &Regex, geneid_pattern : String, gene_matches : &mut BTreeMap<String, u32>, ref_libIds: &mut BTreeMap<String, u32>) {


    let fasta_file = BufReader::new(File::open(fasta_file).expect("Problem opening fastq file"));

    for line in fasta_file.lines() {
        let ln = line.expect("programmer error in reading fasta line by line");

        ref_libIds.extend(
            fasta_re.captures_iter(&ln)
                // iterate over all Matches, which may have multiple capture groups each
                .map(|captures: regex::Captures| {
                    let key = captures.get(1) // of this match, take the first capture group
                        .expect("fasta regex match should have had first capture group")
                        .as_str().to_owned(); // make Owned copy of capture-group contents
                    // add to gene_matches as well
                    gene_matches.insert(key.split("_").nth(0).unwrap().to_string(),  0);
                    (key, 0)
                }
                )
        );
    }
}


fn process_sam(sam_file: &str,
               mismatch_in_pattern: bool,
               mapping_match_pattern: &str,
               gene_matches : &mut BTreeMap<String, u32>,
               ref_libIds: &mut BTreeMap<String, u32>)
               -> (u32,u32) {

//-> (HashMap<String, i32>, u32) {

    // our buffer for the sam parser
    let sam_file = BufReader::new(File::open(sam_file).expect("Problem opening fastq file"))
        .lines();

    let sam_mismatch_re =
        Regex::new(r"MD:Z:([0-9]+)([A-Z]+)[0-9]+").expect("programmer error in mismatch regex");
    let match_string_re = Regex::new(r"([0-9]+)([MIDSH])").expect("programmer error in match regex");
    let mapping_match_re =
        Regex::new(mapping_match_pattern).expect("programmer error in mapping match regexp");

    let mut count_total : u32 = 0;
    let mut count_matched : u32 = 0;
    for l in sam_file {
        let next_line = l.expect("io-error reading from samfile");

        // fast forward the sam header to the beginning of the
        // alignment section - skip the header starting with @
        if next_line.starts_with('@') {
            continue;
        }

        // ----------the basic algorithm starts here ---
        // now split
        let al_arr: Vec<&str> = next_line.trim_right().split("\t").collect();
        //only count the mapped read if 2nd field, the FLAG, indicates an alignment that is neither rev-complementary, nor unmapped, nor a mutliple alignment (FLAG = 0)

        if(al_arr[1] != "0") {
            continue;
        }

        count_total += 1;
        //println!("{}", al_arr[2]);
        //let gene_id = al_arr[2].split("_").nth(0).unwrap();

        let mut found_mismatch = false;
        // the sam file format is so BAD that a certain position of any optional field cannot be
        // predicted for sure, so we need to parse the whole line for the mismatch string
        // at least we know that we have to search from the right end to the left because in the
        // beginning we have mandantory fields (first 11)

        let mut mm_positions: Vec<usize> = Vec::new();
        for caps in sam_mismatch_re.captures_iter(&next_line) {
            let mm_pos: i32 = caps[1].parse().expect("programmer error: cannot parse string to number for iterating");
            mm_positions.push(mm_pos as usize);

            found_mismatch = true;
        }

        // do some prechecks to save computation time...skip the obvious
        let skip = !mismatch_in_pattern && found_mismatch || mismatch_in_pattern && !found_mismatch;
        if !skip {
            // build / expand cigar string, e.g. 20M -> MMMMMMMMMMMMMMMMMMMM, 10M,1I,5D ->
            // MMMMMMMMMMIDDDDD, 20M1D =
            let mut match_string = String::new();
            for caps in match_string_re.captures_iter(&al_arr[5]) {
                //println!("{}", &caps[1]);
                let until_pos: i32 = caps[1].parse().expect("programmer error: cannot convert string to number for iterating");
                for _ in 0..until_pos {
                    match_string.push_str(&caps[2]);
                }
            }

            // now introduce mismatches int the string if needed
            if found_mismatch {
                for pos in mm_positions {
                    // TODO: next line is not compiling
                    match_string.insert_str(pos, "X");
                }
            }

            // now apply input mapping regex


            if mapping_match_re.is_match(&match_string) {
                count_matched += 1;

                match gene_matches.get_mut(al_arr[2].split("_").nth(0).unwrap()) {
                    Some(v) => *v += 1,
                    None => println!("illegal gene id encountered '{}'", &al_arr[2].split("_").nth(0).unwrap())
                }
               //ref_libIds.get(&x).ok_or("illegal gene id encountered").map(|v| v += 1);
                match ref_libIds.get_mut(&al_arr[2].to_owned().clone()) {
                    Some(v) => *v += 1,
                    None => println!("illegal reference lib id encountered '{}'", &al_arr[2].split("_").nth(0).unwrap())
                }
            }
        }

        // --------- end of basic algorithm ---
    }

   // (mapped_geneids, count_total)
    (count_total, count_matched)
}
