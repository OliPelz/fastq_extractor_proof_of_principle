// #![feature(alloc_system)]
// extern crate alloc_system;
extern crate regex;
extern crate clap;

use regex::Regex;
use std::fs::File;
use clap::{Arg, App};
use std::collections::BTreeMap;
use std::io::{BufReader, BufRead, BufWriter, Write};


// print $seqgene "Gene\tCount\tdesigns-present\n";
// print $seqgene $target."\t".$reads{$target}{"genematch"}."\t".$reads{$target}{"targetmatched"}."\n";
// open ($stats, ">", $logfile) or die $!; 
// print $stats "Total\tMatched\n";
// print $stats $counttotal . "\t" . $countmatched . "\n";
// close($stats);


fn main() {
    let matches = App::new("sam_mapper")
        .version("0.0.1")
        .author("Oliver P. <oliverpelz@gmail.com>")
        .about("SAM File mapper")
        .arg(Arg::with_name("MATCHPATTERN")
            .short("m")
            .long("match-pattern")
            .value_name("match_pattern")
            .help("PERL style regexp to match CIGAR strings")
            .takes_value(true))
        .arg(Arg::with_name("FASTA")
            .help("FASTA input file to process")
            .short("f")
            .long("fasta-file")
            .value_name("fasta_input_file")
            .required(true))
        .arg(Arg::with_name("SAM")
            .help("SAM input file to process")
            .short("s")
            .long("sam-file")
            .value_name("sam_input_file")
            .required(true))
        .arg(Arg::with_name("GENEIDFSEPERATOR")
            .short("g")
            .long("geneid-fs-pattern")
            .value_name("geneid_fs_pattern")
            .help("the gene id field seperator pattern, defaults to '_'")
            .takes_value(true))
        .arg(Arg::with_name("LOGFILE")
            .short("l")
            .long("log-file")
            .value_name("log_file")
            .takes_value(true)
            .help("file name of the log file to output"))
        .get_matches();

    let fasta_file_arg = matches.value_of("FASTA").expect("cannot get fasta input file");
    let sam_file_arg = matches.value_of("SAM").expect("cannot get sam input file");
    let out_base_name = sam_file_arg.replace(".sam", "");

    // define some default arguments for non-required values
    let mapping_match_pattern = matches.value_of("MATCHPATTERN").unwrap_or(r"M{20,21}$");
    let geneid_pattern = matches.value_of("GENEIDFSEPERATOR").unwrap_or("_").to_owned();



    //let fasta_re = Regex::new(&format!(r"^>(.+){}", geneid_pattern))
    let fasta_re = Regex::new(r"^>(.+)")
            .expect("programmer error in accession regex");

    let mismatch_in_pattern = mapping_match_pattern.contains('x') ||
                              mapping_match_pattern.contains('X');

    let mut gene_matches = BTreeMap::<String, u32>::new();
    let mut ref_lib_ids = BTreeMap::<String, u32>::new();
    let mut targets_matched = BTreeMap::<String, u32>::new();
    //let mut geneIds = TreeMap::<String, u32>::new();

    // first parse the reference genome from the fasta file
    process_fasta(&fasta_file_arg, &fasta_re, geneid_pattern, &mut gene_matches, &mut ref_lib_ids);

    // now parse the samfile
    //let (mapped_geneids, count_total) =
    let (count_total, count_matched) =
        process_sam(&sam_file_arg, mismatch_in_pattern, &mapping_match_pattern, &mut gene_matches, &mut ref_lib_ids);



    let mut design_out_file =
        BufWriter::new(File::create(format!("{}-designs.txt", out_base_name)).expect("problem opening output file"));


    design_out_file.write_all(b"sgRNA\tCount\n").unwrap();
//    let mut uniqLibIds
    for (k, v) in &ref_lib_ids {
        design_out_file.write_all(k.replace("\"", "").as_bytes()).unwrap();
        design_out_file.write_all(b"\t").unwrap();
        design_out_file.write_all(v.to_string().as_bytes()).unwrap();
        design_out_file.write_all(b"\n").unwrap();

        if *v > 0 {
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

    // write log file
    let log_file_str = format!("{}_log.txt", out_base_name);
    let mut log_file =
        BufWriter::new(File::create(matches.value_of("LOGFILE").unwrap_or(&log_file_str)).expect("cannot create out log file"));
    log_file.write_all(b"Total\tMatched\n").unwrap();
    log_file.write_all(count_total.to_string().as_bytes()).unwrap();
    log_file.write_all(b"\t").unwrap();
    log_file.write_all(count_matched.to_string().as_bytes()).unwrap();
    log_file.write_all(b"\n").unwrap();

}

fn process_fasta(fasta_file: &str, fasta_re: &Regex, geneid_pattern : String, gene_matches : &mut BTreeMap<String, u32>, ref_lib_ids: &mut BTreeMap<String, u32>) {


    let fasta_file = BufReader::new(File::open(fasta_file).expect("Problem opening fastq file"));

    for line in fasta_file.lines() {
        let ln = line.expect("programmer error in reading fasta line by line");

        ref_lib_ids.extend(
            fasta_re.captures_iter(&ln)
                // iterate over all Matches, which may have multiple capture groups each
                .map(|captures: regex::Captures| {
                    let key = captures.get(1) // of this match, take the first capture group
                        .expect("fasta regex match should have had first capture group")
                        .as_str().to_owned(); // make Owned copy of capture-group contents
                    // add to gene_matches as well
                    gene_matches.insert(key.split(&geneid_pattern).nth(0).unwrap().to_string(),  0);
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
               ref_lib_ids: &mut BTreeMap<String, u32>)
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

        if al_arr[1] != "0" {
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

                match gene_matches.get_mut(al_arr[2].split("_").nth(0).expect("gene matches not working")) {
                    Some(v) => *v += 1,
                    None => println!("illegal gene id encountered '{}'", &al_arr[2].split("_").nth(0).expect("illegal gene id cannot split"))
                }
               //ref_lib_ids.get(&x).ok_or("illegal gene id encountered").map(|v| v += 1);
                match ref_lib_ids.get_mut(&al_arr[2].to_owned().clone()) {
                    Some(v) => *v += 1,
                    None => println!("illegal reference lib id encountered '{}'", &al_arr[2].split("_").nth(0).expect("cannot split ref_lib_ids"))
                }
            }
        }

        // --------- end of basic algorithm ---
    }

   // (mapped_geneids, count_total)
    (count_total, count_matched)
}
