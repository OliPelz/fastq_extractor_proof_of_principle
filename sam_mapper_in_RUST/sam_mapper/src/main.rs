extern crate regex;

use regex::Regex;
use std::env;
use std::process;
use std::fs::File;

use std::io::{BufReader, BufRead, BufWriter, Write};

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        println!("error, argument missing: a sam file");
        process::exit(1);
    }
    let sam_file_arg = &args[1];
    let mut next_line = String::new();
    let mut sam_file = BufReader::new(File::open(sam_file_arg).expect("Problem opening fastq file"));
   
// fast forward the sam header to the beginning of the
// alignment section
   let mut alignment = String::from("");
  
   while sam_file.read_line(&mut next_line).unwrap() > 0 {
	//println!("{}" , next_line);
	if ! next_line.starts_with('@') {
	  break;
	}
	next_line.clear();
   }
   

   let mut count_total = 0;

   if ! next_line.is_empty() {
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
