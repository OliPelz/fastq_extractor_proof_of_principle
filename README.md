# fastq_extractor_proof_of_principle
this repo tries to optimize some Perl scripts which are part of the CRISPRAnalyzer R shiny package and which are too slow for production use in web applications. The original code's benchmark is as follow (please note: PERL script not included here)

PERL
```
time perl CRISPR-extract.pl "ACC(.{20,21})G" data/TRAIL-Replicate1.fastq no
```

```bash
real	0m45.006s
user	0m44.191s
sys	0m0.682s
``` 
first benchmark without much optimization.

RUST
```bash
$ time ./extractor_in_RUST/fastq_parser/target/release/fastq_parser /home/olip/Desktop/Oli/data/TRAIL-Replicate1.fastq 
```

output
```bash
real	0m14.334s
user	0m4.267s
sys	0m10.036s
```

C
```bash
time ./extractor default ../../data/TRAIL-Replicate1.fastq  no
```

output
```bash
real	0m4.409s
user	0m3.953s
sys	0m0.445s
```

TODO:  make the regexp parsing multithreaded in RUST and C
