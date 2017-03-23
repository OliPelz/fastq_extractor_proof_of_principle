# fastq_extractor_proof_of_principle


first benchmark without much optimization

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
