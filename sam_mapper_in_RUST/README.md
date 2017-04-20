```bash
sam_mapper --help
```

```bash
Usage:
    ./sam_mapper_in_RUST/sam_mapper/target/release/sam_mapper [OPTIONS]

mapper for CRISPRanalyzer

optional arguments:
  -h,--help             show this help message and exit
  -f,--fasta-file FASTA_FILE
                        Fasta Library Input File
  -s,--sam-file SAM_FILE
                        Sam Input file
  -m,--match-pattern MATCH_PATTERN
                        Mapping match pattern e.g. M{20,21}$
  -g,--geneid-pattern GENEID_PATTERN
                        GeneId pattern to parse, e.g. '_'
  -l,--logfile LOGFILE  Logfile filename
```

for example

```bash
./sam_mapper -f pilotscreen.fasta -s play_it_again.sam
```
