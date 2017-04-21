```bash
fastq_parser --help
```

```bash
USAGE:
    fastq_parser [OPTIONS] --fastq-file <fastq_input_file>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -f, --fastq-file <fastq_input_file>              fastq input file to process
    -p, --pattern <match_pattern>                    PERL style regexp to extract sub sequences
    -c, --reverse-complement <reverse_complement>    set to 'yes' if reverse complement, otherwise (default) set to no

```
for example

```bash
fastq_parser -f ./TRAIL-Replicate1.fastq 
```
