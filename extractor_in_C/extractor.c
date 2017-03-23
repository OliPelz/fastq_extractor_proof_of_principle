 // include later, when we extract fastq.gz files as well
 // #include <zlib.h>  
    #include <string.h>
    #include <stdlib.h>
    #include <stdio.h>  
    #include <regex.h>
      
    int main(int argc, char *argv[])  
    {  
    //    gzFile fp; 
    	FILE * fp; 
	int FASTQ_MAX_LINE=1024;
// structures needed for the regexp search
	char *pattern;
        regex_t    preg;
	int is_complement = 0;
        int        rc;
        size_t     nmatch = 2;
        regmatch_t pmatch[2];
	


	if (argc == 1) {  
            fprintf(stderr, "Usage: %s <regexp> <fastq/fastq.gz file> <is reverse complement yes/no>\n", argv[0]);  
            return 1;  
        }
        pattern = argv[1];
	if(strcmp(pattern, "default") == 0) {
	   pattern = "ACC.\\{20,21\\}G";
	}
// compile the regexp
        if (0 != (rc = regcomp(&preg, pattern, 0))) {
           printf("regcomp() failed, returning nonzero (%d)\n", rc);
           exit(EXIT_FAILURE);
        }


        //fp = gzopen(argv[2], "r"); // STEP 2: open the file handler  
        fp = fopen(argv[2], "r"); // STEP 2: open the file handler  
	if(strcmp(argv[3], "yes")) {
	   is_complement = 1;
	}

	char buffer[FASTQ_MAX_LINE];
        int counttotal = 0;
	int countextracted = 0;
	int linecnt = 1;
	int n = 1;
	int valid_rec = 0;
	char header[256];
	char seq[256];
	char strand[2];
	char quality[256];
	int seq_start;
	int seq_stop;
	FILE *fp_out;
	FILE *fp_stat;
        fp_out  = fopen("/tmp/out.fastq", "w"); 
        fp_stat = fopen("/tmp/statistic.txt", "w");  

	while(fgets(buffer, FASTQ_MAX_LINE, fp) != NULL) {
           n = linecnt++ % 4;
	   if(n == 1) {
		valid_rec = 0;
		strcpy(header, buffer);
	   }
	   else if (n == 2) {
	// do regexp search
		if (0 == (rc = regexec(&preg, buffer, nmatch, pmatch, 0))) {
// cut out sequence taking into account regexp positions -CC -G
// TODO: this need to be set dynamically depending on the regexp  
                      seq_start = pmatch[0].rm_so + 3;
                      seq_stop = pmatch[0].rm_eo; 	        
                      memcpy(&seq, buffer + seq_start,  seq_stop - seq_start - 3);
		      valid_rec = 1;
               }
	   }
	   else if (n == 3 && valid_rec) {
		strcpy(strand, buffer);
	   }
	   else if (n == 0 && valid_rec) {
		// print out all to file
		fputs(header, fp_out);
		fputs(seq, fp_out);
		fputc('\n', fp_out);
		fputs(strand, fp_out);
		memcpy(&quality, buffer + seq_start, seq_stop - seq_start - 1 );
		fputs(quality, fp_out);
		fputs("\n", fp_out);
	   }
        }
	fclose(fp);
	fclose(fp_out);
	fclose(fp_stat);	
    }
