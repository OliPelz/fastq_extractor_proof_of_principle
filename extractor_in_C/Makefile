all: extractor 

extractor: extractor.c kseq.h
	$(CC) -O2 extractor.c -o extractor -lz
	
clean:
		rm -f *.o extractor
