#include <stdio.h>
#include <string.h>
#include <cstdlib>


// usage: pipe samfile | ./bam2
int main(int argc, char **argv)
{
	//FILE *ah = fopen(argv[2], "a"); // append
	//char *len = argv[1];
	char *buf = NULL;
	size_t getlineLen = 0;
	
	while(getline(&buf, &getlineLen, stdin) != -1)
	{
		char *readName = strtok(buf, "\t");
		char *foo = strtok(NULL, "\t");
		char *chr = strtok(NULL, "\t");
		char *start = strtok(NULL, "\t");
		char *q = strtok(NULL, "\t");
		char *cigar = strtok(NULL, "\t");
		char *nameChr = strtok(readName, ":-");
		char *nameStart = strtok(NULL, ":-");
		
		if(strcmp(nameChr, chr) || strcmp(nameStart,start))
		{
			fprintf(stderr, "inconsistent: %s %s %s %s %s %s\n", nameChr, nameStart, chr, start, cigar, q);
		}
		else
		{
			fprintf(stdout, "%s\t%s\n", chr, start);
		}
	}
	//fclose(ah);
	free(buf);
	return(0);
}
