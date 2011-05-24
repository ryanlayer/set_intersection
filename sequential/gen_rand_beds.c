#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>
#include <math.h>

#include "../lib/bed.h"
#include "../lib/set_intersect.h"
#include "../lib/mt.h"
#include "../lib/timer.h"

#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)<(b)?(b):(a))

int main(int argc, char *argv[]) {
	if (argc < 4) {
		fprintf(stderr, "usage: %s <u> <a size> <b size> <p>\n"
						"u\tuniverse\n"	
						"a size\tsize of first set\n"	
						"b size\tsize of second set\n"	
						"p\tprobability element in b intersects a\n",
						argv[0]);
		return 1;
	}

	int chrom_num = 24;

	/***********************REPLACE WITH INPUT FILE************************/	
	char *chrom_names[] = {
		"chr1",  "chr2",  "chr3", "chr4",  "chr5",  "chr6",  "chr7", "chr8",
		"chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
		"chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX",  "chrY"
	};
	/**********************************************************************/	

	struct chr_list *U_list;

	char *U_file = argv[1];

	if (chr_list_from_bed_file(&U_list, chrom_names, chrom_num, U_file) == 1) {
		fprintf(stderr, "Error parsing bed files.\n");
		return 1;
	}

	int U_size;
	struct bed_line *U_array;

	// Move the universe and intervals from linked lists to arrays
	U_size = chr_array_from_list(U_list, &U_array, chrom_num);


	// we need a few random numers.  One to pick the next section of the
	// universe (this will over represent smaller chromosomes)
	unsigned int bits = (int)( ceil(log(chrom_num)/log(2) ) );
	unsigned int mask = (2 << (bits-1)) - 1;

	int i;
	for (i = 0; i < U_size; i++) {
		unsigned int rand_chr = get_rand(chrom_num, mask);
		unsigned int chr_max = U_array[rand_chr].end; 
		unsigned int chr_bits = (int)( ceil(log(chr_max)/log(2) ) );
		unsigned int chr_mask = (2 << (chr_bits-1)) - 1;
		unsigned int rand_start = get_rand(chr_max, chr_mask);

		printf("%s\t%u\t%u\n", chrom_names[rand_chr], rand_start, rand_start + 50);
	}



	/*
	for (i = 0; i < U_size; i++) {
		printf("%u\t%u\t%u\t%u\n",
				U_array[i].chr, 
				U_array[i].start, 
				U_array[i].end, 
				U_array[i].offset);
	}
	*/

	return 0;
}
