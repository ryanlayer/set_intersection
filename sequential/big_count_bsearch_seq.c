#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>

#include "../lib/bed.h"
#include "../lib/set_intersect.h"
#include "../lib/mt.h"
#include "../lib/timer.h"

#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)<(b)?(b):(a))

int main(int argc, char *argv[]) {

	/*
	 * Load the set of quiereis from A, then search B in some number of rounds
	 */
	if (argc < 4) {
		fprintf(stderr, "usage: %s <U> <A> <B> <chunks>\n"
				"U\tUnivese\n"
				"A\tQuery set\n"
				"B\tDB set\n"
				"chunks\tNumber of B to process each round\n",
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

	struct chr_list *U_list, *A_list;

	char *U_file_name = argv[1], *A_file_name = argv[2], *B_file_name = argv[3];

	int chunk_size = atoi(argv[4]);

	if((chr_list_from_bed_file(&U_list, chrom_names, 
					chrom_num, U_file_name) == 1) ||
			(chr_list_from_bed_file(&A_list, chrom_names, 
					chrom_num, A_file_name) == 1) ){
		fprintf(stderr, "Error parsing bed files.\n");
		return 1;
	}

	FILE *B_file = fopen(B_file_name, "r");

	if (B_file == NULL) {
		fprintf(stderr, "Could not open file:%s\n", B_file_name);
		return 1;
	}   

	unsigned int max = add_offsets(U_list, chrom_num);

	if (max == 0) {
		fprintf(stderr, "Max is zero.\n");
		return 1;
	}

	trim(U_list, A_list, chrom_num);

	int A_size, U_size;

	struct bed_line *U_array, *A_array;

	// Move the universe and intervals from linked lists to arrays
	U_size = chr_array_from_list(U_list, &U_array, chrom_num);
	A_size = chr_array_from_list(A_list, &A_array, chrom_num);

	unsigned int *A_start = (unsigned int *) malloc(
			A_size * sizeof(unsigned int));
	unsigned int *A_len = (unsigned int *) malloc(
			A_size * sizeof(unsigned int));

	map_to_start_len_array(A_start, A_len, A_array, A_size, 
						   U_array, U_size);

	unsigned int *R = (unsigned int *) malloc(
			A_size * sizeof(unsigned int));
	bzero(R, A_size * sizeof(unsigned int));

	unsigned int *B_start = (unsigned int *) malloc(
			chunk_size * sizeof(unsigned int));
	unsigned int *B_end = (unsigned int *) malloc(
			chunk_size * sizeof(unsigned int));

	unsigned int B_curr_size;

	// fill B arrays up to either the remaining amount of B, or full B
	while ( map_start_end_from_file(
					B_file, B_start, B_end, chunk_size, &B_curr_size,
					U_array, U_size) ) {
		
		printf("*");

		qsort(B_start, B_curr_size, sizeof(unsigned int), compare_uints);
		qsort(B_end, B_curr_size, sizeof(unsigned int), compare_uints);

		big_count_intersections_bsearch_seq(A_start, A_len, A_size,
				B_start, B_end, B_curr_size, R);
	}

	int i, O = 0;
	for (i = 0; i < A_size; i++) 
		O += R[i];

	printf("O:%u\n", O);

	return 0;
}
