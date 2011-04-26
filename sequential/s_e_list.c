#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

#include "../lib/bed.h"
#include "../lib/set_intersect.h"

#define MIN(a,b) ((a)>(b)?(b):(a))

// See if A contains B
int contains(struct interval_triple *A, struct interval_triple *B) {
	return ( A->start <= B->start ) && ( A->end > B->end);
}

int main(int argc, char *argv[]) {
	//struct timeval t0_start, t0_end, t1_start, t1_end, t2_start, t2_end;
	//gettimeofday(&t0_start,0);
	//gettimeofday(&t1_start,0);
	
	if (argc < 4) {
		fprintf(stderr, "usage: order <u> <a> <b> <N>\n");
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

	struct chr_list *U_list, *A_list, *B_list;

	char *U_file = argv[1], *A_file = argv[2], *B_file = argv[3];

	if((chr_list_from_bed_file(&U_list, chrom_names, chrom_num, U_file) == 1) ||
	   (chr_list_from_bed_file(&A_list, chrom_names, chrom_num, A_file) == 1) ||
	   (chr_list_from_bed_file(&B_list, chrom_names, chrom_num, B_file) == 1) ){

		fprintf(stderr, "Error parsing bed files.\n");
		return 1;
	}

	int max = add_offsets(U_list, chrom_num);

	if (max == 0) {
		fprintf(stderr, "Max is zero.\n");
		return 1;
	}

	trim(U_list, A_list, chrom_num);
	trim(U_list, B_list, chrom_num);


	int A_size, B_size, U_size;

	struct bed_line *U_array, *A_array, *B_array;

	// Move the universe and intervals from linked lists to arrays
	U_size = chr_array_from_list(U_list, &U_array, chrom_num);
	A_size = chr_array_from_list(A_list, &A_array, chrom_num);
	B_size = chr_array_from_list(B_list, &B_array, chrom_num);


	// make one large array to hold these
	/* 
	 * We need to put both A and B into a single array then sort it
	 *
	 * Each interval becomes a triple: 
	 *   key:  offset
	 *   sample:  A (0) or B (1)
	 *   type:  start (0) or  end (1)
	 *   rank: order within
	 *
	 */
	struct interval_triple *A = (struct interval_triple *)
			malloc((A_size)*sizeof(struct interval_triple));
	struct interval_triple *B = (struct interval_triple *)
			malloc((A_size)*sizeof(struct interval_triple));

	map_to_interval_triple(A, A_array, A_size, U_array, U_size, 0 );
	map_to_interval_triple(B, B_array, B_size, U_array, U_size, 1 );

	int i;
	// sort A and B so they can be ranked
	qsort(A, A_size, sizeof(struct interval_triple),
			compare_interval_triples_by_start);
	qsort(B, B_size, sizeof(struct interval_triple),
			compare_interval_triples_by_start);

    // Set ranks
	for (i = 0; i < A_size; i++)
		A[i].rank = i;
	/*
	for (i = 0; i < B_size; i++) 
		B[i].rank = i;

	qsort(B, B_size, sizeof(struct interval_triple),
			compare_interval_triples_by_end);

	for (i = 0; i < A_size; i++) 
		printf("%d\t<%d,%d,%d>\t<%d,%d,%d>\n", i,
				A[i].rank, A[i].start, A[i].end,
				B[i].rank, B[i].start, B[i].end);
	*/


	// Find the number of sublists
	int n_lists = 1;
	for (i = 1; i < A_size; i++ ) {
		if ( contains(&A[i-1], &A[i]) ) 
			++n_lists;
	}
	printf("n_lists:%d\n", n_lists);


	/*
	struct pair *A_start = (struct pair *) malloc(A_size*sizeof(struct pair));
	struct pair *A_end = (struct pair *) malloc(A_size*sizeof(struct pair));

	for (i = 0; i < 2*A_size; i++) {
		printf("=r:%u t:%u k:%u s:%u\n", 
				A[i].rank,
				A[i].type,
				A[i].key,
				A[i].sample);
		if (A[i].type == 0) {
			A_start[ i/2 ].rank = A[i].rank;
			A_start[ i/2 ].key = A[i].key;
			printf("-%d %u, %u\n", i, A_start[ i/2 ].rank, A[i].rank);
		} else {
			A_end[ i/2 ].rank = A[i].rank;
			A_end[ i/2 ].key = A[i].key;
			printf("+%d %u, %u\n", i, A_end[ i/2 ].rank, A[i].rank);
		}
	}

	for (i = 0; i < A_size; i++) {
		printf("%u\t%u\n", A_start[i].rank, A_end[i].rank);
	}

	*/
	return 0;
}
