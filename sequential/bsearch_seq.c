#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

#include "../lib/bed.h"
#include "../lib/set_intersect.h"

#define MIN(a,b) ((a)>(b)?(b):(a))

int main(int argc, char *argv[]) {
	
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

	trim(U_list, A_list, chrom_num);
	trim(U_list, B_list, chrom_num);

	int i;

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
	struct triple *AB = (struct triple *)
			malloc((2*A_size + 2*B_size)*sizeof(struct triple));
	//A and B points to AB, A to the beging and B to the interior, after A
	struct triple *A = AB;
	struct triple *B = AB + 2*A_size;

	map_intervals(A, A_array, A_size, U_array, U_size, 0 );
	map_intervals(B, B_array, B_size, U_array, U_size, 1 );

	// sort A and B so they can be ranked
	qsort(A, 2*A_size, sizeof(struct triple), compare_triple_lists);
	qsort(B, 2*B_size, sizeof(struct triple), compare_triple_lists);

	int *A_len = (int *) malloc(A_size * sizeof(int));
	int *B_len = (int *) malloc(B_size * sizeof(int));

	// Set ranks
	for (i = 0; i < 2*A_size; i++)
		A[i].rank = i/2;
	for (i = 0; i < 2*B_size; i++) 
		B[i].rank = i/2;

	// Get lengths
	for (i = 0; i < A_size; i++)
		A_len[i] = A[i*2 + 1].key - A[i*2].key;
	for (i = 0; i < B_size; i++)
		B_len[i] = B[i*2 + 1].key - B[i*2].key;

	struct interval *A_r = (struct interval *) malloc (
			A_size * sizeof(struct interval));
	struct interval *B_r = (struct interval *) malloc (
			B_size * sizeof(struct interval));


	for (i = 0; i < A_size; i++) {
		A_r[i].start = A[i*2].key;
		A_r[i].end =  A[i*2 + 1].key;
	}
	//qsort(A_r, A_size, sizeof(struct interval), compare_interval_by_start);

	for (i = 0; i < B_size; i++) {
		B_r[i].start = B[i*2].key;
		B_r[i].end = B[i*2 + 1].key;
	}
	//qsort(B_r, B_size, sizeof(struct interval), compare_interval_by_start);

	/*
	int c = 0;
	for (i = 0; i < A_size; i++) {
		// Search for the left-most interval in B with the start in A
		int lo = -1, hi = B_size, mid;
		while ( hi - lo > 1) {
			mid = (hi + lo) / 2;

			if ( B_r[mid].start < A_r[i].start ) 
				lo = mid;
			else
				hi = mid;

		}

		int left = hi;
		// Small hack to make our property hold
		if ( B_r[hi].start == A_r[i].start)
			left++;

		lo = -1;
		hi = B_size;
		while ( hi - lo > 1) {
			mid = (hi + lo) / 2;

			if ( B_r[mid].start < A_r[i].end ) 
				lo = mid;
			else
				hi = mid;
		}

		int right = hi;
		if ( B_r[hi].start == A_r[i].end)
			right++;


		int k;

		c += (right - left) + 
				( (left > 0)  && (A_r[i].start < B_r[left - 1].end) );

	}
	*/

	int c = count_intersections_bsearch(A_r, A_size, B_r, B_size);


	printf("%d\n",c);

	return 0;
}
