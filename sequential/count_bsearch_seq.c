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

int interval_triple_bsearch_end( struct interval_triple *A_end, 
								 int A_size,
								 unsigned int key) {
	int lo = -1, hi = A_size, mid;
	while ( hi - lo > 1) {
		mid = (hi + lo) / 2;

		if ( A_end[mid].end < key)
			lo = mid;
		else
			hi = mid;
	}

	return hi;
}

int interval_triple_bsearch_start( struct interval_triple *A_start, 
								   int A_size,
								   unsigned int key) {
	int lo = -1, hi = A_size, mid;
	while ( hi - lo > 1) {
		mid = (hi + lo) / 2;

		if ( A_start[mid].start < key )
			lo = mid;
		else
			hi = mid;
	}

	return hi;
}


int main(int argc, char *argv[]) {
	if (argc < 4) {
		fprintf(stderr, "usage: %s <u> <a> <b>\n", argv[0]);
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


	int i;

	int A_size, B_size, U_size;

	struct bed_line *U_array, *A_array, *B_array;

	// Move the universe and intervals from linked lists to arrays
	U_size = chr_array_from_list(U_list, &U_array, chrom_num);
	A_size = chr_array_from_list(A_list, &A_array, chrom_num);
	B_size = chr_array_from_list(B_list, &B_array, chrom_num);


	struct interval_triple *A = (struct interval_triple *)
			malloc(A_size * sizeof(struct interval_triple));
	struct interval_triple *A_end = (struct interval_triple *)
			malloc(A_size * sizeof(struct interval_triple));

	struct interval_triple *B = (struct interval_triple *)
			malloc(B_size * sizeof(struct interval_triple));

	map_to_interval_triple(A, A_array, A_size, U_array, U_size, 0 );
	map_to_interval_triple(B, B_array, B_size, U_array, U_size, 1 );

	// sort A and B so they can be ranked
	start();
	qsort(A, A_size, sizeof(struct interval_triple), compare_interval_triples_by_start);
	qsort(B, B_size, sizeof(struct interval_triple), compare_interval_triples_by_start);
	stop();
	unsigned long sort_seq = report();

	for (i = 0; i < A_size; i++)
		A[i].rank = i;

	memcpy(A_end, A, A_size * sizeof(struct interval_triple));
	qsort(A_end, A_size, sizeof(struct interval_triple), compare_interval_triples_by_end);

	// Weight each interval by the number of intervals that contain it
	unsigned int *A_w = (unsigned int *) malloc(A_size * sizeof(unsigned int));
	bzero(A_w, A_size * sizeof(unsigned int));

	int j;
	for (i = 0; i < A_size - 1; i++) {
		j = i + 1;
		while (A[i].end >= A[j].end) {
			A_w[j] = A_w[j] + 1;
			++j;
		}
	}

	int O = 0;

	for (i = 0; i < B_size; i++) {
		int start_pos = interval_triple_bsearch_end(A_end, A_size, B[i].start);
		int end_pos = interval_triple_bsearch_end(A, A_size, B[i].end);

		if ( B[i].start <= A_end[start_pos].start ) 
			O += A_w[ A_end[start_pos].rank ] + 1;

		if ( B[i].end >= A_end[start_pos].start ) 

				/*
		if (A_end[start_pos].rank != A[end_pos].rank)
		printf("%d\t%u,%u\t%u,%u,%u\t%u,%u,%u\n",
				i,
				B[i].start,	
				B[i].end,	
				A_end[start_pos].rank,
				A_end[start_pos].start,
				A_end[start_pos].end,
				A[end_pos].rank,
				A[end_pos].start,
				A[end_pos].end);
				*/
	}


		/*

	start();
	int O = count_intersections_bsearch(
			A_start, A_len, A_size,
			B_start, B_len, B_size );
	stop();
	unsigned long count_seq = report();

	unsigned long total = sort_seq + count_seq;


	fprintf(stderr,"O:%d\n", O);
	printf("t:%ld\tsort:%ld,%G\tsearch:%ld,%G\n",
			total,
			sort_seq, (double)sort_seq / (double)total,
			count_seq, (double)count_seq / (double)total
	  );


	*/
	return 0;
}
