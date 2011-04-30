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
	qsort(A, A_size, sizeof(struct interval_triple),
			compare_interval_triples_by_start);
	/*
	qsort(B, B_size, sizeof(struct interval_triple),
			compare_interval_triples_by_start);
	*/
	/*  Looks like we don't need ranks
	for (i = 0; i < A_size; i++)
		A[i].rank = i;
	*/

	memcpy(A_end, A, A_size * sizeof(struct interval_triple));
	qsort(A_end, A_size, sizeof(struct interval_triple),
			compare_interval_triples_by_end);
	stop();
	unsigned long sort_seq = report();


	/* It appears ranking is not needed
	// Weight each interval by the number of intervals that contain it
	unsigned int *A_w = (unsigned int *) malloc(A_size * sizeof(unsigned int));
	bzero(A_w, A_size * sizeof(unsigned int));

	int j;
	for (i = 0; i < A_size - 1; i++) {
		j = i + 1;
		while ( (A[i].end >= A[j].end) && (j < A_size) ) {
			A_w[j] = A_w[j] + 1;
			++j;
		}
	}

	for (i = 0; i < A_size; i++)
		printf("%d\t%d (%u, %u)\t%d (%u, %u)\n", 
				i, A[i].rank, A[i].start, A[i].end,
				A_end[i].rank, A_end[i].start, A_end[i].end);
	*/
	int O = 0;

	start();
	for (i = 0; i < B_size; i++) {
		// Find the position of the last interval to END before the query
		// interval STARTS
		int a = interval_triple_bsearch_end(A_end, A_size, B[i].start);
		// The key in this search is the start of the query, and the search
		// space is the set of interval ends.
		// The result, a,  is either the 
		// (1) insert position:
		//		the seiries of  intervals A_end[0], ..., A_end[a-1] end before
		//		the query starts, and the intervals 
		//		A_end[a], ..., A_end[A_size - 1] end after the query starts.
		//		In this case a is also the size of the set 
		//		{A_end[0], ... ,A_end[a-1]}
		// (2) position of a match:
		//		since we are talking about closed intervals, the interval 
		//		A_end[a] does not end before the query starts, so 
		//		the seiries of intervals the end before the query starts is 
		//		is A_end[0], ..., A_end[a-1], and a is still the size of that
		//		set

		// Find the position of the first interval to START after the query
		// ENDS
		int b = interval_triple_bsearch_start(A, A_size, B[i].end);
		// The key in this search is the end of the query, and the search space
		// is the set of interval starts
		// The result, b, is either the
		// (1) insert position:
		//		in this case, A[b] starts after the query ends, so the series of
		//		intervals A[b], ..., A[A_size - 1] start after the querty and,
		//		and the size of that set is A_size - b
		// (2) position of a match:
		//		int this case, A[b + 1] starts after the query ends, and the
		//		size of that set is A_size - b - 1.  It is possible  that
		//		A[b] == A[b + 1] == A[b + 2] and so on, in which case we need a
		//		little while loop to move b past all of these

		//int start_pos = MAX(0,a);
		//int end_pos = MIN(A_size - 1,b);

		int num_cant_before = a; 
		//int num_cant_after = A_size - b;

		while ( A[b].start == B[i].end )
			++b;

		int num_cant_after = A_size - b;

		int num_left = A_size - num_cant_before - num_cant_after;
		O += num_left;
		/*
		printf("can b:%d\tcant a:%d\tcan:%d\t",
				num_cant_before,
				num_cant_after,
				num_left);

		printf("a:%d\tb:%d\t", a, b);

		printf("q:(%u,%u)\tSE:%d (%u,%u)\tSS:%d (%u,%u)\n", 
				B[i].start,
				B[i].end,
				A_end[start_pos].rank,
				A_end[start_pos].start,A_end[start_pos].end,
				A[end_pos].rank, 
				A[end_pos].start, A[end_pos].end
				); 
		*/
	}
	stop();
	unsigned long count_seq = report();

	unsigned long total = sort_seq + count_seq;

	fprintf(stderr,"O:%d\n", O);

	printf("%d,%d,%d\tT:%ld\t"
			"sort:%ld,%G\t"
			"search:%ld,%G\n",
			A_size,
			B_size,
			A_size + B_size,
			total,
			sort_seq, (double)sort_seq / (double)total,
			count_seq, (double)count_seq / (double)total
		  );

	return 0;
}
