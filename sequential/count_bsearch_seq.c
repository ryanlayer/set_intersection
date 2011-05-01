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

	// sort A so it can be searched by start
	start();
	qsort(A, A_size, sizeof(struct interval_triple),
			compare_interval_triples_by_start);
	/*
	qsort(B, B_size, sizeof(struct interval_triple),
			compare_interval_triples_by_start);
	*/

	// copy A to A_end, sort A_end so it can be searched by end
	memcpy(A_end, A, A_size * sizeof(struct interval_triple));
	qsort(A_end, A_size, sizeof(struct interval_triple),
			compare_interval_triples_by_end);
	stop();
	unsigned long sort_seq = report();

	start();
	int O = count_intersections_bsearch_seq( A, A_end, A_size, B, B_size );
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
