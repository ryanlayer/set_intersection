#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <omp.h>
#include <string.h>
#include <math.h>

#include "../lib/bed.h"
#include "../lib/set_intersect.h"
#include "../lib/mt.h"
#include "../lib/timer.h"

#define MIN(a,b) ((a)>(b)?(b):(a))

int main(int argc, char *argv[]) {
	if (argc < 5) {
		fprintf(stderr, "usage: %s <u> <a> <b> <N> <p>\n", argv[0]);
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
	int reps = atoi(argv[4]);
	int n = atoi(argv[5]);


	if((chr_list_from_bed_file(&U_list, chrom_names, chrom_num, U_file) == 1) ||
	   (chr_list_from_bed_file(&A_list, chrom_names, chrom_num, A_file) == 1) ||
	   (chr_list_from_bed_file(&B_list, chrom_names, chrom_num, B_file) == 1) ){

		fprintf(stderr, "Error parsing bed files.\n");
		return 1;
	}

	unsigned int max = add_offsets(U_list, chrom_num);
	unsigned int bits = (int)( ceil(log(max)/log(2) ) );
	unsigned int mask = (2 << (bits-1)) - 1;

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

	start();
	// sort A so it can be searched by start
	start();
	qsort(A, A_size, sizeof(struct interval_triple),
			compare_interval_triples_by_start);
	qsort(B, B_size, sizeof(struct interval_triple),
			compare_interval_triples_by_start);

	unsigned int *A_len = (unsigned int *) malloc(
			A_size * sizeof(unsigned int));
	unsigned int *B_len = (unsigned int *) malloc(
			B_size * sizeof(unsigned int));

	// Set sized
	for (i = 0; i < A_size; i++)
		A_len[i] = A[i].end - A[i].start;
	for (i = 0; i < B_size; i++) 
		B_len[i] = B[i].end - B[i].start;

	memcpy(A_end, A, A_size * sizeof(struct interval_triple));
	qsort(A_end, A_size, sizeof(struct interval_triple),
			compare_interval_triples_by_end);
	stop();
	//unsigned long sort_seq = report();

	start();
	int O = count_intersections_bsearch_seq( A, A_end, A_size, B, B_size );
	stop();
	//unsigned long count_seq = report();

	//init_genrand((unsigned) time(NULL));
	//init_genrand(2);

	/*
	unsigned long rand_total_time = 0,
				  sort_total_time = 0,
				  intersect_total_time = 0;
				  */

	/*
	unsigned long rand_total_time_p[n],
				  sort_total_time_p[n],
				  intersect_total_time_p[n];
	bzero(rand_total_time_p, n * sizeof(unsigned long));
	bzero(sort_total_time_p, n * sizeof(unsigned long));
	bzero(intersect_total_time_p, n * sizeof(unsigned long));
	*/

	omp_set_num_threads(n);

	unsigned int seed = (unsigned) time(NULL);
	int *x = (int *) calloc(n, sizeof(int) );
	int th_id, j;

	int flag[n];
	bzero(flag, n * sizeof(int));

	unsigned int t[n];
	bzero(t, n * sizeof(unsigned int));

	unsigned long *mt_r;
	int mti_r;

	struct interval_triple *A_r;
	struct interval_triple *A_end_r;
	struct interval_triple *B_r;

	#pragma omp parallel for private(th_id, mt_r, mti_r,i,A_r,B_r,A_end_r)
	for (j = 0; j < reps; j++) {
		th_id = omp_get_thread_num();

		if (flag[th_id] == 0) {

			A_r = (struct interval_triple *)
					malloc(A_size * sizeof(struct interval_triple));
			A_end_r = (struct interval_triple *)
					malloc(A_size * sizeof(struct interval_triple));
			B_r = (struct interval_triple *)
					malloc(B_size * sizeof(struct interval_triple));

			mt_r = (unsigned long *) malloc (N * sizeof(unsigned long));
			init_genrand_omp(seed + th_id, mt_r, &mti_r);
			flag[th_id] = 1;
		}

		//start();
		for (i = 0; i < A_size; i++)
			A_r[i].start = get_rand_omp(max, mask, mt_r, &mti_r);

		for (i = 0; i < B_size; i++)
			B_r[i].start = get_rand_omp(max, mask, mt_r, &mti_r);
		//stop();
		//rand_total_time_p[th_id] += report();

		//start();
		qsort(A_r, A_size, sizeof(struct interval_triple),
			compare_interval_triples_by_start);
		qsort(B_r, B_size, sizeof(struct interval_triple),
			compare_interval_triples_by_start);

		for (i = 0; i < A_size; i++) 
			A_r[i].end = A_r[i].start + A_len[i];
		for (i = 0; i < B_size; i++)
			B_r[i].end = B_r[i].start + B_len[i];

		memcpy(A_end_r, A_r, A_size * sizeof(struct interval_triple));
		qsort(A_end_r, A_size, sizeof(struct interval_triple),
				compare_interval_triples_by_end);
		//stop();

		//sort_total_time_p[th_id] += report();

		//start();
		int r = count_intersections_bsearch_seq(
				A_r, A_end_r, A_size, B_r, B_size );
		//stop();
		//intersect_total_time_p[th_id] += report();

		t[th_id] += r;
		if (r >= O)
			x[th_id] = x[th_id] + 1;		
	}
	stop();

	int X = 0;

	for (i = 0; i < n; i++)
		X += x[i];

	unsigned int all_t = 0;
	for (i = 0; i < n; i++)
		all_t += t[i];


	double p = ( (double)(X + 1) ) / ( (double)(reps + 1) );

	fprintf(stderr, "O:%d\tp:%G\tE:%G\n", O, p, (double)all_t/(double)reps);
	printf("%d,%d,%d\tt:%lu\n", A_size, B_size, A_size + B_size, report());

	/*
	unsigned long rand_total_time_t = 0, 
				sort_total_time_t = 0,
			   	intersect_total_time_t = 0;

	for (i = 0; i < n; i++) {
		rand_total_time_t += rand_total_time_p[i];
		sort_total_time_t += sort_total_time_p[i];
		intersect_total_time_t += intersect_total_time_p[i];
	}

	double rand_total_time_a = (double)rand_total_time_t / (double)n,
			sort_total_time_a = (double)sort_total_time_t / (double)n,
			intersect_total_time_a = (double)intersect_total_time_t / (double)n;


	double  rand_avg_time = rand_total_time_a / reps,
			sort_avg_time = sort_total_time_a / reps,
			intersect_avg_time = intersect_total_time_a / reps;

	double total_avg_time = rand_avg_time + sort_avg_time + intersect_avg_time;

	double  rand_prop_time = rand_avg_time/total_avg_time,
			sort_prop_time = sort_avg_time/total_avg_time,
			intersect_prop_time = intersect_avg_time/total_avg_time;

	fprintf(stderr, "p:%G\tO:%d\n", p, O);
	printf("%d,%d,%d\tt:%G\tr:%G,%G\ts:%G,%G\ti:%G,%G\n", 
			A_size,
			B_size,
			A_size + B_size,
			//total_avg_time,
			rand_total_time_a + sort_total_time_a + intersect_total_time_a,
			rand_avg_time, rand_prop_time,
			sort_avg_time, sort_prop_time,
			intersect_avg_time, intersect_prop_time);
			*/

	return 0;
}
