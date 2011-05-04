#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <mpi.h>
#include <math.h>
#include <string.h>

#include "../lib/bed.h"
#include "../lib/set_intersect.h"
#include "../lib/mt.h"
#include "../lib/timer.h"

#define MIN(a,b) ((a)>(b)?(b):(a))

int main(int argc, char *argv[]) {
	MPI_Status status;
	int rank, size, seen, buf;

	MPI_Init (&argc, &argv);    /* starts MPI */
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);  /* get current process id */
	MPI_Comm_size (MPI_COMM_WORLD, &size);  /* get number of processes */

	if (argc < 5) {
		fprintf(stderr, "usage: num_sim_scan_seq <u> <a> <b> <N>\n");
		return 1;
	}

	unsigned long seed = 0;

	start();
	if (rank == 0) {
		seed = (unsigned long) time(NULL);
		int j;
		for (j = 1; j < size; j++)
			MPI_Send(&seed, 1, MPI_UNSIGNED, j, 0, MPI_COMM_WORLD);

	} else {
		MPI_Recv(&seed, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, &status);
	}

	fprintf(stderr, "r:%d\ts:%lu\n", rank, seed);

	init_genrand(seed + rank);

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

	// sort A so it can be searched by start
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

	int O = count_intersections_bsearch_seq( A, A_end, A_size, B, B_size );


	int j, x = 0;
	for (j = 0; j < (reps / size); j++) {
		for (i = 0; i < A_size; i++)
			A[i].start = get_rand(max, mask);

		for (i = 0; i < B_size; i++)
			B[i].start = get_rand(max, mask);

		qsort(A, A_size, sizeof(struct interval_triple),
			compare_interval_triples_by_start);
		qsort(B, B_size, sizeof(struct interval_triple),
			compare_interval_triples_by_start);

		for (i = 0; i < A_size; i++) 
			A[i].end = A[i].start + A_len[i];
		for (i = 0; i < B_size; i++)
			B[i].end = B[i].start + B_len[i];

		memcpy(A_end, A, A_size * sizeof(struct interval_triple));
		qsort(A_end, A_size, sizeof(struct interval_triple),
				compare_interval_triples_by_end);
		int r = count_intersections_bsearch_seq( A, A_end, A_size, B, B_size );

		if (r >= O)
			++x;
	}

	if (rank == 0) {
		seen = 0;
		printf("%d\tgot:%u\tfrom:%d\n", rank, x, rank);
		while (seen < (size - 1)) {
			MPI_Recv(&buf, 1, MPI_UNSIGNED, MPI_ANY_SOURCE, 0,
			MPI_COMM_WORLD, &status);
			++seen;
			x += buf;
			printf("%d\tgot:%u\tfrom:%d\t%d\n", rank, buf, status.MPI_SOURCE,
					x);
		}

		double p = ((double)x + 1) / (double)(reps/size * size );
		printf("O:%d\tp:%f\n", O, p);
	} else {
		MPI_Send(&x, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD);
	}

	stop();
	printf("%d,%d,%d\tt:%lu\tr:%d\n", A_size, B_size, A_size + B_size, report(),rank);
	MPI_Finalize();

	return 0;
}
