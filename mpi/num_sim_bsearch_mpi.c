#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <mpi.h>

#include "../lib/bed.h"
#include "../lib/set_intersect.h"
#include "../lib/mt.h"
#include "../lib/timer.h"

#define MIN(a,b) ((a)>(b)?(b):(a))

int main(int argc, char *argv[]) {
	if (argc < 5) {
		fprintf(stderr, "usage: num_sim_scan_seq <u> <a> <b> <N>\n");
		return 1;
	}

	MPI_Status status;
	int rank, size;

	MPI_Init (&argc, &argv);    /* starts MPI */
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);  /* get current process id */
	MPI_Comm_size (MPI_COMM_WORLD, &size);  /* get number of processes */

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

	unsigned int *A_len = (unsigned int *) malloc(
			A_size * sizeof(unsigned int));
	unsigned int *B_len = (unsigned int *) malloc(
			B_size * sizeof(unsigned int));

	unsigned int *A_start = (unsigned int *) malloc(
			A_size * sizeof(unsigned int));
	unsigned int *B_start = (unsigned int *) malloc(
			B_size * sizeof(unsigned int));

	// Set sized
	for (i = 0; i < A_size; i++)
		A_start[i] = A[i*2].key;
	for (i = 0; i < B_size; i++) 
		B_start[i] = B[i*2].key;

	// Get lengthsrank = i/2;
	for (i = 0; i < A_size; i++)
		A_len[i] = A[i*2 + 1].key - A[i*2].key;
	for (i = 0; i < B_size; i++)
		B_len[i] = B[i*2 + 1].key - B[i*2].key;


	unsigned int O = count_intersections_bsearch(
			A_start, A_len, A_size,
			B_start, B_len, B_size );

	init_genrand((unsigned) time(NULL) + rank);

	int x = 0, j;

	unsigned int *A_rand = (unsigned int *) malloc(
			A_size * sizeof(unsigned int));
	unsigned int *B_rand = (unsigned int *) malloc(
			B_size * sizeof(unsigned int));

	for (j = 0; j < (reps / size) + 1; j++) {
		//start();
		int k;
		for (k = 0; k < A_size; k++)
			A_rand[k] = genrand_int32() % max;
		for (k = 0; k < B_size; k++) 
			B_rand[k] = genrand_int32() % max;
		//stop();
		//printf("r:%ld\t", report());
		//rand_total_time += report();

		//start();
		qsort(A_rand, A_size, sizeof(unsigned int), compare_uints);
		qsort(B_rand, B_size, sizeof(unsigned int), compare_uints);
		//stop();
		//printf("s:%ld\t", report());
		//sort_total_time += report();

		//start();
		unsigned int r = count_intersections_bsearch(
				A_rand, A_len, A_size,
				B_rand, B_len, B_size );
		//stop();

		//printf("i:%ld\t", report());
		//intersect_total_time += report();
		//
		//printf("i:%d\tr:%u\n", rank, r);

		if (r >= O)
			++x;
	}

	if (rank == 0) {
		int seen = 0, buf;
		printf("%d\tgot:%u\tfrom:%d\n", rank, x, rank);
		while (seen < (size - 1)) {
			MPI_Recv(&buf, 1, MPI_UNSIGNED, MPI_ANY_SOURCE, 0,
			MPI_COMM_WORLD, &status);
			++seen;
			x += buf;
			printf("%d\tgot:%u\tfrom:%d\t%d\n", rank, buf, status.MPI_SOURCE,
					x);
		}

		double p = ((double)x + 1) / ((double)reps + 1);
		printf("O:%d\tp:%f\n", O, p);
	} else {
		MPI_Send(&x, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD);
	}

	MPI_Finalize();

	return 0;
}
