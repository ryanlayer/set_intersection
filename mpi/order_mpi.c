#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <mpi.h>

#include "../lib/bed.h"
#include "../lib/set_intersect.h"

#define MIN(a,b) ((a)>(b)?(b):(a))

int main(int argc, char *argv[]) {

	//MPI_Status status;
	int rank, size;

	MPI_Init (&argc, &argv);    /* starts MPI */
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);  /* get current process id */
	MPI_Comm_size (MPI_COMM_WORLD, &size);  /* get number of processes */

	if (argc < 4) {
		fprintf(stderr, "usage: order <u> <a> <b> <N>\n");
		MPI_Finalize();
		return 1;
	}

	fprintf(stderr, "%d\n", rank);

	// Going to assign each thread to take N/size of the simulations.  Instead
	// of having the master send out the results of the original intersection,
	// we are going to make everyone do that.  We a thread is done, it will
	// send its results, for both each pair and for the number of pairs

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
	 */
	struct triple *AB = (struct triple *)
			malloc((2*A_size + 2*B_size)*sizeof(struct triple));
	//A and B points to AB, A to the beging and B to the interior, after A
	struct triple *A = AB;
	struct triple *B = AB + 2*A_size;

	map_intervals(A, A_array, A_size, U_array, U_size, 0 );
	map_intervals(B, B_array, B_size, U_array, U_size, 1 );

	int j;

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


	qsort(AB, 2*A_size + 2*B_size, sizeof(struct triple), compare_triple_lists);

	// find the intersecting ranks there are atmost A + B pairs
	int *pairs = (int *) malloc( 2 * (A_size + B_size) * sizeof(int));
	int num_pairs = find_intersecting_ranks(AB, A_size, B_size, pairs);

	int *A_r = (int *) malloc (A_size * sizeof(int));
	int *B_r = (int *) malloc (B_size * sizeof(int));
	int *R = (int *) calloc(num_pairs, sizeof(int));

	int r = 0;

	srand((unsigned)time(NULL) + rank);

	for (j = 0; j < (reps / size); j++) {

		// Fill and sort A and B
		for (i = 0; i < A_size; i++)
			A_r[i] = rand() % max;
		qsort(A_r, A_size, sizeof(int), compare_ints);

		for (i = 0; i < B_size; i++)
			B_r[i] = rand() % max;
		qsort(B_r, B_size, sizeof(int), compare_ints);


		//int x = check_observed_ranks(pairs, A_r, A_len, B_r, B_len, 
		check_observed_ranks(pairs, A_r, A_len, B_r, B_len, 
				num_pairs, R);

		int o = count_intersecitons_scan(A_r, A_len, A_size, B_r, B_len,
				B_size);

		if (o >= num_pairs)
			r++;
	}

	printf("%d\to:%d\tr:%d\tN:%d\n", rank, num_pairs, r, 
			(reps / size) );

	MPI_Finalize();
	return 0;
}
