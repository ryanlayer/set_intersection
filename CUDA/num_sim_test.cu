#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <cuda.h>
#include <cutil.h>

#include "../lib/bed.h"
#include "../lib/set_intersect.h"
#include "../lib/mt.h"
#include "../lib/timer.h"
#include "set_intersect_cuda.h"

#define MIN(a,b) ((a)>(b)?(b):(a))

int main(int argc, char *argv[]) {
	if (argc < 4) {
		fprintf(stderr, "usage: num_sim_scan_seq <u> <a> <b> <N>\n");
		return 1;
	}

	int chrom_num = 24;
	fprintf(stderr, "*");
	cudaFree(NULL);
	fprintf(stderr, "*");

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


	int O = count_intersections_bsearch(
			A_start, A_len, A_size,
			B_start, B_len, B_size );

	init_genrand((unsigned) time(NULL));

	unsigned int *A_start_d, *A_len_d, *B_start_d, *B_len_d, *R_d;

	cudaMalloc((void **)&A_start_d, (A_size)*sizeof(unsigned int));
	cudaMalloc((void **)&A_len_d, (A_size)*sizeof(unsigned int));
	cudaMalloc((void **)&B_start_d, (B_size)*sizeof(unsigned int));
	cudaMalloc((void **)&B_len_d, (B_size)*sizeof(unsigned int));
	cudaMalloc((void **)&R_d, (A_size)*sizeof(unsigned int));

	cudaMemcpy(A_len_d, A_len, (A_size) * sizeof(unsigned int), 
			cudaMemcpyHostToDevice);
	cudaMemcpy(B_len_d, B_len, (B_size) * sizeof(unsigned int), 
			cudaMemcpyHostToDevice);

	int block_size = 256;
	dim3 dimBlock(block_size);

	int grid_size = ( A_size + block_size - 1) / (block_size * 1);
	dim3 dimGridSearch( grid_size );

	cudaError_t err;


	int j, R, x = 0;
	for (j = 0; j < reps; j++) {

		for (i = 0; i < A_size; i++)
			A_start[i] = genrand_int32();
		for (i = 0; i < B_size; i++) 
			B_start[i] = genrand_int32();

		qsort(A_start, A_size, sizeof(unsigned int), compare_uints);
		qsort(B_start, B_size, sizeof(unsigned int), compare_uints);

		int r = count_intersections_bsearch(
				A_start, A_len, A_size,
				B_start, B_len, B_size );

		int t = count_intersections_scan(
				A_start, A_len, A_size,
				B_start, B_len, B_size );

		int u = count_intersections_brute_force(
				A_start, A_len, A_size,
				B_start, B_len, B_size );

		cudaMemcpy(A_start_d, A_start, (A_size) * sizeof(unsigned int), 
				cudaMemcpyHostToDevice);
		cudaMemcpy(B_start_d, B_start, (B_size) * sizeof(unsigned int), 
				cudaMemcpyHostToDevice);

		intersection_brute_force <<<dimGridSearch, dimBlock >>>
				( A_start_d, A_len_d, A_size, B_start_d, B_len_d, B_size, R_d);

		cudaThreadSynchronize();
		err = cudaGetLastError();
		if(err != cudaSuccess)
			fprintf(stderr, "intersect search: %s.\n", 
					cudaGetErrorString(err) );

		parallel_sum( R_d, block_size, A_size, 1024);
		cudaMemcpy(&R, R_d, 1 * sizeof(unsigned int), cudaMemcpyDeviceToHost);


		printf("%d,%d,%d,%u\t", r,t,u,R);
	}
	return 0;
}
