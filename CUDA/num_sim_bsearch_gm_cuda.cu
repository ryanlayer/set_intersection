#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cutil.h>
#include <sys/time.h>


#include "../lib/bed.h"
#include "../lib/set_intersect.h"
#include "radixsort.h"
//#include "gpu.hpp"
#include "random.hpp"
#include "../lib/timer.h"

#include "set_intersect_cuda.h"

int main(int argc, char *argv[]) {


	if (argc < 6) {
			fprintf(stderr, "usage: %s <u> <a> <b> <N> "
					"<intersect N> <sum N> <device>\n"
					"e.g., order U.bed A.bed B.bed 1 1024 0\n",
					argv[0]);
			return 1;
	}

	int chrom_num = 24;

	CUDA_SAFE_CALL( cudaSetDevice( atoi(argv[7] ) ) );
	//CUDA_SAFE_CALL( cudaFree(NULL) );


	/***********************REPLACE WITH INPUT FILE************************/	
	char *chrom_names[] = {
		"chr1",  "chr2",  "chr3", "chr4",  "chr5",  "chr6",  "chr7", "chr8",
		"chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
		"chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX",  "chrY"
	};
	/**********************************************************************/	

	struct chr_list *U, *A, *B;

	char *U_file = argv[1], *A_file = argv[2], *B_file = argv[3];
	int reps = atoi(argv[4]);
	int inter_threads = atoi(argv[5]);
	int sum_threads = atoi(argv[6]);

	if	( ( chr_list_from_bed_file(&U, chrom_names, chrom_num, U_file) == 1) ||
		  ( chr_list_from_bed_file(&A, chrom_names, chrom_num, A_file) == 1) ||
		  ( chr_list_from_bed_file(&B, chrom_names, chrom_num, B_file) == 1) ) {
		fprintf(stderr, "Error parsing bed files.\n");
		return 1;
	}

	unsigned int max = add_offsets(U, chrom_num);

	trim(U, A, chrom_num);
	trim(U, B, chrom_num);

	int A_size, B_size, U_size;

	struct bed_line *U_array, *A_array, *B_array;

	U_size = chr_array_from_list(U, &U_array, chrom_num);
	A_size = chr_array_from_list(A, &A_array, chrom_num);
	B_size = chr_array_from_list(B, &B_array, chrom_num);

	start();

	int block_size = 256;
	dim3 dimBlock(block_size);

	int grid_size = ( A_size + block_size - 1) / (block_size * 1);
	dim3 dimGridSearch( grid_size );

	cudaError_t err;
	dim3 dimGridAR( ceil(float(A_size)/float(dimBlock.x)));
	dim3 dimGridBR( ceil(float(B_size)/float(dimBlock.x)));

	unsigned int *A_start_h =
			(unsigned int *) malloc( (A_size) * sizeof(unsigned int));
	unsigned int *A_len_h =
			(unsigned int *) malloc( (A_size) * sizeof(unsigned int));
	unsigned int *B_start_h =
			(unsigned int *) malloc( (B_size) * sizeof(unsigned int));
	//unsigned int *B_end_h =
			//(unsigned int *) malloc( (B_size) * sizeof(unsigned int));
	unsigned int *B_len_h =
			(unsigned int *) malloc( (B_size) * sizeof(unsigned int));

    /*
	 * In CUDA we can sort key value pairs, 
	 * the key can be the offset, and the value can be the length
	 */
	set_start_len( U_array, U_size, A_array, A_start_h, A_len_h, A_size );
	//set_start_end( U_array, U_size, B_array, B_start_h, B_end_h, B_size );
	set_start_len( U_array, U_size, B_array, B_start_h, B_len_h, B_size );

	// Move A and B to deviceB
	unsigned int *A_start_d, *A_len_d, *B_start_d, *B_end_d, *B_len_d;
	cudaMalloc((void **)&A_start_d, (A_size)*sizeof(unsigned int));
	cudaMalloc((void **)&A_len_d, (A_size)*sizeof(unsigned int));
	cudaMalloc((void **)&B_start_d, (B_size)*sizeof(unsigned int));
	cudaMalloc((void **)&B_len_d, (B_size)*sizeof(unsigned int));
	cudaMalloc((void **)&B_end_d, (B_size)*sizeof(unsigned int));

	cudaMemcpy(A_start_d, A_start_h, (A_size) * sizeof(unsigned int), 
			cudaMemcpyHostToDevice);
	cudaMemcpy(A_len_d, A_len_h, (A_size) * sizeof(unsigned int),
			cudaMemcpyHostToDevice);
	cudaMemcpy(B_start_d, B_start_h, (B_size) * sizeof(unsigned int), 
			cudaMemcpyHostToDevice);
	cudaMemcpy(B_len_d, B_len_h, (B_size) * sizeof(unsigned int),
			cudaMemcpyHostToDevice);

	set_end <<<dimGridBR, dimBlock>>> (B_start_d, B_end_d, B_len_d, B_size);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "Set end B: %s.\n", cudaGetErrorString( err) );


	// R will hold the results of the intersection, for each interval A[i],
	// R[i] will be the number of intervals in B that A[i] intersects,
	unsigned int *R_d;
	cudaMalloc((void **)&R_d, (A_size)*sizeof(unsigned int));
	//unsigned long memup_time = report();


	// Sort A
	nvRadixSort::RadixSort radixsortA_start(A_size, true);
	/*
	radixsortA_start.sort((unsigned int*)A_start_d, (unsigned int*)A_len_d, 
			A_size, 32);
	*/
	// Sort B by start
	nvRadixSort::RadixSort radixsortB_start(B_size, true);
	radixsortB_start.sort((unsigned int*)B_start_d, 0, B_size, 32);

	// Sort B by end
	nvRadixSort::RadixSort radixsortB_end(B_size, true);
	radixsortB_end.sort((unsigned int*)B_end_d, 0, B_size, 32);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "Sort: %s.\n", cudaGetErrorString( err) );

	count_bsearch_cuda <<<dimGridSearch, dimBlock >>> (
			A_start_d, A_len_d, A_size,
			B_start_d, B_end_d, B_size,
			R_d,
			1);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "Search: %s.\n", cudaGetErrorString( err) );

	parallel_sum(R_d, block_size, A_size, sum_threads);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "Sum: %s.\n", cudaGetErrorString( err) );

	unsigned int O;
	cudaMemcpy(&O, R_d, 1 * sizeof(unsigned int), cudaMemcpyDeviceToHost);

	srand(time(NULL));	

	RNG_rand48 A_start_r(rand());
	RNG_rand48 B_start_r(rand());

	unsigned int *A_start_r_d, *B_start_r_d;
	unsigned int R, t=0;

	int i, r = 0;
	for  (i = 0; i < reps; i ++) {
		A_start_r.generate(A_size);

		cudaThreadSynchronize();
		err = cudaGetLastError();
		if(err != cudaSuccess)
			fprintf(stderr, "Rand A: %s.\n", cudaGetErrorString( err) );

		B_start_r.generate(B_size);

		cudaThreadSynchronize();
		err = cudaGetLastError();
		if(err != cudaSuccess)
			fprintf(stderr, "Rand B: %s.\n", cudaGetErrorString( err) );

		A_start_r_d = (unsigned int *)A_start_r.get_random_numbers();
		B_start_r_d = (unsigned int *)B_start_r.get_random_numbers();

		normalize_rand <<<dimGridAR, dimBlock>>> (A_start_r_d, max, A_size);

		cudaThreadSynchronize();
		err = cudaGetLastError();
		if(err != cudaSuccess)
			fprintf(stderr, "Norm A: %s.\n", cudaGetErrorString( err) );

		normalize_rand <<<dimGridBR, dimBlock>>> (B_start_r_d, max, B_size);

		cudaThreadSynchronize();
		err = cudaGetLastError();
		if(err != cudaSuccess)
			fprintf(stderr, "Norm B: %s.\n", cudaGetErrorString( err) );

		radixsortA_start.sort((unsigned int*)A_start_r_d, 0, A_size, 32);

		cudaThreadSynchronize();
		err = cudaGetLastError();
		if(err != cudaSuccess)
			fprintf(stderr, "Sort start A: %s.\n", cudaGetErrorString( err) );

		radixsortB_start.sort((unsigned int*)B_start_r_d, 0, B_size, 32);

		cudaThreadSynchronize();
		err = cudaGetLastError();
		if(err != cudaSuccess)
			fprintf(stderr, "Sort start B: %s.\n", cudaGetErrorString( err) );

		set_end <<<dimGridBR, dimBlock>>> (B_start_r_d, B_end_d, B_len_d, 
				B_size);

		cudaThreadSynchronize();
		err = cudaGetLastError();
		if(err != cudaSuccess)
			fprintf(stderr, "Set end B: %s.\n", cudaGetErrorString( err) );

		radixsortB_end.sort((unsigned int*)B_end_d, 0, B_size, 32);

		cudaThreadSynchronize();
		err = cudaGetLastError();
		if(err != cudaSuccess)
			fprintf(stderr, "Sort start B: %s.\n", cudaGetErrorString( err) );
		/*

		unsigned int *A_start_t = (unsigned int *)
				malloc(A_size * sizeof(unsigned int));
		unsigned int *B_start_t = (unsigned int *)
				malloc(B_size * sizeof(unsigned int));
		unsigned int *B_end_t = (unsigned int *)
				malloc(B_size * sizeof(unsigned int));
		unsigned int *B_len_t = (unsigned int *)
				malloc(B_size * sizeof(unsigned int));

		cudaMemcpy(A_start_t, A_start_r_d, A_size * sizeof(unsigned int),
				cudaMemcpyDeviceToHost);
		cudaMemcpy(B_start_t, B_start_r_d, B_size * sizeof(unsigned int),
				cudaMemcpyDeviceToHost);
		cudaMemcpy(B_end_t, B_end_d, B_size * sizeof(unsigned int), 
				cudaMemcpyDeviceToHost);
		cudaMemcpy(B_len_t, B_len_d, B_size * sizeof(unsigned int), 
				cudaMemcpyDeviceToHost);

		int z;
		for (z = 0; (z < B_size) && (z < A_size); z++)
			printf("s:%u\ts:%u\tl:%u\te:%u\n",
					A_start_t[z],B_start_t[z], B_end_t[z], B_len_t[z]);



		*/


		count_bsearch_cuda <<<dimGridSearch, dimBlock >>> (
				A_start_r_d, A_len_d, A_size,
				B_start_r_d, B_end_d, B_size,
				R_d,
				1);

		cudaThreadSynchronize();
		err = cudaGetLastError();
		if(err != cudaSuccess)
			fprintf(stderr, "Search: %s.\n", cudaGetErrorString( err) );

		parallel_sum(R_d, block_size, A_size, sum_threads);

		cudaThreadSynchronize();
		err = cudaGetLastError();
		if(err != cudaSuccess)
			fprintf(stderr, "Sum: %s.\n", cudaGetErrorString( err) );

		cudaMemcpy(&R, R_d, 1 * sizeof(unsigned int), cudaMemcpyDeviceToHost);
		
		//fprintf(stderr, "%u\n", R);
		t+=R;
		if (R >= O)
			++r;
	}

	stop();

	double p = ( (double)(r + 1) ) / ( (double)(reps + 1) );
	fprintf(stderr,"O:%d\tp:%G\tE:%G\n", O, p,(double)t/(double)reps);
	printf("%d,%d,%d\tt:%lu\n", A_size, B_size, A_size + B_size, report());

	/*

	double  rand_avg_time = ( (double) rand_total_time) / reps,
			sort_avg_time = ( (double) sort_total_time) / reps,
			intersect_avg_time = ( (double)  intersect_total_time) / reps;

	double total_avg_time = rand_avg_time + sort_avg_time + intersect_avg_time;

	double  rand_prop_time = rand_avg_time/total_avg_time,
			sort_prop_time = sort_avg_time/total_avg_time,
			intersect_prop_time = intersect_avg_time/total_avg_time;

	printf("t:%G\tr:%G,%G\ts:%G,%G\ti:%G,%G\n", 
			total_avg_time,
			rand_avg_time, rand_prop_time,
			sort_avg_time, sort_prop_time,
			intersect_avg_time, intersect_prop_time);

	printf("%d,%d,%d\tt:%G\tr:%G,%G\ts:%G,%G\ti:%G,%G\n", 
			A_size,
			B_size,
			A_size + B_size,
			total_avg_time,
			rand_avg_time,
			rand_prop_time,
			sort_avg_time,
			sort_prop_time,
			intersect_avg_time,
			intersect_prop_time);
	*/

	cudaFree(A_start_d);
	cudaFree(A_len_d);
	cudaFree(B_start_d);
	cudaFree(B_len_d);
	cudaFree(B_end_d);
	cudaFree(R_d);

	return 0;
}
