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


	if (argc < 7) {
		fprintf(stderr, "usage: %s <u> <a> <b> "
				"<inter N> <sum N> <chuck size> <device>\n"
				"e.g., order U.bed A.bed B.bed 1 1024 1000 0\n",
				argv[0]);
		return 1;
	}

	int chrom_num = 24;

	//CUDA_SAFE_CALL( cudaFree(NULL) );
	CUDA_SAFE_CALL( cudaSetDevice( atoi(argv[7] ) ) );


	/***********************REPLACE WITH INPUT FILE************************/	
	char *chrom_names[] = {
		"chr1",  "chr2",  "chr3", "chr4",  "chr5",  "chr6",  "chr7", "chr8",
		"chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
		"chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX",  "chrY"
	};
	/**********************************************************************/	

	struct chr_list *U, *A;

	char *U_file_name = argv[1], *A_file_name = argv[2], *B_file_name = argv[3];
	int inter_threads = atoi(argv[4]);
	int sum_threads = atoi(argv[5]);
	int chunk_size = atoi(argv[6]);


	if	( ( chr_list_from_bed_file(&U, chrom_names, chrom_num,
					U_file_name) == 1) ||
		  ( chr_list_from_bed_file(&A, chrom_names, chrom_num,
					A_file_name) == 1) ) {
		fprintf(stderr, "Error parsing bed files.\n");
		return 1;
	}

	FILE *B_file = fopen(B_file_name, "r");

	if (B_file == NULL) {
		fprintf(stderr, "Could not open file:%s\n", B_file_name);
		return 1;
	}

	unsigned int max = add_offsets(U, chrom_num);

	trim(U, A, chrom_num);

	int A_size, U_size;

	struct bed_line *U_array, *A_array;

	U_size = chr_array_from_list(U, &U_array, chrom_num);
	A_size = chr_array_from_list(A, &A_array, chrom_num);

	start();
	unsigned int *A_start_h = 
		(unsigned int *) malloc( (A_size) * sizeof(unsigned int));
	unsigned int *A_len_h = 
		(unsigned int *) malloc( (A_size) * sizeof(unsigned int));

	unsigned int *B_start_h = 
		(unsigned int *) malloc( (chunk_size) * sizeof(unsigned int));
	unsigned int *B_end_h = 
		(unsigned int *) malloc( (chunk_size) * sizeof(unsigned int));

	map_to_start_len_array(A_start_h, A_len_h, A_array, A_size, U_array,
			U_size);

	// Move A and B to deviceB
	unsigned int *A_start_d, *A_len_d, *B_start_d, *B_end_d;
	cudaMalloc((void **)&A_start_d, (A_size)*sizeof(unsigned int));
	cudaMalloc((void **)&A_len_d, (A_size)*sizeof(unsigned int));
	cudaMalloc((void **)&B_start_d, (chunk_size)*sizeof(unsigned int));
	cudaMalloc((void **)&B_end_d, (chunk_size)*sizeof(unsigned int));

	cudaMemcpy(A_start_d, A_start_h, (A_size) * sizeof(unsigned int), 
			cudaMemcpyHostToDevice);
	cudaMemcpy(A_len_d, A_len_h, (A_size) * sizeof(unsigned int),
			cudaMemcpyHostToDevice);
	/*
	cudaMemcpy(B_start_d, B_start_h, (B_size) * sizeof(unsigned int), 
			cudaMemcpyHostToDevice);
	cudaMemcpy(B_end_d, B_end_h, (B_size) * sizeof(unsigned int),
			cudaMemcpyHostToDevice);
	*/

	// R will hold the results of the intersection, for each interval A[i],
	// R[i] will be the number of intervals in B that A[i] intersects,
	unsigned int *R_d;
	cudaMalloc((void **)&R_d, (A_size)*sizeof(unsigned int));
	cudaMemset(R_d, 0, (A_size)*sizeof(unsigned int));

	int block_size = 256;
	dim3 dimBlock(block_size);

	int grid_size = ( A_size + block_size - 1) / (block_size * 1);
	dim3 dimGridSearch( grid_size );

	cudaError_t err;

	// Sort A
	/*
	nvRadixSort::RadixSort radixsortA(A_size, false);
	radixsortA.sort((unsigned int*)A_start_d, (unsigned int*)A_len_d, 
			A_size, 32);
	*/
	// Sort B by start
	nvRadixSort::RadixSort radixsortB_start(chunk_size, true);
	//radixsortB_start.sort((unsigned int*)B_start_d, 0, B_size, 32);

	// Sort B by end
	nvRadixSort::RadixSort radixsortB_end(chunk_size, true);
	//radixsortB_end.sort((unsigned int*)B_end_d, 0, B_size, 32);

	unsigned int B_curr_size;
	unsigned int line = 0;
	int x=1;
	while ( map_start_end_from_file( B_file, B_start_h, B_end_h, 
									 chunk_size, &B_curr_size,
									 U_array, U_size) ) {
		line += B_curr_size;
		//printf("%d\n",x++);
		cudaMemcpy(B_start_d, B_start_h, (B_curr_size) * sizeof(unsigned int), 
				cudaMemcpyHostToDevice);
		cudaMemcpy(B_end_d, B_end_h, (B_curr_size) * sizeof(unsigned int),
			cudaMemcpyHostToDevice);

		if (B_curr_size < chunk_size) {
			nvRadixSort::RadixSort last_radixsortB_start(B_curr_size, true);
			nvRadixSort::RadixSort last_radixsortB_end(B_curr_size, true);
			last_radixsortB_start.sort((unsigned int*)B_start_d, 0, 
					B_curr_size, 32);
			last_radixsortB_end.sort((unsigned int*)B_end_d, 0,
					B_curr_size, 32);
		} else {
			radixsortB_start.sort((unsigned int*)B_start_d, 0, B_curr_size, 32);
			radixsortB_end.sort((unsigned int*)B_end_d, 0, B_curr_size, 32);
		}

		big_count_bsearch_cuda <<<dimGridSearch, dimBlock >>> (
				A_start_d, A_len_d, A_size,
				B_start_d, B_end_d, B_curr_size,
				R_d,
				1);

		cudaThreadSynchronize();
		err = cudaGetLastError();
		if(err != cudaSuccess)
			fprintf(stderr, "Sort: %s.\n", cudaGetErrorString( err) );
	}

	parallel_sum(R_d, block_size, A_size, sum_threads);

	unsigned int O;
	cudaMemcpy(&O, R_d, 1 * sizeof(unsigned int), cudaMemcpyDeviceToHost);
	stop();
	unsigned long total = report();

	printf("%d,%d,%d\tO:%d\t\tt:%ldc:%d\n",
				A_size, line, A_size + line, O, total, chunk_size);

	cudaFree(A_start_d);
	cudaFree(B_start_d);
	cudaFree(R_d);

	return 0;
}
