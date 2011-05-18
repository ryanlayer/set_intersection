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
		fprintf(stderr, "usage: %s <u> <a> <b> "
				"<inter N> <sum N> <chuck size>\n"
				"e.g., order U.bed A.bed B.bed 1 1024 0\n",
				argv[0]);
		return 1;
	}

	int chrom_num = 24;

	CUDA_SAFE_CALL( cudaFree(NULL) );

	/***********************REPLACE WITH INPUT FILE************************/	
	char *chrom_names[] = {
		"chr1",  "chr2",  "chr3", "chr4",  "chr5",  "chr6",  "chr7", "chr8",
		"chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
		"chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX",  "chrY"
	};
	/**********************************************************************/	

	struct chr_list *U, *A, *B;

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

	unsigned int max = add_offsets(U, chrom_num);

	trim(U, A, chrom_num);

	int A_size, U_size;

	struct bed_line *U_array, *A_array;

	U_size = chr_array_from_list(U, &U_array, chrom_num);
	A_size = chr_array_from_list(A, &A_array, chrom_num);

	unsigned int *A_start_h = 
		(unsigned int *) malloc( (A_size) * sizeof(unsigned int));
	unsigned int *A_len_h = 
		(unsigned int *) malloc( (A_size) * sizeof(unsigned int));

	unsigned int *B_start_h = 
		(unsigned int *) malloc( (chunk_size) * sizeof(unsigned int));
	unsigned int *B_end_h = 
		(unsigned int *) malloc( (chunk_size) * sizeof(unsigned int));

	map_to_start_len_array(A_start, A_len, A_array, A_size, U_array, U_size);

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
	nvRadixSort::RadixSort radixsortB_start(B_size, true);
	//radixsortB_start.sort((unsigned int*)B_start_d, 0, B_size, 32);

	// Sort B by end
	nvRadixSort::RadixSort radixsortB_end(B_size, true);
	//radixsortB_end.sort((unsigned int*)B_end_d, 0, B_size, 32);

	while ( map_start_end_from_file( B_file, B_start_d, B_end_d, 
									 chunk_size, &B_curr_size,
									 U_array, U_size) ) {

		cudaMemcpy(B_start_d, B_start_h, (B_curr_size) * sizeof(unsigned int), 
				cudaMemcpyHostToDevice);
		cudaMemcpy(B_end_d, B_end_h, (B_curr_size) * sizeof(unsigned int),
			cudaMemcpyHostToDevice);

		radixsortB_start.sort((unsigned int*)B_start_d, 0, B_curr_size, 32);
		radixsortB_end.sort((unsigned int*)B_end_d, 0, B_curr_size, 32);

		count_bsearch_cuda <<<dimGridSearch, dimBlock >>> (
				A_start_d, A_len_d, A_size,
				B_start_d, B_end_d, B_size,
				R_d,
				1);
	}

	parallel_sum(R_d, block_size, A_size, sum_threads);


	unsigned int O;
	cudaMemcpy(&O, R_d, 1 * sizeof(unsigned int), cudaMemcpyDeviceToHost);

	/*
	unsigned long total = memup_time + sort_time +
			search_time + sum_time + memdown_time;
	*/

	fprintf(stderr,"O:%d\n", O);

	/*
	printf("%d,%d,%d\tT:%ld\t"
			"up:%ld,%G\t"
			"sort:%ld,%G\t"
			"search:%ld,%G\t"
			"sum:%ld,%G\t"
			"down:%ld,%G\n",
			A_size,
			B_size,
			A_size + B_size,
			total,
			memup_time, (double)memup_time / (double)total,
			sort_time, (double)sort_time / (double)total,
			search_time, (double)search_time / (double)total,
			sum_time, (double)sum_time / (double)total,
			memdown_time, (double)memdown_time / (double)total
		  );
	*/

	cudaFree(A_start_d);
	cudaFree(B_start_d);
	cudaFree(R_d);

	return 0;
}
