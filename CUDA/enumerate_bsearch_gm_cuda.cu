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

	cudaError_t err;

	cudaFree(NULL);

	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "cudaFree: %s.\n", cudaGetErrorString( err) );

	if (argc < 6) {
		fprintf(stderr, "usage: order <u> <a> <b> <reps> <inter N> <sum N>\n"
						"e.g., order U.bed A.bed B.bed 10000 1 1024\n");
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

	unsigned int *A_key_h = 
		(unsigned int *) malloc( (A_size) * sizeof(unsigned int));
	unsigned int *A_val_h = 
		(unsigned int *) malloc( (A_size) * sizeof(unsigned int));

	unsigned int *B_key_h = 
		(unsigned int *) malloc( (B_size) * sizeof(unsigned int));
	unsigned int *B_val_h = 
		(unsigned int *) malloc( (B_size) * sizeof(unsigned int));


	/*
	 * In CUDA we can sort key value pairs, 
	 * the key can be the offset, and the value can be the length
	 */
	set_start_len( U_array, U_size,
				   A_array, A_key_h, A_val_h, A_size );

	set_start_len( U_array, U_size,
				   B_array, B_key_h, B_val_h, B_size );

	// Move A and B to deviceB
	start();
	unsigned int *A_key_d, *A_val_d, *B_key_d, *B_val_d;
	CUDA_SAFE_CALL( 
		cudaMalloc((void **)&A_key_d, (A_size)*sizeof(unsigned int)));
	CUDA_SAFE_CALL(
		cudaMalloc((void **)&A_val_d, (A_size)*sizeof(unsigned int)));
	CUDA_SAFE_CALL(
		cudaMalloc((void **)&B_key_d, (B_size)*sizeof(unsigned int)));
	CUDA_SAFE_CALL(
		cudaMalloc((void **)&B_val_d, (B_size)*sizeof(unsigned int)));

	CUDA_SAFE_CALL( 
		cudaMemcpy(A_key_d, A_key_h, (A_size) * sizeof(unsigned int), 
			cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL( 
		cudaMemcpy(A_val_d, A_val_h, (A_size) * sizeof(unsigned int),
			cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL( 
		cudaMemcpy(B_key_d, B_key_h, (B_size) * sizeof(unsigned int), 
			cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL( 
		cudaMemcpy(B_val_d, B_val_h, (B_size) * sizeof(unsigned int),
			cudaMemcpyHostToDevice));
	stop();
	unsigned long memup_time = report();

	int block_size = 256;
	dim3 dimBlock(block_size);

	unsigned int *P1_d,*P2_d;
	start();
	CUDA_SAFE_CALL(
			cudaMalloc((void **)&P1_d, (A_size + B_size)*sizeof(unsigned int)));
	CUDA_SAFE_CALL(cudaMemset(P1_d, 1, (A_size + B_size)*sizeof(unsigned int)));
	CUDA_SAFE_CALL(
			cudaMalloc((void **)&P2_d, (A_size + B_size)*sizeof(unsigned int)));
	CUDA_SAFE_CALL(cudaMemset(P2_d, 1, (A_size + B_size)*sizeof(unsigned int)));
	stop();
	cudaThreadSynchronize();
	unsigned long memset_time = report();


	// *_key_d holds the start position, and *_val_d holds the length,
	// the end position is *_key_d + *_val_d
	//
	// Each thread will search |reps| items in A, we will keep the blocksize
	// fixed at 256, but we will need to adjust the grid size 
	
	int grid_size = ( A_size + block_size - 1) / (block_size * 1);
	dim3 dimGridSearch( grid_size );

	// Sort A
	start();
	nvRadixSort::RadixSort radixsortA(A_size, false);
	radixsortA.sort((unsigned int*)A_key_d, (unsigned int*)A_val_d, 
			A_size, 32);

	// Sort B
	nvRadixSort::RadixSort radixsortB(B_size, false);
	radixsortB.sort((unsigned int*)B_key_d, (unsigned int*)B_val_d, 
			B_size, 32);
	cudaThreadSynchronize();
	stop();

	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "Sort: %s.\n", cudaGetErrorString( err) );
	unsigned long sort_time = report();

	unsigned int *P1_h = (unsigned int *) calloc( 
			(A_size + B_size), sizeof(unsigned int));
	unsigned int *P2_h = (unsigned int *) calloc( 
			(A_size + B_size), sizeof(unsigned int));

	start();
	enumerate_b_search_gm <<<dimGridSearch, 
							 dimBlock >>> ( A_key_d, A_val_d, A_size,
										   B_key_d, B_val_d, B_size,
										   P1_d, P2_d, 1);

	cudaThreadSynchronize();
	stop();
	unsigned long search_time = report();

	unsigned int O;
	start();
	cudaMemcpy(P1_h, P1_d, (A_size + B_size)*sizeof(unsigned int), 
			cudaMemcpyDeviceToHost);
	cudaMemcpy(P2_h, P2_d, (A_size + B_size)*sizeof(unsigned int), 
			cudaMemcpyDeviceToHost);
	stop();
	unsigned long memdown_time = report();

	unsigned long total = memup_time + memset_time + 
			sort_time + search_time + memdown_time;

	/*
	int i;
	for (i = 0; i < A_size + B_size; i++)
		if ( (P1_h[i] != 16843009) && (P2_h[i] != 16843009) )
			printf("%ld\t%ld\n", P1_h[i] , P2_h[i]);
	*/


	printf("t:%ld\t"
			"up:%ld,%G\t"
			"set:%ld,%G\t"
			"sort:%ld,%G\t"
			"search:%ld,%G\t"
			"down:%ld,%G\n",
			total,
			memup_time, (double)memup_time / (double)total,
			memset_time, (double)memset_time / (double)total,
			sort_time, (double)sort_time / (double)total,
			search_time, (double)search_time/ (double)total,
			memdown_time, (double)memdown_time / (double)total
		  );



	cudaFree(P1_d);
	cudaFree(P2_d);
	cudaFree(A_key_d);
	cudaFree(B_key_d);

	return 0;
}
