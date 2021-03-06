#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cutil.h>
#include <sys/time.h>
#include <time.h>


#include "bed.h"
#include "set_intersect.h"
#include "radixsort.h"
#include "gpu.hpp"
#include "random.hpp"

#include "order_kernel.cu"

void set_start_len( struct bed_line *U_array,
					int U_size,
					struct bed_line *A_array,
					unsigned int *A_key_h,
					unsigned int *A_val_h,
					int A_size )
{
	int i, j, k = 0;
	for (i = 0; i < A_size; i++) {
		int start = -1, offset = -1;
		for (j = 0; j < U_size; j++) {
			if ( ( U_array[j].chr == A_array[i].chr) &&
				 ( U_array[j].start <= A_array[i].end) &&
				 ( U_array[j].end >= A_array[i].start) ) {
				start = U_array[j].start;
				offset = U_array[j].offset;
				break;
			}
		}
		A_key_h[k] = A_array[i].start - start + offset;
		A_val_h[k] = A_array[i].end -A_array[i].start;
		++k;
	}
}

void parallel_sum( int *R_d,
				   int block_size,
				   int A_size,
				   int n)
{
	unsigned int left = A_size;
	//int n = 1024;
	while (left > 1) {

		dim3 dimGridR( left / (block_size * n) + 1);
		//dim3 dimGridR( left / blocksize  + 1);
		dim3 dimBlockR( block_size );
		size_t sm_size = dimBlockR.x * sizeof(int); 

		my_reduce <<<dimGridR, dimBlockR, sm_size>>> (R_d, left, n);

		cudaThreadSynchronize();
		cudaError_t err;
		err = cudaGetLastError();
		if(err != cudaSuccess)
			fprintf(stderr, "My Reduce0: %s.\n", cudaGetErrorString( err) );

		left = dimGridR.x;
	}
}

int main(int argc, char *argv[]) {

	if (argc < 4) {
		fprintf(stderr, "usage: order <u> <a> <b> <N>\n");
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
	unsigned int *A_key_d, *A_val_d, *B_key_d, *B_val_d;
	cudaMalloc((void **)&A_key_d, (A_size)*sizeof(unsigned int));
	cudaMalloc((void **)&A_val_d, (A_size)*sizeof(unsigned int));
	cudaMalloc((void **)&B_key_d, (B_size)*sizeof(unsigned int));
	cudaMalloc((void **)&B_val_d, (B_size)*sizeof(unsigned int));

	cudaMemcpy(A_key_d, A_key_h, (A_size) * sizeof(unsigned int), 
			cudaMemcpyHostToDevice);
	cudaMemcpy(A_val_d, A_val_h, (A_size) * sizeof(unsigned int),
			cudaMemcpyHostToDevice);
	cudaMemcpy(B_key_d, B_key_h, (B_size) * sizeof(unsigned int), 
			cudaMemcpyHostToDevice);
	cudaMemcpy(B_val_d, B_val_h, (B_size) * sizeof(unsigned int),
			cudaMemcpyHostToDevice);

	// Sort A
	nvRadixSort::RadixSort radixsortA(A_size, false);
	radixsortA.sort((unsigned int*)A_key_d, (unsigned int*)A_val_d, 
			A_size, 32);

	// Sort B
	nvRadixSort::RadixSort radixsortB(B_size, false);
	radixsortB.sort((unsigned int*)B_key_d, (unsigned int*)B_val_d, 
			B_size, 32);

	cudaThreadSynchronize();
	cudaError_t err;
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "Rand A: %s.\n", cudaGetErrorString( err) );


	int block_size = 256;
	dim3 dimBlock(block_size);

	// I want a thread per element in A
	dim3 dimGridA( ceil(float(A_size)/float(dimBlock.x)));

	// R will hold the results of the intersection, for each interval A[i],
	// R[i] will be the number of intervals in B that A[i] intersects,
	int *R_d;
	cudaMalloc((void **)&R_d, (A_size)*sizeof(int));

	// *_key_d holds the start position, and *_val_d holds the length,
	// the end position is *_key_d + *_val_d
	intersection_b_search <<<dimGridA, dimBlock>>> (
			A_key_d, A_val_d, A_size,
			B_key_d, B_val_d, B_size,
			R_d, 1);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "intersect search: %s.\n", cudaGetErrorString( err) );

	int n = 1024;
	parallel_sum( R_d, block_size, A_size, n);

	int O;
	cudaMemcpy(&O, R_d, 1 * sizeof(int), cudaMemcpyDeviceToHost);
	printf("O:%d\n", O);

	cudaFree(A_key_d);
	cudaFree(B_key_d);

	/************** SIMULATION *********************/
	srand(time(NULL));	

	RNG_rand48 A_r(rand());
	RNG_rand48 B_r(rand());
	dim3 dimGridAR( ceil(float(A_size)/float(dimBlock.x)));
	dim3 dimGridBR( ceil(float(B_size)/float(dimBlock.x)));
	nvRadixSort::RadixSort radixsortAR(A_size, true);
	nvRadixSort::RadixSort radixsortBR(B_size, true);

	int i, R, r = 0;

	for  (i = 0; i < reps; i ++) {
		A_r.generate(A_size);

		cudaThreadSynchronize();
		err = cudaGetLastError();
		if(err != cudaSuccess)
			fprintf(stderr, "Rand A: %s.\n", cudaGetErrorString( err) );


		B_r.generate(B_size);

		cudaThreadSynchronize();
		err = cudaGetLastError();
		if(err != cudaSuccess)
			fprintf(stderr, "Rand B: %s.\n", cudaGetErrorString( err) );


		unsigned int *A_r_d = (unsigned int *)A_r.get_random_numbers();
		unsigned int *B_r_d = (unsigned int *)B_r.get_random_numbers();

		normalize_rand <<<dimGridAR, dimBlock>>> (A_r_d, max, A_size);

		cudaThreadSynchronize();
		err = cudaGetLastError();
		if(err != cudaSuccess)
			fprintf(stderr, "Norm A: %s.\n", cudaGetErrorString( err) );

		normalize_rand <<<dimGridBR, dimBlock>>> (B_r_d, max, B_size);

		cudaThreadSynchronize();
		err = cudaGetLastError();
		if(err != cudaSuccess)
			fprintf(stderr, "Norm B: %s.\n", cudaGetErrorString( err) );

		radixsortAR.sort((unsigned int*)A_r_d, 0, A_size, 32);

		cudaThreadSynchronize();
		err = cudaGetLastError();
		if(err != cudaSuccess)
			fprintf(stderr, "Sort A: %s.\n", cudaGetErrorString( err) );

		radixsortBR.sort((unsigned int*)B_r_d, 0, B_size, 32);

		cudaThreadSynchronize();
		err = cudaGetLastError();
		if(err != cudaSuccess)
			fprintf(stderr, "Sort B: %s.\n", cudaGetErrorString( err) );

		intersection_b_search <<<dimGridA, dimBlock>>> (
				A_r_d, A_val_d, A_size,
				B_r_d, B_val_d, B_size,
				R_d);

		cudaThreadSynchronize();
		err = cudaGetLastError();
		if(err != cudaSuccess)
			fprintf(stderr, "intersect search: %s.\n", 
					cudaGetErrorString(err) );

		parallel_sum( R_d, block_size, A_size, n);

		cudaMemcpy(&R, R_d, 1 * sizeof(int), cudaMemcpyDeviceToHost);

		if (R >= O)
			++r;
	}

	printf("p:%f\n", (double)(r + 1) / (double)(reps + 1));

	return 0;
}
