#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cutil.h>
#include <sys/time.h>


#include "../lib/bed.h"
#include "../lib/set_intersect.h"
#include "radixsort.h"
#include "random.hpp"
#include "../lib/timer.h"

#include "set_intersect_cuda.h"

int main(int argc, char *argv[]) {


	if (argc < 4) {
			fprintf(stderr, "usage: %s <D size> <I size> <Q size> <seed>\n",
					argv[0]);
			return 1;
	}

	//CUDA_SAFE_CALL( cudaSetDevice( atoi(argv[7] ) ) );
	CUDA_SAFE_CALL( cudaFree(NULL) );

	int D_size = atoi(argv[1]);
	int I_size = atoi(argv[2]);
	int Q_size = atoi(argv[3]);
	int seed = atoi(argv[4]);

	RNG_rand48 D_r(seed);
	D_r.generate(D_size);
	unsigned int *D_d = (unsigned int *)D_r.get_random_numbers();

	RNG_rand48 Q_r(seed);
	Q_r.generate(Q_size);
	unsigned int *Q_d = (unsigned int *)Q_r.get_random_numbers();

	nvRadixSort::RadixSort sort_D_d(D_size, true);
	sort_D_d.sort((unsigned int*)D_d, 0, D_size, 32);

	unsigned int *R_d;
	cudaMalloc((void **)&R_d, (Q_size)*sizeof(unsigned int));

	/*
	unsigned int *D_h = (unsigned int *)malloc( D_size * sizeof(unsigned int));
	cudaMemcpy(D_h, D_d, (D_size) * sizeof(unsigned int), 
			cudaMemcpyDeviceToHost);
	int i;
	for (i = 0; i < D_size; i++)
		printf("%d\n", D_h[i]);
	*/

	int block_size = 256;
	dim3 dimBlock(block_size);

	int grid_size = ( Q_size + block_size - 1) / (block_size * 1);
	dim3 dimGrid( grid_size );

	binary_search_n <<<dimGrid, dimBlock>>> (D_d, D_size, Q_d, Q_size, R_d);

	/*
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

	set_start_len( U_array, U_size, A_array, A_start_h, A_len_h, A_size );
	//set_start_end( U_array, U_size, B_array, B_start_h, B_end_h, B_size );
	set_start_len( U_array, U_size, B_array, B_start_h, B_len_h, B_size );

	// Move A and B to deviceB
	unsigned int *A_start_d, *A_len_d, *B_start_d, *B_end_d, *B_len_d;
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


	cudaFree(A_start_d);
	cudaFree(A_len_d);
	cudaFree(B_start_d);
	cudaFree(B_len_d);
	cudaFree(B_end_d);
	cudaFree(R_d);
	*/

	return 0;
}
