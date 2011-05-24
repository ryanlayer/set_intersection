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
	cudaError_t err;

	RNG_rand48 D_r(seed);
	D_r.generate(D_size);
	unsigned int *D_d = (unsigned int *)D_r.get_random_numbers();

	RNG_rand48 Q_r(seed);
	Q_r.generate(Q_size);
	unsigned int *Q_d = (unsigned int *)Q_r.get_random_numbers();
	
	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "rand errors: %s.\n", cudaGetErrorString( err) );


	nvRadixSort::RadixSort sort_D_d(D_size, true);
	start();
	sort_D_d.sort((unsigned int*)D_d, 0, D_size, 32);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "binary_search_i: %s.\n", cudaGetErrorString( err) );

	stop();
	printf("x:%lu\n", report());



	unsigned int *D_h = (unsigned int *)malloc( D_size * sizeof(unsigned int));
	cudaMemcpy(D_h, D_d, (D_size) * sizeof(unsigned int),
			cudaMemcpyDeviceToHost);

	int block_size = 256;
	dim3 dimBlock(block_size);

	int grid_size = ( Q_size + block_size - 1) / (block_size * 1);
	dim3 dimGrid( grid_size );

	unsigned int *R_d;
	cudaMalloc((void **)&R_d, (Q_size)*sizeof(unsigned int));

	start();
	binary_search_n <<<dimGrid, dimBlock>>> (D_d, D_size, Q_d, Q_size, R_d);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "binary_search_n: %s.\n", cudaGetErrorString( err) );

	stop();
	printf("n:%lu\n", report());

	unsigned int *R_h_n = (unsigned int *)malloc( I_size * sizeof(unsigned int));
	cudaMemcpy(R_h_n, R_d, (I_size) * sizeof(unsigned int), 
			cudaMemcpyDeviceToHost);

	start();
	binary_search_i <<< dimGrid, dimBlock, I_size * sizeof(unsigned int) >>> (
			D_d, D_size, Q_d, Q_size, R_d, I_size);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "binary_search_i: %s.\n", cudaGetErrorString( err) );

	stop();
	printf("i:%lu\n", report());

	unsigned int *R_h_i = (unsigned int *)malloc( I_size * sizeof(unsigned int));
	cudaMemcpy(R_h_i, R_d, (I_size) * sizeof(unsigned int), 
			cudaMemcpyDeviceToHost);

	int index_grid_size = ( I_size + block_size - 1) / (block_size * 1);
	dim3 index_dimGrid( index_grid_size );

	unsigned int *I_d;
	cudaMalloc((void **)&I_d, (I_size)*sizeof(unsigned int));

	start();
	gen_index <<<index_dimGrid,dimBlock>>> ( D_d, D_size, I_d, I_size);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "index: %s.\n", cudaGetErrorString( err) );

	stop();

	printf("p:%lu\n", report());

	start();
	binary_search_p <<< dimGrid, dimBlock, I_size * sizeof(unsigned int) >>> (
			D_d, D_size, Q_d, Q_size, R_d, I_d, I_size);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "binary_search_i: %s.\n", cudaGetErrorString( err) );

	stop();
	printf("s:%lu\n", report());




	return 0;
}
