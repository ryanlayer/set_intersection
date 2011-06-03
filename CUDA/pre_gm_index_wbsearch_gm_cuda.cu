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


	if (argc < 5) {
			fprintf(stderr, "usage: %s <D size> <I size> <Q size> <seed> <device>\n",
					argv[0]);
			return 1;
	}

	CUDA_SAFE_CALL( cudaSetDevice( atoi(argv[5] ) ) );
	//CUDA_SAFE_CALL( cudaFree(NULL) );

	int D_size = atoi(argv[1]);
	int I_size = atoi(argv[2]);
	int Q_size = atoi(argv[3]);
	int seed = atoi(argv[4]);
	cudaError_t err;

	//{{{ gen Q and D
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
	//}}}

	//{{{ sort D
	start();
	nvRadixSort::RadixSort sort_D_d(D_size, true);
	sort_D_d.sort((unsigned int*)D_d, 0, D_size, 32);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "sort d: %s.\n", cudaGetErrorString( err) );

	stop();
	unsigned long sort_d_time = report();
	//}}}

	unsigned int *D_h = (unsigned int *)malloc( D_size * sizeof(unsigned int));
	cudaMemcpy(D_h, D_d, (D_size) * sizeof(unsigned int),
			cudaMemcpyDeviceToHost);

	int block_size = 256;
	dim3 dimBlock(block_size);

	int grid_size = ( Q_size + block_size - 1) / (block_size * 1);
	dim3 dimGrid( grid_size );

	unsigned int *R_d;
	cudaMalloc((void **)&R_d, (Q_size)*sizeof(unsigned int));

	//{{{binary_search_n 
	start();
	binary_search_n <<<dimGrid, dimBlock>>> (D_d, D_size, Q_d, Q_size, R_d);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "binary_search_n: %s.\n", cudaGetErrorString( err) );

	stop();
	unsigned long search_noindex_1_time = report();

	start();
	binary_search_n <<<dimGrid, dimBlock>>> (D_d, D_size, Q_d, Q_size, R_d);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "binary_search_n: %s.\n", cudaGetErrorString( err) );

	stop();
	unsigned long search_noindex_2_time = report();
	//}}}


	unsigned int *R_h_n = (unsigned int *)malloc( I_size * sizeof(unsigned int));
	cudaMemcpy(R_h_n, R_d, (I_size) * sizeof(unsigned int), 
			cudaMemcpyDeviceToHost);

	int index_grid_size = ( I_size + block_size - 1) / (block_size * 1);
	dim3 index_dimGrid( index_grid_size );

	unsigned int *I_d;
	cudaMalloc((void **)&I_d, (I_size)*sizeof(unsigned int));

	//{{{ index
	start();
	gen_index <<<index_dimGrid,dimBlock>>> ( D_d, D_size, I_d, I_size);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "index: %s.\n", cudaGetErrorString( err) );
	stop();
	unsigned long index_time = report();

	unsigned int *I_h = (unsigned int *) malloc(I_size*sizeof(unsigned int));
	cudaMemcpy(I_h, I_d, (I_size) * sizeof(unsigned int),
			cudaMemcpyDeviceToHost);

	//{{{binary_search_n 
	start();
	binary_search_n <<<dimGrid, dimBlock>>> (I_d, I_size, Q_d, Q_size, R_d);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "binary_search_n: %s.\n", cudaGetErrorString( err) );

	stop();
	printf("%lu\t", report());

	start();
	binary_search_n <<<dimGrid, dimBlock>>> (D_d, D_size, Q_d, Q_size, R_d);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "binary_search_n: %s.\n", cudaGetErrorString( err) );

	stop();
	//unsigned long search_noindex_sorted_2_time = report();
	printf("%lu\n", report());
	//}}}


	///}}}
	//{{{ binary_search_gp
	start();
	binary_search_gp <<< dimGrid, dimBlock, I_size * sizeof(unsigned int) >>> (
			D_d, D_size, Q_d, Q_size, R_d, I_d, I_size);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "binary_search_gp: %s.\n", cudaGetErrorString( err) );

	stop();
	printf("%lu\n", report());
	//}}}
/*

	//{{{ sort Q
	start();
	nvRadixSort::RadixSort sort_Q_d(Q_size, true);
	sort_Q_d.sort((unsigned int*)Q_d, 0, Q_size, 32);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "sort q: %s.\n", cudaGetErrorString( err) );

	stop();
	unsigned long sort_q_time = report();
	//}}}
	
	*/
	
	
	/*
	printf(
			"no_index_search "
			"pre_g_index_search "
			"no_index_sort_search\n"
			"%lu "
			"%lu "
			"%lu\n",
			search_noindex_2_time,
			total_pre_g_index_time,
			search_noindex_sorted_2_time + sort_q_time);
	*/

	/*
	printf(
			"no_index_search\t%lu\n"
			"index_search\t%lu\n"
			"pre_index_search\t%lu\n"
			"no_index_sort\t%lu\n",
			search_noindex_2_time,
			search_index_time,
			total_pre_index_time,
			search_noindex_sorted_2_time + sort_q_time);
	*/
	/*
	printf( "sort_d:%lu\n"
			"sort_q:%lu\n"
			"no_index_search_1:%lu\n"
			"no_index_search_2:%lu\n"
			"index_search:%lu,%f\n"
			"pre_index_search:%lu,%f\ti:%lus:%lu\n"
			"no_index_search_sorted_1:%lu\t%lu\n"
			"no_index_search_sorted_2:%lu\t%lu\n"
			"index_sorted_search:%lu,%f\n"
			"pre_index_sorted_search:%lu,%f\ti:%lus:%lu\n",
			sort_d_time,
			sort_q_time,
			search_noindex_1_time,
			search_noindex_2_time,
			search_index_time, (double)search_noindex_2_time / (double)search_index_time,
			total_pre_index_time, (double)search_noindex_2_time / (double)total_pre_index_time,
			index_time,
			search_pre_index_time,
			search_noindex_sorted_1_time,search_noindex_sorted_1_time + sort_q_time,
			search_noindex_sorted_2_time,search_noindex_sorted_2_time + sort_q_time,

			search_index_sorted_time, (double)search_noindex_2_time / (double)search_index_sorted_time,
			total_pre_index_sorted_time, (double)search_noindex_2_time / (double)total_pre_index_sorted_time,
			index_sorted_time,
			search_pre_index_sorted_time
	);
	*/
	return 0;
}
