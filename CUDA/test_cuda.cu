#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cutil.h>
#include <sys/time.h>

#include "set_intersect_cuda.h"

int main(int argc, char *argv[]) {

	cudaFree(NULL);

	int block_size = 256;
	int n=1, size=1000;
	dim3 dimBlock(block_size);

	unsigned int *R_d, *R_h;
	cudaMalloc((void **)&R_d, (size)*sizeof(unsigned int));
	R_h = (unsigned int *) malloc ((size)*sizeof(unsigned int));


	int i;
	for (i = 0; i < size; i++)
		R_h[i] = 1;

	cudaMemcpy(R_d, R_h, size* sizeof(unsigned int), cudaMemcpyHostToDevice);

	parallel_sum( R_d, block_size, size, n );

	cudaMemcpy(R_h, R_d, size* sizeof(unsigned int), cudaMemcpyDeviceToHost);

	for (i = 0; i < 5; i++)
		printf("%i\n", R_h[i]);

	cudaFree(R_d);

	return 0;
}
