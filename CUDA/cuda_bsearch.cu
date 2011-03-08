#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cutil.h>
#include <sys/time.h>
#include <time.h>


#include "bed.h"
//#include "rand_model.h"
//#include "order.h"
#include "set_intersect.h"
#include "radixsort.h"
#include "gpu.hpp"
#include "random.hpp"

#include "order_kernel.cu"

#include "reduction.h"

//extern void set_ranks(int *val_d, int size);

int main(int argc, char *argv[]) {
	//struct timeval t0_start, t0_end, t1_start, t1_end, t2_start, t2_end;
	//gettimeofday(&t0_start,0);
	//gettimeofday(&t1_start,0);

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

	int i;

	int A_size, B_size, U_size;

	struct bed_line *U_array, *A_array, *B_array;

	U_size = chr_array_from_list(U, &U_array, chrom_num);
	A_size = chr_array_from_list(A, &A_array, chrom_num);
	B_size = chr_array_from_list(B, &B_array, chrom_num);

	/*
	unsigned int *AB_key_h = (unsigned int *) malloc( 
			(A_size + B_size) * sizeof(unsigned int));
	unsigned int *AB_val_h = (unsigned int *) malloc( 
			(A_size + B_size) * sizeof(unsigned int));
	*/
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
	int j,k = 0;
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

	k=0;
	for (i = 0; i < B_size; i++) {
		int start = -1, offset = -1;
		for (j = 0; j < U_size; j++) {
			if ( ( U_array[j].chr == B_array[i].chr) &&
				 ( U_array[j].start <= B_array[i].end) &&
				 ( U_array[j].end >= B_array[i].start) ) {
				start = U_array[j].start;
				offset = U_array[j].offset;
				break;
			}
		}
		B_key_h[k] = B_array[i].start - start + offset;
		B_val_h[k] = B_array[i].end - B_array[i].start;
		++k;
	}

	unsigned int *A_key_d, *A_val_d;
	unsigned int *B_key_d, *B_val_d;
	//unsigned int *AB_key_d, *AB_val_d;
	//cudaMalloc((void **)&A_key_d, (2*A_size)*sizeof(int));
	//cudaMalloc((void **)&A_val_d, (2*A_size)*sizeof(int));
	//cudaMalloc((void **)&B_key_d, (B_size)*sizeof(int));
	//cudaMalloc((void **)&B_val_d, (B_size)*sizeof(int));
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


	nvRadixSort::RadixSort radixsortA(A_size, false);
	nvRadixSort::RadixSort radixsortB(B_size, false);

	radixsortA.sort((unsigned int*)A_key_d, (unsigned int*)A_val_d, 
			A_size, 32);
	

	radixsortB.sort((unsigned int*)B_key_d, (unsigned int*)B_val_d, 
			B_size, 32);

	cudaThreadSynchronize();
	cudaError_t err;
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "Rand A: %s.\n", cudaGetErrorString( err) );


	/*
	cudaMemcpy(A_key_h, A_key_d, A_size * sizeof(unsigned int), 
			cudaMemcpyDeviceToHost);
	cudaMemcpy(A_val_h, A_val_d, A_size * sizeof(unsigned int),
			cudaMemcpyDeviceToHost);

	cudaMemcpy(B_key_h, B_key_d, B_size * sizeof(unsigned int), 
			cudaMemcpyDeviceToHost);
	cudaMemcpy(B_val_h, B_val_d, B_size * sizeof(unsigned int),
			cudaMemcpyDeviceToHost);

	for (i = 0; (i < A_size) && (i < B_size); i++)
			printf("%d,%d\t%d,%d\n", 
				A_key_h[i], A_val_h[i],
				B_key_h[i], B_val_h[i]);
	*/

	int blocksize = 256;
	dim3 dimBlock(blocksize);
	// I want a thread per element in A
	dim3 dimGridA( ceil(float(A_size)/float(dimBlock.x)));

	int *R_di, *R_h, *R_do;
	cudaMalloc((void **)&R_di, (A_size)*sizeof(int));
	cudaMalloc((void **)&R_do, (A_size)*sizeof(int));
	R_h = (int *) malloc( (B_size) * sizeof(int));

	// *_key_d holds the start position, and *_val_d holds the length,
	// the end position is *_key_d + *_val_d
	intersection_b_search <<<dimGridA, dimBlock>>> (
			A_key_d, A_val_d, A_size,
			B_key_d, B_val_d, B_size,
			R_di);

	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(err != cudaSuccess)
		fprintf(stderr, "Rand A: %s.\n", cudaGetErrorString( err) );

	//reduce<int>(A_size, A_size, 

	cudaMemcpy(R_h, R_di, A_size * sizeof(int), cudaMemcpyDeviceToHost);

	for (i = 0; i < A_size; i++)
		printf("%d\n", R_h[i]);


	// set the ranks for the two sets at the same time we will store the size
	// of each interval
	/*
	int *A_len_d, *B_len_d;
	cudaMalloc((void **)&A_len_d, (A_size)*sizeof(int));
	cudaMalloc((void **)&B_len_d, (B_size)*sizeof(int));
	int blocksize = 256;

	dim3 dimBlock(blocksize);

	dim3 dimGridA( ceil(float(2*A_size)/float(dimBlock.x)));
	set_ranks_lens <<<dimGridA, dimBlock>>> (A_val_d, A_key_d, A_len_d, A_size);

	dim3 dimGridB( ceil(float(2*B_size)/float(dimBlock.x)));
	set_ranks_lens <<<dimGridB, dimBlock>>> (B_val_d, B_key_d, B_len_d, B_size);

	// move B over
	cudaMemcpy(AB_key_d + (2*A_size), B_key_d, (2*B_size) * sizeof(int),
			cudaMemcpyDeviceToDevice);
	cudaMemcpy(AB_val_d + (2*A_size), B_val_d, (2*B_size) * sizeof(int),
			cudaMemcpyDeviceToDevice);

	// sort them together
	nvRadixSort::RadixSort radixsortAB(2*A_size + 2*B_size, false);
	radixsortAB.sort((unsigned int*)AB_key_d, (unsigned int*)AB_val_d, 
			2*A_size + 2*B_size, 32);

	// move the results back to the host
	cudaMemcpy(AB_key_h, AB_key_d, (2*A_size + 2*B_size) * sizeof(int),
			cudaMemcpyDeviceToHost);
	cudaMemcpy(AB_val_h, AB_val_d, (2*A_size + 2*B_size) * sizeof(int),
			cudaMemcpyDeviceToHost);
	
	// find the intersecting ranks
	int *pairs_h = (int *) malloc( 2 * (A_size + B_size) * sizeof(int));
	int num_pairs = 0;
	int rankA = -1, rankB = -1;
	bool inB = false, inA = false;
	for (i = 0; i < (2*A_size + 2*B_size); i++) {

		if ( AB_val_h[i] & 1) { //B
			rankB = AB_val_h[i] >> 2;
			inB = !(AB_val_h[i] & 2);
		} else  {//A
			rankA = AB_val_h[i] >> 2;
			inA = !(AB_val_h[i] & 2);
		}

		if (inA && inB) {
			pairs_h[2*num_pairs] = rankA;
			pairs_h[2*num_pairs + 1] = rankB;
			++num_pairs;
		}
	}

	// Don't need these anymore
	cudaFree(AB_key_d);
	cudaFree(B_key_d);
	cudaFree(AB_val_d);
	cudaFree(B_val_d);
	free(AB_val_h);

	if (num_pairs == 0) {
		fprintf(stderr, "No intersections\n");
		return 1;
	}

	//Move rank_pairs to device, we will use these and the len vectors for the
	//simulations
	int *pairs_d;
	cudaMalloc((void **)&pairs_d, (2 * num_pairs) * sizeof(int));
	cudaMemcpy(pairs_d, pairs_h, (2 * num_pairs) * sizeof(int),
			cudaMemcpyHostToDevice);

	gettimeofday(&t1_end,0);
	fprintf(stderr, "setup:%ld\t", 
		(t1_end.tv_sec - t1_start.tv_sec)*1000000 + 
		t1_end.tv_usec - t1_start.tv_usec);


	srand(time(NULL));	

	RNG_rand48 A_r(rand());
	RNG_rand48 B_r(rand());
	dim3 dimGridAR( ceil(float(A_size)/float(dimBlock.x)));
	dim3 dimGridBR( ceil(float(B_size)/float(dimBlock.x)));
	dim3 dimGridT( ceil(float(num_pairs)/float(dimBlock.x)));
	nvRadixSort::RadixSort radixsortAR(A_size, true);
	nvRadixSort::RadixSort radixsortBR(B_size, true);

	int *R_d;
	cudaMalloc((void **)&R_d, (num_pairs)*sizeof(int));
	cudaMemset(R_d, 0, num_pairs);


	// Monte Carlo sims
	gettimeofday(&t1_start,0);

	cudaError_t err;
	for (i = 0; i < reps; i++) {
		//gettimeofday(&t2_start,0); //start
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


		int *A_r_d = (int *)A_r.get_random_numbers();
		int *B_r_d = (int *)B_r.get_random_numbers();

		//cudaThreadSynchronize();  //end
		//gettimeofday(&t2_end,0);
		//fprintf(stderr, "rand:%ld\t", 
			//(t2_end.tv_sec - t2_start.tv_sec)*1000000 + 
			//t2_end.tv_usec - t2_start.tv_usec);


		//gettimeofday(&t2_start,0); //start

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

		//cudaThreadSynchronize();
		//gettimeofday(&t2_end,0); //end
		//fprintf(stderr, "norm:%ld\t", 
			//(t2_end.tv_sec - t2_start.tv_sec)*1000000 + 
			//t2_end.tv_usec - t2_start.tv_usec);


		//gettimeofday(&t2_start,0); //start

		radixsortAR.sort((unsigned int*)A_r_d, 0, A_size, 32);

		cudaThreadSynchronize();
		err = cudaGetLastError();
		if(err != cudaSuccess)
			fprintf(stderr, "Sort A: %s.\n", cudaGetErrorString( err) );

		radixsortBR.sort((unsigned int*)B_r_d, 0, B_size, 32);

		//cudaThreadSynchronize();
		//gettimeofday(&t2_end,0); //end
		//fprintf(stderr, "sort:%ld\t", 
			//(t2_end.tv_sec - t2_start.tv_sec)*1000000 + 
			//t2_end.tv_usec - t2_start.tv_usec);

		cudaThreadSynchronize();
		err = cudaGetLastError();
		if(err != cudaSuccess)
			fprintf(stderr, "Sort B: %s.\n", cudaGetErrorString( err) );


		//gettimeofday(&t2_start,0);  //start

		test_pairs <<<dimGridT, dimBlock>>> (A_r_d, B_r_d, A_len_d, B_len_d, 
				pairs_d, R_d, num_pairs);

		cudaThreadSynchronize();
		err = cudaGetLastError();
		if(err != cudaSuccess)
			fprintf(stderr, "Test: %d %d %d %s.\n", dimGridT.x, dimGridT.y,
					dimGridT.z,
					cudaGetErrorString( err) );

		//cudaThreadSynchronize();
		//gettimeofday(&t2_end,0);  //end
		//fprintf(stderr, "test:%ld\t", 
			//(t2_end.tv_sec - t2_start.tv_sec)*1000000 + 
			//t2_end.tv_usec - t2_start.tv_usec);

		err = cudaGetLastError();
		if(err != cudaSuccess)
			fprintf(stderr, "Cuda error: %s.\n", cudaGetErrorString( err) );
	}

	*/
	return 0;
}
