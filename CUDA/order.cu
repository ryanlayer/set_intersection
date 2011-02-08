#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cutil.h>
#include <sys/time.h>
#include <time.h>


#include "bed.h"
#include "rand_model.h"
#include "order.h"

#include "radixsort.h"

#include "gpu.hpp"
#include "random.hpp"

#include "order_kernel.cu"

//#define INTERVAL_START	(1 << 0) /* 0x01 */
//#define INTERVAL_END	(1 << 1) /* 0x02 */
//#define SAMPLE A		(1 << 2) /* 0x04 */
//#define SAMPLE B		(1 << 3) /* 0x08 */



//extern void set_ranks(int *val_d, int size);

int main(int argc, char *argv[]) {
	struct timeval t0_start, t0_end, t1_start, t1_end, t2_start, t2_end;
	gettimeofday(&t0_start,0);
	gettimeofday(&t1_start,0);

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

	trim(U, A, chrom_num);
	trim(U, B, chrom_num);

	// Calculate the offsets for each iterval
	int i, c = 0, max = 0;
	for (i = 0; i < chrom_num; i++) {
		struct interval_node *curr = U[i].head;
		while (curr != NULL) {
			curr->offset = c;

			int end = c + curr->end - curr->start;
			if (end > max)
				max = end;

			c += curr->end - curr->start;
			curr = curr->next;
		}
	}

	int A_size, B_size, U_size;

	struct bed_line *U_array, *A_array, *B_array;

	U_size = chr_array_from_list(U, &U_array, chrom_num);
	A_size = chr_array_from_list(A, &A_array, chrom_num);
	B_size = chr_array_from_list(B, &B_array, chrom_num);


	int *AB_key_h = (int *) malloc( (2*A_size + 2*B_size) * sizeof(int));
	int *AB_val_h = (int *) malloc( (2*A_size + 2*B_size) * sizeof(int));
	int *A_key_h = AB_key_h;
	int *A_val_h = AB_val_h;
	int *B_key_h = AB_key_h + 2*A_size;
	int *B_val_h = AB_val_h + 2*A_size;

	/*
	 * In CUDA we can sort key value pairs, the key can be the offset, but the
	 * the value must represent several values: rank, type, sample.  To
	 * accomplish this we will take the rank, shift it 2 bits, the least
	 * significant bit would then be used to identify the sample (A=0 or B=1)
	 * and the next most significant bit would identify if the particular value
	 * was a start or end (s=0 e=1)
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
		A_val_h[k] = 0;// 0 0  <start, A>
		++k;
		A_key_h[k] = A_array[i].end - start + offset;
		A_val_h[k] = 2;// 1 0 <end, A>
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
		B_val_h[k] = 1;// 0 1  <start, B>
		++k;
		B_key_h[k] = B_array[i].end - start + offset;
		B_val_h[k] = 3;// 1 1  <end, B>
		++k;
	}


	int *A_key_d, *A_val_d;
	int *B_key_d, *B_val_d;
	int *AB_key_d, *AB_val_d;
	//cudaMalloc((void **)&A_key_d, (2*A_size)*sizeof(int));
	//cudaMalloc((void **)&A_val_d, (2*A_size)*sizeof(int));
	cudaMalloc((void **)&B_key_d, (2*B_size)*sizeof(int));
	cudaMalloc((void **)&B_val_d, (2*B_size)*sizeof(int));
	cudaMalloc((void **)&AB_key_d, (2*A_size + 2*B_size)*sizeof(int));
	cudaMalloc((void **)&AB_val_d, (2*A_size + 2*B_size)*sizeof(int));
	A_key_d = AB_key_d;
	A_val_d = AB_val_d;

	cudaMemcpy(A_key_d, A_key_h, (2*A_size) * sizeof(int), 
			cudaMemcpyHostToDevice);
	cudaMemcpy(A_val_d, A_val_h, (2*A_size) * sizeof(int),
			cudaMemcpyHostToDevice);

	cudaMemcpy(B_key_d, B_key_h, (2*B_size) * sizeof(int), 
			cudaMemcpyHostToDevice);
	cudaMemcpy(B_val_d, B_val_h, (2*B_size) * sizeof(int),
			cudaMemcpyHostToDevice);

	nvRadixSort::RadixSort radixsortA(2*A_size, false);
	nvRadixSort::RadixSort radixsortB(2*B_size, false);
 
	radixsortA.sort((unsigned int*)A_key_d, (unsigned int*)A_val_d, 
			2*A_size, 32);
	radixsortB.sort((unsigned int*)B_key_d, (unsigned int*)B_val_d, 
			2*B_size, 32);

	// set the ranks for the two sets at the same time we will store the size
	// of each interval
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

	gettimeofday(&t1_end,0);
	fprintf(stderr, "sim:%ld\t", 
		(t1_end.tv_sec - t1_start.tv_sec)*1000000 + 
		t1_end.tv_usec - t1_start.tv_usec);


	/*
	int *R_h = (int *) malloc( (num_pairs)*sizeof(int) );
	cudaMemcpy(R_h, R_d, num_pairs* sizeof(int), cudaMemcpyDeviceToHost);

	for (i = 0; i < num_pairs; i++)
		printf("%d\n", R_h[i]);
		*/

	/*
	cudaFree(A_len_d);
	cudaFree(B_len_d);
	cudaFree(R_d);
	*/

	gettimeofday(&t0_end,0);
	fprintf(stderr, "total:%ld\n", 
		(t0_end.tv_sec - t0_start.tv_sec)*1000000 + 
		t0_end.tv_usec - t0_start.tv_usec);

	return 0;
}
