#include <stdio.h>
#include <stdlib.h>
#include "set_intersect_cuda.h"

//{{{void parallel_sum( int *R_d,
/**
 * @param R_d Address of element array on device
 * @param block_size Number of threads per block
 * @param Rd_size Number of elemens in R_d
 * @param n Number of elemens each thread handles
 */
void parallel_sum( unsigned int *R_d,
				   int block_size,
				   int Rd_size,
				   int n)
{
	unsigned int left = Rd_size;
	//int n = 1024;
	while (left > 1) {

		//dim3 dimGridR( left / (block_size * n) + 1);
		//dim3 dimGridR( left / blocksize  + 1);
		int grid_size = ( left + block_size*n - 1) / (block_size * n);
		dim3 dimGridR( grid_size);

		dim3 dimBlockR( block_size );
		size_t sm_size = dimBlockR.x * sizeof(int); 

		my_reduce <<<dimGridR, dimBlockR, sm_size>>> (R_d, left, n);

		cudaThreadSynchronize();
		cudaError_t err;
		err = cudaGetLastError();
		if(err != cudaSuccess)
			fprintf(stderr, "My Reduce: %s.\n", cudaGetErrorString( err) );

		left = dimGridR.x;
	}
}
//}}}

//{{{__global__ void my_reduce( int *gdata,
__global__
void my_reduce( unsigned int *gdata,
				unsigned int size,
				unsigned int n )
{
	extern __shared__ int sdata[];


	/* v1 load:  need N threads
	unsigned int tid = threadIdx.x;
	unsigned int i = blockIdx.x *  blockDim.x  + tid;

	if (i < size)
		sdata[tid] = gdata[i];
	else
		sdata[tid] = 0;
	__syncthreads();
	*/

	/* v2 load:  need N/2 threads 
	unsigned int tid = threadIdx.x;
	unsigned int i = blockIdx.x * (2 * blockDim.x) + threadIdx.x;
	if (i < size)
		sdata[tid] = gdata[i];
	else
		sdata[tid] = 0;

	if (i + blockDim.x < size)
		sdata[tid] += gdata[i + blockDim.x];
	else
		sdata[tid] += 0;

	__syncthreads();
	*/

	/* v3 load: need N/n threads */
	unsigned int tid = threadIdx.x;
	unsigned int i = blockIdx.x * ( 2 * blockDim.x ) + tid;
	unsigned int grid_size = blockDim.x * ( 2 * gridDim.x);

	sdata[tid] = 0;

	while ( i < (n * grid_size) ) {
		if (i < size)
			sdata[tid] += gdata[i];

		if ( (i + blockDim.x) < size)
			sdata[tid] += gdata[i + blockDim.x];
		i += grid_size;
	}
	__syncthreads();


	/* v1 calc
	unsigned int s;
	for (s = 1; s < blockDim.x; s*=2) {
		if (tid % (2*s) == 0)
			sdata[tid] += sdata[tid + s];

		__syncthreads();
	}
	*/
	/* v2 calc
	unsigned int s;
	for (s = 1; s < blockDim.x; s*=2) {
		int index = 2 * s * tid;
		if (index < blockDim.x)
			sdata[index] += sdata[index + s];

		__syncthreads();
	}
	*/

	/* v3 calc */
	unsigned int s;
	for (s = blockDim.x / 2; s > 0; s >>= 1) {
		if (tid < s)
			sdata[tid] += sdata[tid + s];

		__syncthreads();
	}


	/* v5 calc
	if (blockDim.x >= 512) {
		if (tid < 256)
			sdata[tid] += sdata[tid + 256];
		__syncthreads();
	}
	if (blockDim.x >= 256) {
		if (tid < 128)
			sdata[tid] += sdata[tid + 128];
		__syncthreads();
	}
	if (blockDim.x >= 128) {
		if (tid < 64)
			sdata[tid] += sdata[tid + 64];
		__syncthreads();
	}

	if (tid < 32) {
		if (blockDim.x >= 64)
			sdata[tid] += sdata[tid + 32];
		if (blockDim.x >= 32)
			sdata[tid] += sdata[tid + 16];
		if (blockDim.x >= 16)
			sdata[tid] += sdata[tid + 8];
		if (blockDim.x >= 8)
			sdata[tid] += sdata[tid + 4];
		if (blockDim.x >= 4)
			sdata[tid] += sdata[tid + 2];
		if (blockDim.x >= 4)
			sdata[tid] += sdata[tid + 2];
		if (blockDim.x >= 2)
			sdata[tid] += sdata[tid + 1];
	} */

	if (tid == 0)
		gdata[blockIdx.x] = sdata[0];
}
//}}}

//{{{ __global__ void binary_search_n( unsigned int *db,
__global__
void binary_search_n( unsigned int *db,
					 int size_db, 
					 unsigned int *q,
					 int size_q, 
					 unsigned int *R )
				     
{
	unsigned int id = (blockIdx.x * blockDim.x) + threadIdx.x;
	if ( id < size_q )
		R[id] = binary_search(db, size_db, q[id] );
}
//}}}

//{{{ __global__ void gen_index( unsigned int *db,
__global__
void gen_index( unsigned int *db,
			    int size_db, 
				unsigned int *I,
				int size_I)
				     
{
	//extern __shared__ unsigned int I[];

	unsigned int id = (blockIdx.x * blockDim.x) + threadIdx.x;
	int i;
	if ( id < size_I) {
		i = ((id + 1)*size_db - (size_I - id + 1))/size_I;
		I[id] = db[i];
	}
}
//}}}

//{{{ __global__ void binary_search_p( unsigned int *db,
__global__
void binary_search_p( unsigned int *db,
					 int size_db, 
					 unsigned int *q,
					 int size_q, 
					 unsigned int *R,
					 unsigned int *I,
					 int size_I)
				     
{
	extern __shared__ unsigned int L[];

	unsigned int id = (blockIdx.x * blockDim.x) + threadIdx.x;
	int c, round = 0;

	while ( ( (blockDim.x * round) + threadIdx.x ) < size_I) {
		c = (blockDim.x*round) + threadIdx.x;
		L[c] = I[c];
		++round;
	}
	__syncthreads();

	if (id < size_q) {
		int key = q[id];
		int b = binary_search(L, size_I, key);

		int new_hi = ( (b+1)*size_db - (size_I - (b+2))) / size_I;
		int new_lo = ( (b  )*size_db - (size_I - (b+1))) / size_I;

		if (b == 0)
			new_lo = -1;
		else if (b == size_I) {
			new_hi = size_db;
			//new_lo = ( (b-1)*size_db + (I_size - (b+1))*lo ) / I_size;
			new_lo = ( (b-1)*size_db - (size_I - (b)) ) / size_I;
		}

		R[id] =  bound_binary_search(db, size_db, key, new_lo, new_hi);
	}
}
//}}}

//{{{ __global__ void binary_search_p( unsigned int *db,
__global__
void binary_search_gp( unsigned int *db,
					 int size_db, 
					 unsigned int *q,
					 int size_q, 
					 unsigned int *R,
					 unsigned int *I,
					 int size_I)
				     
{
	unsigned int id = (blockIdx.x * blockDim.x) + threadIdx.x;

	if (id < size_q) {
		int key = q[id];
		int b = binary_search(I, size_I, key);

		int new_hi = ( (b+1)*size_db - (size_I - (b+2))) / size_I;
		int new_lo = ( (b  )*size_db - (size_I - (b+1))) / size_I;

		if (b == 0)
			new_lo = -1;
		else if (b == size_I) {
			new_hi = size_db;
			new_lo = ( (b-1)*size_db - (size_I - (b)) ) / size_I;
		}

		R[id] =  bound_binary_search(db, size_db, key, new_lo, new_hi);
	}
}
//}}}

//{{{ __global__ void binary_search_i( unsigned int *db,
__global__
void binary_search_i( unsigned int *db,
					 int size_db, 
					 unsigned int *q,
					 int size_q, 
					 unsigned int *R,
					 int size_I)
				     
{
	extern __shared__ unsigned int I[];

	unsigned int id = (blockIdx.x * blockDim.x) + threadIdx.x;
	int i, c, round = 0;

	while ( ( (blockDim.x * round) + threadIdx.x ) < size_I) {
		c = (blockDim.x*round) + threadIdx.x;
		i = ((c + 1)*size_db - (size_I - c + 1))/size_I;
		I[c] = db[i];
		++round;
		//if ( blockIdx.x == 0 )
			//R[c] = I[c];
	}
	__syncthreads();

	if (id < size_q) {
		int key = q[id];
		int b = binary_search(I, size_I, key);

		int new_hi = ( (b+1)*size_db - (size_I - (b+2))) / size_I;
		int new_lo = ( (b  )*size_db - (size_I - (b+1))) / size_I;

		if (b == 0)
			new_lo = -1;
		else if (b == size_I) {
			new_hi = size_db;
			//new_lo = ( (b-1)*size_db + (I_size - (b+1))*lo ) / I_size;
			new_lo = ( (b-1)*size_db - (size_I - (b)) ) / size_I;
		}

		R[id] =  bound_binary_search(db, size_db, key, new_lo, new_hi);
		//R[id] =  id;
	}
}
//}}}

//{{{ __device__ int bound_binary_search( unsigned int *db,
__device__
int bound_binary_search( unsigned int *db,
				   int size_db, 
				   unsigned int s,
				   int lo,
				   int hi) 
{
	int mid;
	while ( hi - lo > 1) {
		mid = (hi + lo) / 2;

		if ( db[mid] < s )
			lo = mid;
		else
			hi = mid;
	}
	return hi;
}
//}}}

//{{{ __device__ int binary_search( unsigned int *db, int size_db, unsigned int
__device__
int binary_search( unsigned int *db,
				   int size_db, 
				   unsigned int s) 
{
	int lo = -1, hi = size_db, mid;
	while ( hi - lo > 1) {
		mid = (hi + lo) / 2;

		if ( db[mid] < s )
			lo = mid;
		else
			hi = mid;
	}
	return hi;
}
//}}}

//{{{__device__ int __min(int a, int b) {
__device__
int __min(int a, int b) 
{
	if (a < b)
		return a;
	else
		return b;
}
//}}}

//{{{ __global__ void intersection_b_search_sm_2 ( unsigned int *A_start,
__global__
void intersection_b_search_sm_2 ( unsigned int *A_start,
								  unsigned int *A_len,
								  int A_size,
								  unsigned int *B_start,
								  unsigned int *B_len,
								  int B_size,
								  unsigned int *R,
								  int N )
{
	extern __shared__ unsigned int S[];
	//__shared__ int db_min;
	//__shared__ int db_max;
	
	// Move N elements into shared memory
	int start = blockDim.x * N;
	int num_moved = 0;
	while ( ( threadIdx.x + (num_moved * blockDim.x) <=  N ) &&
			( start + threadIdx.x + (num_moved * blockDim.x) < B_size) ) {
		S[ threadIdx.x + (num_moved * blockDim.x) ] =
				B_start[start + threadIdx.x + (num_moved * blockDim.x)];
		num_moved++;
	}
	__syncthreads();

}
//}}}

//{{{ __global__ void intersection_b_search_sm ( unsigned int *A_start,
/*
 *   We want each thread to take an element in A
 *   Need to load the portion of B that contains all elements searched by A
 *   It is possible that the portion of B need will not fit into A.  In that
 *   case we must split A in half, and re-run on two pieces.
 */
__global__
void intersection_b_search_sm ( unsigned int *A_start,
							 unsigned int *A_len,
							 int A_size,
							 unsigned int *B_start,
							 unsigned int *B_len,
							 int B_size,
							 unsigned int *R,
							 int n )
{
	extern __shared__ unsigned int S_start[];
	__shared__ int db_min;
	__shared__ int db_max;


	unsigned int id = (blockIdx.x * blockDim.x) + threadIdx.x;
	R[id] = 0;

	if (threadIdx.x == 0) 
		db_min = binary_search(B_start, B_size, A_start[id]);

	if (threadIdx.x == 32)
		db_max = binary_search(B_start, B_size, 
				A_start[  __min(blockIdx.x * blockDim.x + blockDim.x - 1, 
								A_size - 1) ] +
				A_len[  __min(blockIdx.x * blockDim.x + blockDim.x - 1, 
								A_size - 1) ] );
				
	__syncthreads();

	//exapnd the region by one in each direction
	int S_size = db_max - db_min + 1 + 2;
	
	if (db_min > 0)
		db_min--;

	/* 
	 * the number of elements from db_d that need to be moved into SM is equal
	 * to db_max - db_min.  If this is smaller than blockDim.x, then threads
	 * will have to move more than one element over 
	 */

	int num_moved = 0;
	while ( threadIdx.x + (num_moved * blockDim.x) <=  S_size ) {
		S_start[ threadIdx.x + (num_moved * blockDim.x) ] =
				B_start[db_min + threadIdx.x + (num_moved * blockDim.x)];
		num_moved++;
	}
	__syncthreads();


	unsigned int i = id;
	unsigned int grid_size = blockDim.x * gridDim.x;

	while ( i < (n * grid_size) ) {

		if (i < A_size) {
			//R[id] = 0;

			unsigned int start = A_start[i];
			unsigned int end = start + A_len[i];

			int left = binary_search(S_start, S_size, start);

			int right = binary_search(S_start, S_size, end);

			int range_start, range_end;

			if ( start == S_start[left] )
				range_start = left;
			else if ( (left > 0) &&
					  ( start <= 
						S_start[left - 1] + B_len[db_min + left - 1] ) )
				range_start = left - 1;
			else 
				range_start = left;

			if ( ( right < S_size ) &&  
				 ( end == S_start[right] ) ) 
				range_end = right;
			else
				range_end = right - 1;

			R[i] = range_end - range_start + 1;
		}
		i += grid_size;
	}
}
//}}}

//{{{ __global__ void intersection_b_search ( unsigned int *A_start,
__global__
void intersection_b_search ( unsigned int *A_start,
							 unsigned int *A_len,
							 int A_size,
							 unsigned int *B_start,
							 unsigned int *B_len,
							 int B_size,
							 unsigned int *R,
							 int n )
{
	unsigned int id = (blockIdx.x * blockDim.x) + threadIdx.x;

	unsigned int i = id;
	unsigned int grid_size = blockDim.x * gridDim.x;

	//R[i] = blockIdx.x;
	while ( i < (n * grid_size) ) {

		if (i < A_size) {
			//R[id] = 0;

			unsigned int start = A_start[i];
			unsigned int end = start + A_len[i];

			int left = binary_search(B_start, B_size, start);

			int right = binary_search(B_start, B_size, end);

			int range_start, range_end;

			if ( A_start[i] == B_start[left] )
				range_start = left;
			else if ( (left > 0) &&
					  ( A_start[i] <= B_start[left - 1] + B_len[left - 1] ) )
				range_start = left - 1;
			else 
				range_start = left;

			if ( ( right < B_size ) &&  
				 ( A_start[i] + A_len[i] == B_start[right] ) ) 
				range_end = right;
			else
				range_end = right - 1;

			R[i] = range_end - range_start + 1;
		}
		i += grid_size;
	}
}
//}}}

//{{{ __global__ void enumerate_b_search_gm ( unsigned int *A_start,
__global__
void enumerate_b_search_gm ( unsigned int *A_start,
							 unsigned int *A_len,
							 int A_size,
							 unsigned int *B_start,
							 unsigned int *B_len,
							 int B_size,
							 unsigned int *P1,
							 unsigned int *P2,
							 int n )
{
	unsigned int id = (blockIdx.x * blockDim.x) + threadIdx.x;

	unsigned int i = id;
	unsigned int grid_size = blockDim.x * gridDim.x;

	//R[i] = blockIdx.x;
	while ( i < (n * grid_size) ) {

		if (i < A_size) {
			//R[id] = 0;

			unsigned int start = A_start[i];
			unsigned int end = start + A_len[i];

			int left = binary_search(B_start, B_size, start);

			int right = binary_search(B_start, B_size, end);

			int range_start, range_end;

			if ( A_start[i] == B_start[left] )
				range_start = left;
			else if ( (left > 0) &&
					  ( A_start[i] <= B_start[left - 1] + B_len[left - 1] ) )
				range_start = left - 1;
			else 
				range_start = left;

			if ( ( right < B_size ) &&  
				 ( A_start[i] + A_len[i] == B_start[right] ) ) 
				range_end = right;
			else
				range_end = right - 1;

			//R[i] = range_end - range_start + 1;
			int j;
			for (j = range_start; j <= range_end; j++) {
				//  This is at worst, the i + j pair
				P1[i + j] = i;
				P2[i + j] = j;
			}

		}
		i += grid_size;
	}
}
//}}}

//{{{ __global__ void set_ranks_lens(int *vald, int *keyd, int *lend, int size)
__global__
void set_ranks_lens( int *vald,
					 int *keyd,
					 int *lend,
					 int size)
{
	int id = (blockIdx.x * blockDim.x) + threadIdx.x;

	if (id < (size*2)) 
		vald[id] = ( (id / 2) << 2 ) + vald[id];

	if (id < size)
		lend[id] = keyd[id*2 + 1] - keyd[id*2];
}
//}}}

//{{{ __global__ void normalize_rand(unsigned int *setd, unsigned int max, int
__global__
void normalize_rand( unsigned int *setd,
					 unsigned int max,
					 int size)
{
	int id = (blockIdx.x * blockDim.x) + threadIdx.x;

	if (id < (size)) 
		setd[id] = setd[id] % max;
}
//}}}

//{{{ __global__ void test_pairs (int *A, int *B, int *A_len, int *B_len, int
__global__
void test_pairs ( int *A,
				  int *B,
				  int *A_len,
				  int *B_len,
				  int *pairs_d,
				  int *R,
				  int size)
{
	int id = (blockIdx.x * blockDim.x) + threadIdx.x;

	if (id < size) {
		int A_rank = pairs_d[id*2];
		int B_rank = pairs_d[id*2 + 1];

		int A_start = A[A_rank];
		int A_end = A[A_rank] + A_len[A_rank];
		int B_start = B[B_rank];
		int B_end = B[B_rank] + B_len[B_rank];

		R[id] += (A_start <= B_end) && (A_end >= B_start);
	}
}
//}}}

//{{{ void intersection_brute_force ( unsigned int *A_start,
__global__
void intersection_brute_force ( unsigned int *A_start,
							 unsigned int *A_len,
							 int A_size,
							 unsigned int *B_start,
							 unsigned int *B_len,
							 int n,
							 unsigned int *R,
							 int offset)
{
	unsigned int id = (blockIdx.x * blockDim.x) + threadIdx.x;
	int i, c = 0;


	//unsigned int i = id;
	//unsigned int grid_size = blockDim.x * gridDim.x;

	if (id < A_size) {
		unsigned int A_s = A_start[id];
		unsigned int A_e = A_start[id] + A_len[id];

		for (i = offset; i < offset + n; i++) {
			unsigned int B_s = B_start[i];
			unsigned int B_e = B_start[i] + B_len[i];

			c += (A_s <= B_e) && (A_e >= B_s);

		}
		//R[id] = c;
		R[id] += c;
	}
}
//}}}

//{{{ __global__ void count_bsearch_cuda (	unsigned int *A_start,
/*
 * @param A_start list of start positions to query, does not need to be sorted
 * @param A_len list of lengths that correspond to A_start
 * @param A_size size of A_start and A_len
 * @param B_start list of sorted start positions to be queried
 * @param B_end list of sorted end positions to be queired 
 * @param B_size size of B_start and B_end
 * @param R number of intersections for each interval in A
 * @param n number of intervals per thread
 */
__global__
void count_bsearch_cuda (	unsigned int *A_start,
							unsigned int *A_len,
							int A_size,
							unsigned int *B_start,
							unsigned int *B_end,
							int B_size,
							unsigned int *R,
							int n)
{
	unsigned int id = (blockIdx.x * blockDim.x) + threadIdx.x;

	unsigned int i = id;
	unsigned int grid_size = blockDim.x * gridDim.x;

	//R[i] = blockIdx.x;
	while ( i < (n * grid_size) ) {

		if (i < A_size) {
			//R[id] = 0;

			unsigned int start = A_start[i];
			unsigned int end = start + A_len[i];

			int cant_before = binary_search(B_end, B_size, start);
			int cant_after = binary_search(B_start, B_size, end);

			while ( end == B_start[cant_after] )
				++cant_after;

			cant_after = A_size - cant_after;	

			R[i] = A_size - cant_before - cant_after;
		}
		i += grid_size;
	}
}
//}}}

//{{{ __global__ void big_count_bsearch_cuda (	unsigned int *A_start,
/*
 * @param A_start list of start positions to query, does not need to be sorted
 * @param A_len list of lengths that correspond to A_start
 * @param A_size size of A_start and A_len
 * @param B_start list of sorted start positions to be queried
 * @param B_end list of sorted end positions to be queired 
 * @param B_size size of B_start and B_end
 * @param R number of intersections for each interval in A
 * @param n number of intervals per thread
 */
__global__
void big_count_bsearch_cuda (	unsigned int *A_start,
							unsigned int *A_len,
							int A_size,
							unsigned int *B_start,
							unsigned int *B_end,
							int B_size,
							unsigned int *R,
							int n)
{
	unsigned int id = (blockIdx.x * blockDim.x) + threadIdx.x;

	unsigned int i = id;
	unsigned int grid_size = blockDim.x * gridDim.x;

	//R[i] = blockIdx.x;
	while ( i < (n * grid_size) ) {

		if (i < A_size) {
			//R[id] = 0;

			unsigned int start = A_start[i];
			unsigned int end = start + A_len[i];

			int cant_before = binary_search(B_end, B_size, start);
			int cant_after = binary_search(B_start, B_size, end);

			while ( end == B_start[cant_after] )
				++cant_after;

			cant_after = A_size - cant_after;	

			R[i] = R[i] + A_size - cant_before - cant_after;
		}
		i += grid_size;
	}
}
//}}}

//{{{ __global__ void normalize_rand(unsigned int *setd, unsigned int max, int
__global__
void set_end( unsigned int *start,
			  unsigned int *end,
			  unsigned int *len,
			  int size)
{
	int id = (blockIdx.x * blockDim.x) + threadIdx.x;

	if (id < (size)) 
		end[id] = start[id] + len[id];
}
//}}}
