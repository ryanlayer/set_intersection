#ifndef __ORDER_KERNEL_H__
#define __ORDER_KERNEL_H__

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


//{{{__device__ int binary_search(int *db, int size_db, int s) {
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

//{{{ __global__ void set_ranks_lens(int *vald, int *keyd, int *lend, int size)
__global__
void set_ranks_lens(int *vald, int *keyd, int *lend, int size)
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
void normalize_rand(unsigned int *setd, unsigned int max, int size)
{
	int id = (blockIdx.x * blockDim.x) + threadIdx.x;

	if (id < (size)) 
		setd[id] = setd[id] % max;
}
//}}}

//{{{ __global__ void test_pairs (int *A, int *B, int *A_len, int *B_len, int
__global__
void test_pairs (int *A, int *B, int *A_len, int *B_len, int *pairs_d, int *R, int size) {
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

#endif
