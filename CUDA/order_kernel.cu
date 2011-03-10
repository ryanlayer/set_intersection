#ifndef __ORDER_KERNEL_H__
#define __ORDER_KERNEL_H__

__global__
void intersection_b_search(
		unsigned int *A_start, unsigned int *A_len, int A_size,
		unsigned int *B_start, unsigned int *B_len, int B_size,
		int *R)
{
	int id = (blockIdx.x * blockDim.x) + threadIdx.x;

	if (id < A_size) {
		//R[id] = 0;

		unsigned int start = A_start[id];
		unsigned int end = start + A_len[id];

		int lo = -1, hi = B_size, mid;
		while ( hi - lo > 1) {
			mid = (hi + lo) / 2;

			if ( B_start[mid] < start )
				lo = mid;
			else
				hi = mid;
		}

		int left = hi;

		lo = -1;
		hi = B_size;
		while ( hi - lo > 1) {
			mid = (hi + lo) / 2;

			if ( B_start[mid] < end )
				lo = mid;
			else
				hi = mid;
		}

		int right = hi;

		int range_start, range_end;

		/* v1 */
		if ( A_start[id] == B_start[left] )
			range_start = left;
		else if ( (left > 0) &&
				  ( A_start[id] <= B_start[left - 1] + B_len[left - 1] ) )
			range_start = left - 1;
		else 
			range_start = left;

		if ( ( right < B_size ) &&  
			 ( A_start[id] + A_len[id] == B_start[right] ) ) 
			range_end = right;
		else
			range_end = right - 1;

		R[id] = range_end - range_start + 1;
		/*
		R[id] = (right - left) + ( (left > 0) && 
				(start < (B_start[left - 1] + B_len[left - 1]) ) );
		*/
	}
}
__global__
void set_ranks_lens(int *vald, int *keyd, int *lend, int size)
{
	int id = (blockIdx.x * blockDim.x) + threadIdx.x;

	if (id < (size*2)) 
		vald[id] = ( (id / 2) << 2 ) + vald[id];

	if (id < size)
		lend[id] = keyd[id*2 + 1] - keyd[id*2];
}

__global__
void normalize_rand(unsigned int *setd, unsigned int max, int size)
{
	int id = (blockIdx.x * blockDim.x) + threadIdx.x;

	if (id < (size)) 
		setd[id] = setd[id] % max;
}

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

__global__
void my_reduce( int *gdata,
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

#endif
