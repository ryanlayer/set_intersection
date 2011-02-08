#ifndef __ORDER_KERNEL_H__
#define __ORDER_KERNEL_H__

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
void normalize_rand(int *setd, int max, int size)
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
#endif
