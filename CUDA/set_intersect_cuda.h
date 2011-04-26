#ifndef __SET_INTERSECT_CUDA_H__
#define __SET_INTERSECT_CUDA_H__

void parallel_sum( unsigned int *R_d,
				   int block_size,
				   int Rd_size,
				   int n);

__global__
void my_reduce( unsigned int *gdata,
				unsigned int size,
				unsigned int n );
__device__
int binary_search( unsigned int *db,
				   int size_db, 
				   unsigned int s);
__device__
int __min(int a, int b);

/*
__global__
void intersection_b_search_sm_2 ( unsigned int *A_start,
								  unsigned int *A_len,
								  int A_size,
								  unsigned int *B_start,
								  unsigned int *B_len,
								  int B_size,
								  unsigned int *R,
								  int N );
								  */
__global__
void intersection_b_search_sm ( unsigned int *A_start,
							 unsigned int *A_len,
							 int A_size,
							 unsigned int *B_start,
							 unsigned int *B_len,
							 int B_size,
							 unsigned int *R,
							 int n );
__global__
void intersection_b_search ( unsigned int *A_start,
							 unsigned int *A_len,
							 int A_size,
							 unsigned int *B_start,
							 unsigned int *B_len,
							 int B_size,
							 unsigned int *R,
							 int n );

__global__
void enumerate_b_search_gm ( unsigned int *A_start,
							 unsigned int *A_len,
							 int A_size,
							 unsigned int *B_start,
							 unsigned int *B_len,
							 int B_size,
							 unsigned int *P1,
							 unsigned int *P2,
							 int n );

__global__
void set_ranks_lens( int *vald,
					 int *keyd,
					 int *lend,
					 int size);

__global__
void normalize_rand( unsigned int *setd,
					 unsigned int max,
					 int size);

__global__
void test_pairs ( int *A,
				  int *B,
				  int *A_len,
				  int *B_len,
				  int *pairs_d,
				  int *R,
				  int size);

__global__
void intersection_brute_force ( unsigned int *A_start,
							 unsigned int *A_len,
							 int A_size,
							 unsigned int *B_start,
							 unsigned int *B_len,
							 int B_size,
							 unsigned int *R);
#endif
