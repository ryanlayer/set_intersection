
#ifndef __SET_INTERSECT__OMP_H__
#define __SET_INTERSECT__OMP_H__

int count_intersections_bsearch_omp( unsigned int *A_start,
								 unsigned int *A_len,
								 int A_size,
								 unsigned int *B_start,
								 unsigned int *B_len,
								 int B_size,
								 int p);
#endif
