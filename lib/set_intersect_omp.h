
#ifndef __SET_INTERSECT__OMP_H__
#define __SET_INTERSECT__OMP_H__
int count_intersections_bsearch_omp( struct interval_triple *A, 
									 struct interval_triple *A_end, 
									 int A_size,
									 struct interval_triple *B, 
									 int B_size,
									 int p);
#endif
