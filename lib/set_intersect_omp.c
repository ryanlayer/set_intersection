#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include "../lib/set_intersect.h"


//{{{ int count_intersections_bsearch_omp( struct interval *A_r,
int count_intersections_bsearch_omp( struct interval_triple *A, 
									 struct interval_triple *A_end, 
									 int A_size,
									 struct interval_triple *B, 
									 int B_size,
									 int p)
{
	int i;
	int *c = (int *) calloc(p, sizeof(int) );

	omp_set_num_threads(p);

	#pragma omp parallel for
	for (i = 0; i < B_size; i++) {

		int t_id = omp_get_thread_num();

		int a = interval_triple_bsearch_end(A_end, A_size, B[i].start);
		int b = interval_triple_bsearch_start(A, A_size, B[i].end);

		int num_cant_before = a; 

		while ( A[b].start == B[i].end )
			++b;

		int num_cant_after = A_size - b;

		int num_left = A_size - num_cant_before - num_cant_after;

		c[t_id] += num_left;
	}

	int t_c = 0;
	for (i = 0; i < p; i++)
		t_c += c[i];

	return t_c;
}
// }}}

