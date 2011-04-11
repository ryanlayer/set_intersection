#include <omp.h>
#include <stdio.h>
#include <stdlib.h>


//{{{ int count_intersections_bsearch_omp( struct interval *A_r,
int count_intersections_bsearch_omp( unsigned int *A_start,
								 unsigned int *A_len,
								 int A_size,
								 unsigned int *B_start,
								 unsigned int *B_len,
								 int B_size,
								 int p)
{
	int i;
	int *c = (int *) calloc(p, sizeof(int) );

	omp_set_num_threads(p);


	#pragma omp parallel for
	for (i = 0; i < A_size; i++) {
		//printf("id\t%d\n", omp_get_thread_num());
		// Search for the left-most interval in B with the start in A
		int lo = -1, hi = B_size, mid;
		while ( hi - lo > 1) {
			mid = (hi + lo) / 2;

			if ( B_start[mid] < A_start[i] ) 
				lo = mid;
			else
				hi = mid;

		}

		int left = hi;

		lo = -1;
		hi = B_size;
		while ( hi - lo > 1) {
			mid = (hi + lo) / 2;

			if ( B_start[mid] < (A_start[i] + A_len[i]) ) 
				lo = mid;
			else
				hi = mid;
		}

		int right = hi;

		int range_start, range_end;

		if ( ( A_start[i] == B_start[left] ) ) {
			range_start = left;
		} else if ( ( left > 0 ) &&
					( A_start[i] <= B_start[left - 1] + B_len[left - 1]) ) {
			range_start = left - 1;
		} else {
			range_start = left;
		}


		if ( (right < B_size) &&
			 ( A_start[i] + A_len[i] == B_start[right] ) ) {
			range_end = right;
		} else {
			range_end = right - 1;
		} 

		c[p] += range_end - range_start + 1;
	}

	int t_c = 0;
	for (i = 0; i < p; i++)
		t_c += c[i];

	return t_c;
}
// }}}
