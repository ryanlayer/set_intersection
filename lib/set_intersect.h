#include "bed.h"

#ifndef __SET_INTERSECT_H__
#define __SET_INTERSECT_H__

struct triple {
	int key, sample, type, rank;
};

int compare_triple_lists (const void *a, const void *b);

int compare_ints (const void *a, const void *b);

int count_intersecitons_scan( int *A, 
							  int *A_len, 
							  int A_size,
							  int *B, 
							  int *B_len,
							  int B_size );
int add_offsets( struct chr_list *U_list, 
				  int chrom_num );

int count_intersections( struct triple *AB,
						 int A_size,
						 int B_size);

int find_intersecting_ranks( struct triple *AB,
							 int A_size,
							 int B_size,
							 int *pairs);

int check_observed_ranks( int *pairs,
						  int *A_r,
						  int *A_len,
						  int *B_r,
						  int *B_len,
						  int num_pairs,
						  int *R );

void map_intervals( struct triple *A, 
					struct bed_line *A_array,
					int A_size, 
					struct bed_line *U_array, 
					int U_size,
					int sample);

#endif
