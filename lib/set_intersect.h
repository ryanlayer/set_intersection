#include "bed.h"

#ifndef __SET_INTERSECT_H__
#define __SET_INTERSECT_H__

struct triple {
	unsigned int key, sample, type, rank;
};

void set_start_len( struct bed_line *U_array,
					int U_size,
					struct bed_line *A_array,
					unsigned int *A_key_h,
					unsigned int *A_val_h,
					int A_size );

int compare_triple_lists (const void *a, const void *b);

int compare_ints (const void *a, const void *b);

int count_intersections_bsearch( unsigned int *A_start,
								 unsigned int *A_len,
								 int A_size,
								 unsigned int *B_start,
								 unsigned int *B_len,
								 int B_size );

int count_intersections_scan( unsigned int *A, 
							  unsigned int *A_len, 
							  int A_size,
							  unsigned int *B, 
							  unsigned int *B_len,
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
						  unsigned int *A_r,
						  unsigned int *A_len,
						  unsigned int *B_r,
						  unsigned int *B_len,
						  int num_pairs,
						  int *R );

void map_intervals( struct triple *A, 
					struct bed_line *A_array,
					int A_size, 
					struct bed_line *U_array, 
					int U_size,
					int sample);

#endif
