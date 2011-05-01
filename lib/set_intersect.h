#include "bed.h"

#ifndef __SET_INTERSECT_H__
#define __SET_INTERSECT_H__

struct triple {
	unsigned int key, sample, type, rank;
};

struct pair
{
	unsigned int key, rank;
};

struct interval_triple
{
	unsigned int start, end, rank;
};

struct interval_pair
{
	unsigned int start, end;
};

void set_start_len( struct bed_line *U_array,
					int U_size,
					struct bed_line *A_array,
					unsigned int *A_start,
					unsigned int *A_len,
					int A_size );

void set_start_end( struct bed_line *U_array,
					int U_size,
					struct bed_line *A_array,
					unsigned int *A_start,
					unsigned int *A_end,
					int A_size );


int compare_interval_triples_start_to_end (const void *key, const void *b);

int compare_interval_triples_by_start (const void *a, const void *b);

int compare_interval_triples_by_end (const void *a, const void *b);

int compare_interval_pairs_by_start (const void *a, const void *b);

int compare_pair_lists (const void *a, const void *b);

int compare_triple_lists (const void *a, const void *b);

int compare_uints (const void *a, const void *b);

int count_intersections_bsearch( unsigned int *A_start,
								 unsigned int *A_len,
								 int A_size,
								 unsigned int *B_start,
								 unsigned int *B_len,
								 int B_size );

int count_intersections_brute_force( unsigned int *A, 
							  unsigned int *A_len, 
							  int A_size,
							  unsigned int *B, 
							  unsigned int *B_len,
							  int B_size );

int enumerate_intersections_bsearch( unsigned int *A_start,
								 unsigned int *A_len,
								 int A_size,
								 unsigned int *B_start,
								 unsigned int *B_len,
								 int B_size,
								 unsigned int *pairs);

int count_intersections_scan( unsigned int *A, 
							  unsigned int *A_len, 
							  int A_size,
							  unsigned int *B, 
							  unsigned int *B_len,
							  int B_size );

int count_intersections_sweep( struct triple *AB, 
							   int A_size,
							   int B_size );

unsigned int add_offsets( struct chr_list *U_list, 
				  int chrom_num );

int find_intersecting_ranks( struct triple *AB,
							 int A_size,
							 int B_size,
							 unsigned int *pairs);

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

void map_to_interval_triple( struct interval_triple *A, 
					struct bed_line *A_array,
					int A_size, 
					struct bed_line *U_array, 
					int U_size,
					int sample);

void map_to_interval_pair( struct interval_pair *A, 
					struct bed_line *A_array,
					int A_size, 
					struct bed_line *U_array, 
					int U_size,
					int sample);

#endif
