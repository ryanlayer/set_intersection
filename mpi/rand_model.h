#ifndef __RAND_MODEL_H__
#define __RAND_MODEL_H__

#include "bed.h"

struct interval 
{
	char type;
	int offset;
	int sample;
	int incr;
};

int compare_interval_offset (const void *a, const void *b);

int push_intervals(int chrom_num, int *j,  struct chr_list *universe,
		struct chr_list *interval_set, struct interval *intervals,
		int sample_num); 

int get_intersections(struct interval *intervals, int total_size);

double* test_intersection(char *universe_file_name, char *source_file_name, 
						 char *target_file_name, int iters);

#endif
