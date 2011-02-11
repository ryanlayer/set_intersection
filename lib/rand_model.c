#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <malloc.h>
#include <limits.h>
#include <stdlib.h>
#include <search.h>
#include <sys/time.h>
#include <time.h>
#include "bed.h"
#include "rand_model.h"

int compare_interval_offset (const void *a, const void *b) 
{   
	struct interval *a_i = (struct interval *)a;
	struct interval *b_i = (struct interval *)b;

	return a_i->offset - b_i->offset;
}

/**
 * Take a set of intervals in a chr_list, and project it into a continous array
 * using a mapping defined by the universe
 *
 * @param chrom_num the number of chromosomes
 * @param j the start loction within intervals to begin adding
 * @param universe the genomic range being considered and mappings to a
 * continous array
 * @param interval_set the intervals to add
 * @param intervals the array intervals are mapped to
 * @param sample_num the id for the sample being pushed
 * @return the number of intervals pushed, this may not equal the total number
 * of intervals if some in the set are not within the universe
 */
int push_intervals(int chrom_num, int *j,  struct chr_list *universe,
		struct chr_list *interval_set, struct interval *intervals,
		int sample_num)
{
	int i, c = 0;
	for (i = 0; i < chrom_num; i++) {
		struct interval_node *curr_chrm;

		struct interval_node *curr = interval_set[i].head;
		while (curr != NULL) {
			curr_chrm = universe[i].head;
			int offset = -1, start = -1;
			while (curr_chrm != NULL) {
				/*
				printf( "-- %d %d %d %d %d\n", curr_chrm->start, curr_chrm->end,
						curr->start >= curr_chrm->start, 
						curr->end <= curr_chrm->end,
							(curr->start >= curr_chrm->start) &&
							(curr->end <= curr_chrm->end)	
						);
				*/
				/*
				if (	(curr->start >= curr_chrm->start) &&
						(curr->end <= curr_chrm->end) ) {
				*/
				if (	(curr->start >= curr_chrm->start) &&
						(curr->end <= curr_chrm->end) ) {
					offset = curr_chrm->offset;
					start = curr_chrm->start;

					break;
				}

				curr_chrm = curr_chrm->next;
			}

			if (curr_chrm != NULL) {
				intervals[*j].sample = sample_num;
				intervals[*j].type = 's';
				intervals[*j].offset = curr->start - start + offset;
				if (sample_num == 0)
					intervals[*j].incr = 1;
				else 
					intervals[*j].incr = 0;
				*j = *j + 1;

				intervals[*j].sample = sample_num;
				intervals[*j].type = 'e';
				intervals[*j].offset = curr->end - start + offset;
				if (sample_num == 0)
					intervals[*j].incr = -1;
				else 
					intervals[*j].incr = 0;

				*j = *j + 1;


				++c;

			} else {
				fprintf(stderr, "** %s %d %d\n", universe[i].name, curr->start,
						curr->end);
			}
			curr = curr->next;
		}

	}

	return c;
}

int get_intersections(struct interval *intervals, int total_size) 
{

	qsort(intervals, 2 * total_size, sizeof(struct interval), 
			compare_interval_offset);

	int i, j = 0;

	int num_curr_intervals = 0, last_set = -1;
	for ( i = 0; i < total_size * 2; i++) {


		if (intervals[i].type == 's') 
			++num_curr_intervals;
		else
			--num_curr_intervals;

		if (num_curr_intervals > 1) 
			++j;

		/*
		printf("%d\t", intervals[i].offset);
		if (intervals[i].type == 's') 
			++num_curr_intervals;
		else
			--num_curr_intervals;

		if (num_curr_intervals > 1) {
			printf("!");
			++j;
		}

		printf("\n");
		*/
	}

	return j;
}

double* test_intersection(char *universe_file_name, char *source_file_name, 
						 char *target_file_name, int iters)
{

	struct timeval t_start,t_end;

	int chrom_num = 24;

	/***********************REPLACE WITH INPUT FILE************************/	
	char *chrom_names[] = {
		"chr1",  "chr2",  "chr3", "chr4",  "chr5",  "chr6",  "chr7", "chr8",
		"chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
		"chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX",  "chrY"
	};
	/**********************************************************************/	

	struct chr_list universe[chrom_num], source[chrom_num], target[chrom_num];

	// initialize chrom lists
	int i;
	for (i = 0; i < chrom_num; i++) {
		universe[i] = new_chr_list(chrom_names[i]);
		source[i] = new_chr_list(chrom_names[i]);
		target[i] = new_chr_list(chrom_names[i]);
	}

	// chr_lists need to be sorted before used
	qsort(universe, chrom_num, sizeof(struct chr_list), compare_chr_lists);
	qsort(source, chrom_num, sizeof(struct chr_list), compare_chr_lists);
	qsort(target, chrom_num, sizeof(struct chr_list), compare_chr_lists);

	FILE *universe_file = fopen(universe_file_name, "r");
	FILE *source_file = fopen(source_file_name, "r");
	FILE *target_file = fopen(target_file_name, "r");

	if (	(universe_file == NULL) || (source_file == NULL) || 
			(target_file == NULL) ) {
		fprintf(stderr, "%s\n", strerror(errno));
		return 0;
	}

	parse_bed_file(universe_file, universe, chrom_num);
	parse_bed_file(source_file, source, chrom_num);
	parse_bed_file(target_file, target, chrom_num);

	fclose(universe_file);
	fclose(source_file);
	fclose(target_file);

	trim(universe, source, chrom_num);
	trim(universe, target, chrom_num);

	// Calculate the offsets for each iterval
	int c = 0;
	int max = 0;
	for (i = 0; i < chrom_num; i++) {
		struct interval_node *curr = universe[i].head;
		while (curr != NULL) {
			curr->offset = c;

			int end = c + curr->end - curr->start;
			if (end > max)
				max = end;

			c += curr->end - curr->start;
			curr = curr->next;
		}
	}

	int total_size = 0, target_size = 0, source_size = 0;

	/* 
	 * Get the total number of intervals in the source and target sets, we also
	 * store the number of intervals individually to be used later for
	 * randomization
	 */
	for (i = 0; i < chrom_num; i++) {
		total_size += source[i].size;
		total_size += target[i].size;
		target_size += target[i].size;
		source_size += source[i].size;
	}

	/* 
	 * get an array of just the target intervals, each random permutation will
	 * consist of these intervals and a set of randomly generated intervals
	 */
	struct interval *target_intervals = (struct interval *) malloc( 
			2 * target_size * sizeof(struct interval) );
	i = 0;

	int pushed_targets = push_intervals(chrom_num, &i,  universe, target,
			target_intervals, 0);

	/* 
	 * we need will permute the intervals in target, to manage this we will
	 * create and array of the interval sizes in source
	 */
	int *rand_sizes = (int *) malloc( source_size * sizeof(int) );
	int j = 0;
	for (i = 0; i < chrom_num; i++) {
		struct interval_node *curr = source[i].head;
		while (curr != NULL) {
			rand_sizes[j++] = curr->end - curr->start;			
			curr = curr->next;
		}
	}

	/* 
	 * set up an array with target and source to find the observed number of 
	 * intersections
	 */
	struct interval *intervals = (struct interval *) malloc( 
			2 * total_size * sizeof(struct interval) );

	for (i = 0; i < 2 * target_size; i++)
		intervals[i] = target_intervals[i];

	int pushed_sources = push_intervals(chrom_num, &i,  universe, source,
			intervals, 1);

	gettimeofday(&t_start,0);
	int obs = get_intersections(intervals, total_size);
	gettimeofday(&t_end,0);

	double p_of_source = (double)obs / (double)source_size;
	double p_of_target = (double)obs / (double)target_size;

	int r = 0;
	int sum = 0;

	// do some random stuff
	struct interval rand_start, rand_end;
	rand_start.type = 's';
	rand_end.type = 'e';
	rand_start.sample = 1;
	rand_end.sample = 1;

	for (i = 0; i < iters; i++) {
		for (j = 0; j < 2 * target_size; j++)
			intervals[j] = target_intervals[j];

		for (j = 0; j < source_size; j++) {
			rand_start.offset = rand() % (max - rand_sizes[j]);
			rand_end.offset = rand_start.offset + rand_sizes[j];
			intervals[ 2 * target_size + 2 * j ] = rand_start;
			intervals[ 2 * target_size + 2 * j + 1 ] = rand_end;
		}

		int t = get_intersections(intervals, total_size);

		sum += t;

		if ( t >= obs )
			++r;
	}

	double p = ( (double)(r + 1) ) /  ( (double)(iters + 1) );
	double mean = ( (double)(sum) ) /  ( (double)(iters) );

	double *ret = (double *) malloc(5 * sizeof(double));

	ret[0] = obs;
	ret[1] = mean;
	ret[2] = p;
	ret[3] = p_of_source;
	ret[4] = p_of_target;
	
	free_chr_list(universe, chrom_num);
	free_chr_list(source, chrom_num);
	free_chr_list(target, chrom_num);
	free(target_intervals);
	free(rand_sizes);
	free(intervals); 

	return ret;
}
