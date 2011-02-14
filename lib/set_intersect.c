#include "../lib/bed.h"

struct triple {
	int key, sample, type, rank;
};

int compare_triple_lists (const void *a, const void *b) {
	struct triple *a_i = (struct triple *)a;
	struct triple *b_i = (struct triple *)b;
	return a_i->key - b_i->key;
}

int compare_ints (const void *a, const void *b) {
	int *a_i = (int *)a;
	int *b_i = (int *)b;
	return *a_i - *b_i;
}

int count_intersections_bsearch( struct interval *A_r,
								 int A_size,
								 struct interval *B_r,
								 int B_size )
{
	int i, c = 0;
	for (i = 0; i < A_size; i++) {
		// Search for the left-most interval in B with the start in A
		int lo = -1, hi = B_size, mid;
		while ( hi - lo > 1) {
			mid = (hi + lo) / 2;

			if ( B_r[mid].start < A_r[i].start ) 
				lo = mid;
			else
				hi = mid;

		}

		int left = hi;

		lo = -1;
		hi = B_size;
		while ( hi - lo > 1) {
			mid = (hi + lo) / 2;

			if ( B_r[mid].start < A_r[i].end ) 
				lo = mid;
			else
				hi = mid;
		}

		int right = hi;

		/* This is the way to save the intersecting pairs
		for (k = left; k <= right; k++) {
			if ( (k > 0) && (A_r[i].start < B_r[k - 1].end))
			   printf("%d\t%d\n", i, k - 1);	
		}
		*/

		c += (right - left) + 
				( (left > 0)  && (A_r[i].start < B_r[left - 1].end) );

	}

	return c;
}


//{{{ int count_intersections_scan( int *A, 
/*
 * Scan two sorted lists counting overlaps
 */
int count_intersections_scan( int *A, 
							  int *A_len, 
							  int A_size,
							  int *B, 
							  int *B_len,
							  int B_size )
{

	int curr_A = 0, curr_B = 0;
	int A_val = 0, B_val = 0;
	int inA = 0, inB = 0;

	int o = 0;

	while ( (curr_A < A_size) && (curr_B < B_size) ) {

		//printf("A:%d,%d\tB:%d,%d\n", curr_A, A_size, curr_B, B_size);

		// The current values depend on if we are in or not in a segment
		if ( inA )
			A_val = A[curr_A] + A_len[curr_A];
		else
			A_val = A[curr_A];

		if ( inB )
			B_val = B[curr_B] + B_len[curr_B];
		else
			B_val = B[curr_B];

		// Move the pointer
		if ( A_val < B_val ) {
			if (inA)
				++curr_A;

			inA = !inA;
		} else {
			if (inB)
				++curr_B;

			inB = !inB;
		}

		if (inA && inB) 
			++o;
	}

	return o;
}
//}}}

//{{{ int add_offsets( struct chr_list *U_list, 
/*
 * The universe defines the space under consideration.  Each interval in the
 * universe is given an offset in the continguous space.  This offset is used
 * to map intervals in the sample to the continguous space. 
 *
 */
int add_offsets( struct chr_list *U_list, 
				  int chrom_num )
{
	int i, c = 0, max = 0;
	for (i = 0; i < chrom_num; i++) {
		struct interval_node *curr = U_list[i].head;
		while (curr != NULL) {
			curr->offset = c;

			int end = c + curr->end - curr->start;
			if (end > max)
				max = end;

			c += curr->end - curr->start;
			curr = curr->next;
		}
	}

	return max;
}
//}}}

//{{{ int count_intersections( struct triple *AB,
/*
 * Requires AB to be sorted
 */
int count_intersections( struct triple *AB,
						 int A_size,
						 int B_size)
{
	int num_pairs = 0;
	int inB = 0, inA = 0;
	int i;
	for (i = 0; i < (2*A_size + 2*B_size); i++) {
		if ( AB[i].sample == 1) // B
			inB = !(AB[i].type);
		else  // A
			inA = !(AB[i].type);

		if (inA && inB) 
			++num_pairs;
	}

	return num_pairs;
}
//}}}

//{{{ int find_intersecting_ranks( struct triple *AB,
/*
 * Requires AB to be sorted
 */
int find_intersecting_ranks( struct triple *AB,
							 int A_size,
							 int B_size,
							 int *pairs)
{

	int num_pairs = 0;
	int rankA = -1, rankB = -1;
	int inB = 0, inA = 0;
	int i;
	for (i = 0; i < (2*A_size + 2*B_size); i++) {
		if ( AB[i].sample == 1) { //B
			rankB = AB[i].rank;
			inB = !(AB[i].type);
		} else  {//A
			rankA = AB[i].rank;
			inA = !(AB[i].type);
		}

		if (inA && inB) {
			pairs[2*num_pairs] = rankA;
			pairs[2*num_pairs + 1] = rankB;
			++num_pairs;
		}
	}

	return num_pairs;
}
//}}}

//{{{ int check_observed_ranks( int *pairs,
int check_observed_ranks( int *pairs,
						  int *A_r,
						  int *A_len,
						  int *B_r,
						  int *B_len,
						  int num_pairs,
						  int *R )
{
	int x = 0;
	int rankA, rankB;
	int A_start, A_end, B_start, B_end;
	int i;
	for (i = 0; i < num_pairs; i++) {
		rankA = pairs[i*2];	
		rankB = pairs[i*2 + 1];	
		A_start = A_r[rankA];
		A_end = A_r[rankA] + A_len[rankA];
		B_start = B_r[rankB];
		B_end = B_r[rankB] + B_len[rankB];

		/**** DEBUG
		fprintf(stderr,"%d,%d,%d\t%d,%d,%d\t%d\n",
				A_start, A_end, A_len[rankA], B_start, B_end, B_len[rankB],
				(A_start <= B_end) && (A_end >= B_start));
		****/

		R[i] += (A_start <= B_end) && (A_end >= B_start);
		x += (A_start <= B_end) && (A_end >= B_start);
	}

	return x;
}
//}}}

//{{{ void map_intervals( struct triple *A, 
/*
 * Each interval becomes a triple: 
 *   key:  offset
 *   sample:  A (0) or B (1)
 *   type:  start (0) or  end (1)
 *   rank: order within
 *
 */
void map_intervals( struct triple *A, 
					struct bed_line *A_array,
					int A_size, 
					struct bed_line *U_array, 
					int U_size,
					int sample)
{
	int i, j, k = 0;
	for (i = 0; i < A_size; i++) {
		int start = -1, offset = -1;
		// find the universe interval the current interval is in
		// so we can caclulate its place in the continous space
		for (j = 0; j < U_size; j++) {
			if ( ( U_array[j].chr == A_array[i].chr) &&
				 ( U_array[j].start <= A_array[i].end) &&
				 ( U_array[j].end >= A_array[i].start) ) {
				start = U_array[j].start;
				offset = U_array[j].offset;
				//printf("s:%d\to:%d\n", start, offset);
				break;
			}
		}
		A[k].key = A_array[i].start - start + offset;
		A[k].type = 0;// start
		A[k].sample = sample; // A(0) or B(1)
		++k;
		A[k].key = A_array[i].end - start + offset;
		A[k].type = 1;// end
		A[k].sample = sample; // A(0) or B(1)
		++k;
	}
}
//}}}
