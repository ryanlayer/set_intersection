#include "../lib/bed.h"

struct triple {
	unsigned int key, sample, type, rank;
};

int compare_triple_lists (const void *a, const void *b) {
	struct triple *a_i = (struct triple *)a;
	struct triple *b_i = (struct triple *)b;
	return a_i->key - b_i->key;
}

int compare_ints (const void *a, const void *b) {
	unsigned int *a_i = (unsigned int *)a;
	unsigned int *b_i = (unsigned int *)b;
	return *a_i - *b_i;
}


//{{{ void set_start_len( struct bed_line *U_array,
void set_start_len( struct bed_line *U_array,
					int U_size,
					struct bed_line *A_array,
					unsigned int *A_key_h,
					unsigned int *A_val_h,
					int A_size )
{
	int i, j, k = 0;
	for (i = 0; i < A_size; i++) {
		int start = -1, offset = -1;
		for (j = 0; j < U_size; j++) {
			if ( ( U_array[j].chr == A_array[i].chr) &&
				 ( U_array[j].start <= A_array[i].end) &&
				 ( U_array[j].end >= A_array[i].start) ) {

				start = U_array[j].start;
				offset = U_array[j].offset;
				break;
			}
		}
		A_key_h[k] = A_array[i].start - start + offset;
		A_val_h[k] = A_array[i].end -A_array[i].start;
		++k;
	}
}
//}}}

//{{{ int count_intersections_bsearch( struct interval *A_r,
int count_intersections_bsearch( unsigned int *A_start,
								 unsigned int *A_len,
								 int A_size,
								 unsigned int *B_start,
								 unsigned int *B_len,
								 int B_size )
{
	int i, c = 0;
	for (i = 0; i < A_size; i++) {
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

		//printf("i:%d\tl:%d\tr:%d\n", i, left, right);

		/* This is the way to save the intersecting pairs

		// Check to see if the start is in an interval
		int first_hit = 0;
		if ( ( A_start[i] == B_start[left] ) ) {

			++c;

			printf("%d (%u,%u)\t%d (%u,%u) %d\n",
					i, A_start[i], A_start[i] + A_len[i],
					left, B_start[left], B_start[left] + B_len[left],
					c);
			first_hit = 1;

		} else if ( ( left > 0 ) && 
					(A_start[i] <= B_start[left - 1] + B_len[left - 1] ))   {
			++c;

			printf("%d (%u,%u)\t%d (%u,%u) %d\n",
					i, A_start[i], A_start[i] + A_len[i],
					left - 1, B_start[left - 1], 
					B_start[left - 1] + B_len[left - 1],
					c);
		}

		// Check to see if the end is in an interval
		int k;
		for (k = left + first_hit; (k <= right) && (k < B_size); k++) {
			if ( (A_start[i] + A_len[i] >= B_start[k]) ) 
				++c;

				printf("%d (%u,%u)\t%d (%u,%u) %d\n",
						i, A_start[i], A_start[i] + A_len[i],
						k, B_start[k], B_start[k] + B_len[k],
						c);
		}
		*/

		/* v1
		*/
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

		c += range_end - range_start + 1;
	}

	return c;
}
// }}}

//{{{ int count_intersections_scan( int *A, 
/*
 * Scan two sorted lists counting overlaps
 */
int count_intersections_scan( unsigned int *A, 
							  unsigned int *A_len, 
							  int A_size,
							  unsigned int *B, 
							  unsigned int *B_len,
							  int B_size )
{

	int curr_A = 0, curr_B = 0;
	unsigned int A_val = 0, B_val = 0;
	int inA = 0, inB = 0;

	int o = 0;

	//while ( (curr_A < A_size - 1 ) || (curr_B < B_size - 1 ) ) {
	while ( (curr_A < A_size ) || (curr_B < B_size ) ) {

		// The current values depend on if we are in or not in a segment
		if ( inA )
			A_val = A[curr_A] + A_len[curr_A];
		else
			A_val = A[curr_A];

		if ( inB )
			B_val = B[curr_B] + B_len[curr_B];
		else
			B_val = B[curr_B];

		if ( (curr_A < A_size - 1 ) && (curr_B < B_size - 1 ) ) {

			// Move the pointer
			if ( A_val < B_val ) {

				if (inA)
					++curr_A;
				inA = !inA;

			} else if ( A_val > B_val ) {

				if (inB)
					++curr_B;
				inB = !inB;

			} else { // A_val == B_val

				// both are ending (EE)
				if (inA && inB) {
					++curr_A;
					++curr_B;
					inA = !inA;
					inB = !inB;

				// A is ending B is starting (ES)
				} else if (inA && !inB) {
					//++curr_B;
					inB = !inB;
				// A is starting B is ending (SE)
				} else if (!inA && inB) {
					//++curr_A;
					inA = !inA;
				// both are starting (SS)
				} else { 
					inA = !inA;
					inB = !inB;
				}
			}
		} else if (curr_A < A_size - 1 ) {
			if (inA)
				++curr_A;
			inA = !inA;
		} else if (curr_B < B_size - 1 ) {
			if (inB)
				++curr_B;
			inB = !inB;
		}

		if (inA && inB)  {
			printf("%d (%u,%u)\t%d (%u,%u)\n", 
					curr_A, A[curr_A], A[curr_A] + A_len[curr_A],
					curr_B, B[curr_B], B[curr_B] + B_len[curr_B]);
			++o;
		}
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
unsigned int add_offsets( struct chr_list *U_list, 
						  int chrom_num )
{
	int i;
	unsigned int c = 0, max = 0;
	for (i = 0; i < chrom_num; i++) {
		struct interval_node *curr = U_list[i].head;
		while (curr != NULL) {
			curr->offset = c;

			unsigned int end = c + curr->end - curr->start;
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
							 unsigned int *pairs)
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
						  unsigned int *A_r,
						  unsigned int *A_len,
						  unsigned int *B_r,
						  unsigned int *B_len,
						  int num_pairs,
						  int *R )
{
	int x = 0;
	int rankA, rankB;
	unsigned int A_start, A_end, B_start, B_end;
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
		unsigned int start = 0, offset = 0;
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
