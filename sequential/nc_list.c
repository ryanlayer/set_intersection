#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

#include "../lib/bed.h"
#include "../lib/set_intersect.h"

#define MIN(a,b) ((a)>(b)?(b):(a))

struct nc_sublist {
	int start, length;
};

struct nc_list {
	int start, end, sublist;
};

//{{{int contains(struct nc_list *A, struct nc_list *B) {
// See if A contains B
int contains(struct nc_list *A, struct nc_list *B) {
	return ( ( A->start <= B->start ) && ( A->end > B->end));
}
//}}}

//{{{int compare_nc_list_by_start (const void *a, const void *b)
int compare_nc_list_by_start (const void *a, const void *b)
{   
	struct nc_list *a_i = (struct nc_list *)a;
	struct nc_list *b_i = (struct nc_list *)b;

	if (a_i->start < b_i->start)
		return -1;
	else if (a_i->start > b_i->start)
		return 1;
	else if (a_i->end > b_i->end)
		return -1;
	else if (a_i->end < b_i->end)
		return 1;
	else
		return 0;
}
//}}}

//{{{int compare_nc_list_by_sublist_start(const void *a,const void *b)
int compare_nc_list_by_sublist_start(const void *a,const void *b)
{ /* SORT IN SUBLIST ORDER, SECONDARILY BY start */

	struct nc_list *a_i = (struct nc_list *)a;
	struct nc_list *b_i = (struct nc_list *)b;

	if (a_i->sublist < b_i->sublist)
		return -1;
	else if (a_i->sublist > b_i->sublist)
		return 1;
	else if (a_i->start < b_i->start)
		return -1;
	else if (a_i->start > b_i->start)
		return 1;
	else
		return 0;
}
//}}}

//{{{ void map_to_nc_list( struct nc_list *A, 
void map_to_nc_list( struct nc_list *A, 
					struct bed_line *A_array,
					int A_size, 
					struct bed_line *U_array, 
					int U_size,
					int sample)
{
	int i, j;
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
		A[i].start = A_array[i].start - start + offset;
		A[i].end = A_array[i].end - start + offset;
		A[i].sublist = 0;
	}
}
//}}}

int main(int argc, char *argv[]) {
	//struct timeval t0_start, t0_end, t1_start, t1_end, t2_start, t2_end;
	//gettimeofday(&t0_start,0);
	//gettimeofday(&t1_start,0);
	
	if (argc < 4) {
		fprintf(stderr, "usage: order <u> <a> <b> <N>\n");
		return 1;
	}

	int chrom_num = 24;

	/***********************REPLACE WITH INPUT FILE************************/	
	char *chrom_names[] = {
		"chr1",  "chr2",  "chr3", "chr4",  "chr5",  "chr6",  "chr7", "chr8",
		"chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
		"chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX",  "chrY"
	};
	/**********************************************************************/	

	struct chr_list *U_list, *A_list, *B_list;

	char *U_file = argv[1], *A_file = argv[2], *B_file = argv[3];

	if((chr_list_from_bed_file(&U_list, chrom_names, chrom_num, U_file) == 1) ||
	   (chr_list_from_bed_file(&A_list, chrom_names, chrom_num, A_file) == 1) ||
	   (chr_list_from_bed_file(&B_list, chrom_names, chrom_num, B_file) == 1) ){

		fprintf(stderr, "Error parsing bed files.\n");
		return 1;
	}

	int max = add_offsets(U_list, chrom_num);

	if (max == 0) {
		fprintf(stderr, "Max is zero.\n");
		return 1;
	}

	trim(U_list, A_list, chrom_num);
	trim(U_list, B_list, chrom_num);


	int A_size, B_size, U_size;

	struct bed_line *U_array, *A_array, *B_array;

	// Move the universe and intervals from linked lists to arrays
	U_size = chr_array_from_list(U_list, &U_array, chrom_num);
	A_size = chr_array_from_list(A_list, &A_array, chrom_num);
	B_size = chr_array_from_list(B_list, &B_array, chrom_num);


	// make one large array to hold these
	/* 
	 * We need to put both A and B into a single array then sort it
	 *
	 * Each interval becomes a triple: 
	 *   key:  offset
	 *   sample:  A (0) or B (1)
	 *   type:  start (0) or  end (1)
	 *   rank: order within
	 *
	 */
	struct nc_list *A = (struct nc_list *)
			malloc((A_size)*sizeof(struct nc_list));
	struct nc_list *B = (struct nc_list *)
			malloc((A_size)*sizeof(struct nc_list));

	map_to_nc_list(A, A_array, A_size, U_array, U_size, 0 );
	map_to_nc_list(B, B_array, B_size, U_array, U_size, 1 );

	int i;
	// sort A and B so they can be ranked
	qsort(A, A_size, sizeof(struct nc_list), compare_nc_list_by_start);

	for (i = 0; i < A_size; i++)
		printf("i:%d\ts:%d\te:%d\tsl:%d\n", 
				i,
				A[i].start, A[i].end, A[i].sublist);




	// Find the number of sublists
	int n_lists = 1;
	for (i = 1; i < A_size; i++ ) {
		if ( contains(&A[i-1], &A[i]) ) 
			++n_lists;
	}
	printf("n_lists:%d\n", n_lists);


	int p_nlists = n_lists - 1;

	struct nc_sublist *H = (struct nc_sublist *)
			malloc((n_lists + 1)*sizeof(struct nc_sublist));

	A[0].sublist = 0;
	H[0].start = -1;
	H[0].length = 1;
	int parent = 0, i_sublist = 1;
	n_lists = 1;
	i = 1;

	while (i < A_size) {

		if ( ( i_sublist && (A[i].end > A[parent].end) ) || // i not contained
				( (A[i].end == A[parent].end) && // same interval
				  (A[i].start == A[parent].start) ) ) {
			// record the position, within its sublist, the sublist starts
			H[i_sublist].start = H[ A[parent].sublist ].length - 1; 

			// this sublist is over, so move back to the containing sublist
			i_sublist = A[parent].sublist;
			// move parent to the containing sublist parent
			parent = H[ A[parent].sublist ].start;
		} else {

			if( H[i_sublist].length ==0 )
				n_lists++;

			H[i_sublist].length++;
			A[i].sublist = i_sublist;
			parent = i;
			i_sublist = n_lists;
			H[i_sublist].start = parent;
			i++;
		}
	}

	// If we stop in the middle of a sublist
	while (i_sublist > 0) {
		// record parent relative position
		H[ i_sublist].start = H[ A[parent].sublist ].length - 1; 
		i_sublist = A[parent].sublist;
		parent = H[ A[parent].sublist ].start;
	}

	// Set absolute positions
	int total = 0, temp;
	for(i = 0; i < n_lists + 1 ; ++i){
		temp = H[i].length;
		H[i].length = total;
		total += temp;
	}

	/* SUBHEADER.LEN IS NOW START OF THE SUBLIST */

	for(i = 1; i < A_size; i++ ) {
		if( A[i].sublist > A[i-1].sublist )
			H[ A[i].sublist ].start += H[ A[i-1].sublist ].length;
	}

	/* SUBHEADER.START IS NOW ABS POSITION OF PARENT */

	qsort(A, A_size, sizeof(struct nc_list), compare_nc_list_by_sublist_start);


	for (i = 0; i < A_size; i++)
		printf("i:%d\ts:%d\te:%d\tsl:%d\n", 
				i,
				A[i].start, A[i].end, A[i].sublist);

	printf("\n");

	/* AT THIS POINT SUBLISTS ARE GROUPED TOGETHER, READY TO PACK */

	i_sublist = 0;
	H[0].start = 0;
	H[0].length = 0;

	for (i = 0; i < A_size; ++i){
		if (A[i].sublist > i_sublist){
			/* 
			printf("Entering sublist %d (%d,%d)\n", A[i].sublist,
					A[i].start, A[i].end);
			*/
			i_sublist = A[i].sublist;
			parent = H[i_sublist].start;
			/*
			printf("Parent (%d,%d) is at %d, list start is at %d\n",
					A[parent].start, A[parent].end,
					H[i_sublist].start,i); 
			*/
			//A[parent].sublist = i_sublist - 1;
			A[parent].sublist = i_sublist - 1;
			H[i_sublist].length=0;
			H[i_sublist].start = i;
		}

		H[i_sublist].length++;
		A[i].sublist = -1;
	}

	n_lists--;
	H = H + 1;

	for (i = 0; i < A_size; i++)
		printf("i:%d\ts:%d\te:%d\tsl:%d\n", 
				i,
				A[i].start, A[i].end, A[i].sublist);

	printf("\n");
	printf("n_lists:%d\n",n_lists);
	for (i = 0; i < n_lists; i++)
		printf("i:%d\ts:%d\tl:%d\n", 
				i,
				H[i].start, H[i].length);

	/*
	struct pair *A_start = (struct pair *) malloc(A_size*sizeof(struct pair));
	struct pair *A_end = (struct pair *) malloc(A_size*sizeof(struct pair));

	for (i = 0; i < 2*A_size; i++) {
		printf("=r:%u t:%u k:%u s:%u\n", 
				A[i].rank,
				A[i].type,
				A[i].key,
				A[i].sample);
		if (A[i].type == 0) {
			A_start[ i/2 ].rank = A[i].rank;
			A_start[ i/2 ].key = A[i].key;
			printf("-%d %u, %u\n", i, A_start[ i/2 ].rank, A[i].rank);
		} else {
			A_end[ i/2 ].rank = A[i].rank;
			A_end[ i/2 ].key = A[i].key;
			printf("+%d %u, %u\n", i, A_end[ i/2 ].rank, A[i].rank);
		}
	}

	for (i = 0; i < A_size; i++) {
		printf("%u\t%u\n", A_start[i].rank, A_end[i].rank);
	}

	*/
	return 0;
}
