#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

#include "bed.h"

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

int main(int argc, char *argv[]) {
	struct timeval t0_start, t0_end, t1_start, t1_end, t2_start, t2_end;
	gettimeofday(&t0_start,0);
	gettimeofday(&t1_start,0);


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
	int reps = atoi(argv[4]);

	if((chr_list_from_bed_file(&U_list, chrom_names, chrom_num, U_file) == 1) ||
	   (chr_list_from_bed_file(&A_list, chrom_names, chrom_num, A_file) == 1) ||
	   (chr_list_from_bed_file(&B_list, chrom_names, chrom_num, B_file) == 1) ){

		fprintf(stderr, "Error parsing bed files.\n");
		return 1;
	}

	trim(U_list, A_list, chrom_num);
	trim(U_list, B_list, chrom_num);

	// Calculate the offsets for each iterval
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

	int A_size, B_size, U_size;

	struct bed_line *U_array, *A_array, *B_array;

	U_size = chr_array_from_list(U_list, &U_array, chrom_num);
	A_size = chr_array_from_list(A_list, &A_array, chrom_num);
	B_size = chr_array_from_list(B_list, &B_array, chrom_num);


	// make one large array to hold these
	struct triple *AB = (struct triple *)
			malloc((2*A_size + 2*B_size)*sizeof(struct triple));
	struct triple *A = AB;
	struct triple *B = AB + 2*A_size;

	int j,k = 0;
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
		A[k].key = A_array[i].start - start + offset;
		A[k].type = 0;// start
		A[k].sample = 0;// A
		++k;
		A[k].key = A_array[i].end - start + offset;
		A[k].type = 1;// end
		A[k].sample = 0;// A
		++k;
	}

	k=0;
	for (i = 0; i < B_size; i++) {
		int start = -1, offset = -1;
		for (j = 0; j < U_size; j++) {
			if ( ( U_array[j].chr == B_array[i].chr) &&
				 ( U_array[j].start <= B_array[i].end) &&
				 ( U_array[j].end >= B_array[i].start) ) {
				start = U_array[j].start;
				offset = U_array[j].offset;
				break;
			}
		}
		B[k].key = B_array[i].start - start + offset;
		B[k].type = 0;//start
		B[k].sample = 1;// B
		++k;
		B[k].key = B_array[i].end - start + offset;
		B[k].type = 1; //end
		B[k].sample = 1;// B
		++k;
	}

	// sort A and B so they can be ranked
	qsort(A, 2*A_size, sizeof(struct triple), compare_triple_lists);
	qsort(B, 2*B_size, sizeof(struct triple), compare_triple_lists);

	int *A_len = (int *) malloc(A_size * sizeof(int));
	int *B_len = (int *) malloc(B_size * sizeof(int));

	for (i = 0; i < 2*A_size; i++)
		A[i].rank = i/2;

	for (i = 0; i < A_size; i++)
		A_len[i] = A[i*2 + 1].key - A[i*2].key;

	for (i = 0; i < 2*B_size; i++) 
		B[i].rank = i/2;

	for (i = 0; i < B_size; i++)
		B_len[i] = B[i*2 + 1].key - B[i*2].key;

	qsort(AB, 2*A_size + 2*B_size, sizeof(struct triple), compare_triple_lists);

	// find the intersecting ranks
	int *pairs = (int *) malloc( 2 * (A_size + B_size) * sizeof(int));
	int num_pairs = 0;
	int rankA = -1, rankB = -1;
	int inB = 0, inA = 0;
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
			//fprintf(stderr, "%d\t%d\t%d\n", rankA, rankB, num_pairs);
		}
	}

	//sprintf(stderr, "O\t%d\n", num_pairs);

	gettimeofday(&t1_end,0);
/*
	fprintf(stderr, "setup:%ld\t", 
		(t1_end.tv_sec - t1_start.tv_sec)*1000000 + 
		t1_end.tv_usec - t1_start.tv_usec);
*/


	int *A_r = (int *) malloc (A_size * sizeof(int));
	int *B_r = (int *) malloc (B_size * sizeof(int));
	int *R = (int *) calloc(num_pairs, sizeof(int));
	int A_start, A_end, B_start, B_end;

	gettimeofday(&t1_start,0);

	int r = 0;

	srand((unsigned)time(NULL));

	for (j = 0; j < reps; j++) {
		// Fill and sort A and B
		for (i = 0; i < A_size; i++)
			A_r[i] = rand() % max;
		qsort(A_r, A_size, sizeof(int), compare_ints);

		for (i = 0; i < B_size; i++)
			B_r[i] = rand() % max;
		qsort(B_r, B_size, sizeof(int), compare_ints);


		/**** DEBUG
		for (i = 0; i < A_size; i++)
			fprintf(stderr, "A\t%d\t%d\n", A_r[i], A_len[i]);
		for (i = 0; i < B_size; i++)
			fprintf(stderr, "B\t%d\t%d\n", B_r[i], B_len[i]);

		for (i = 0; i < num_pairs; i++) {
			fprintf(stderr, "%d,%d\t", pairs[i*2], pairs[i*2 + 1]);	
		}
		fprintf(stderr, "\n");
		****/

		// See if the observed ranks overlap
		int x = 0;
		for (i = 0; i < num_pairs; i++) {
			rankA = pairs[i*2];	
			rankB = pairs[i*2 + 1];	
			A_start = A_r[rankA];
			A_end = A_r[rankA] + A_len[rankA];
			B_start = B_r[rankB];
			//B_end = B_r[rankB] + A_len[rankB];
			B_end = B_r[rankB] + B_len[rankB];

			/**** DEBUG
			fprintf(stderr,"%d,%d,%d\t%d,%d,%d\t%d\n",
					A_start, A_end, A_len[rankA], B_start, B_end, B_len[rankB],
					(A_start <= B_end) && (A_end >= B_start));
			****/

			R[i] += (A_start <= B_end) && (A_end >= B_start);
			x += (A_start <= B_end) && (A_end >= B_start);
		}

		fprintf(stderr,"\t%d\n", x);

		if (x >= num_pairs)
			r++;
	}

	gettimeofday(&t1_end,0);
	/*
	fprintf(stderr, "sim:%ld\t", 
		(t1_end.tv_sec - t1_start.tv_sec)*1000000 + 
		t1_end.tv_usec - t1_start.tv_usec);
	*/
	gettimeofday(&t0_end,0);
	/*
	fprintf(stderr, "total:%ld\n", 
		(t0_end.tv_sec - t0_start.tv_sec)*1000000 + 
		t0_end.tv_usec - t0_start.tv_usec);
	*/

	/* Print the significance of each intersection
	for (i = 0; i < num_pairs; i++)
		printf("R\t%d\n", R[i]);
	*/

	printf("o:%d\tr:%d\tn:%d\n", num_pairs,r,reps);

	return 0;
}
