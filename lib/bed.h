#ifndef __BED_INTERSECT_H__
#define __BED_INTERSECT_H__
#include <stdio.h>

struct interval {
	int start, end;
};

struct interval_node {
	int start, end, offset;
	struct interval_node *next;
};

struct chr_list {
	char *name;
	int size;
	struct interval_node *head;
};

struct bed_line {
	int chr, start, end, offset;
};


struct interval_node *new_interval_node(int start, int end);

struct chr_list new_chr_list(char *name);

int compare_chr_lists (const void *a, const void *b);

void chr_list_insert_interval(struct chr_list *list, int start, int end);

void chr_list_insert_interval_node(struct chr_list *list, 
		struct interval_node *new_interval_node);

void parse_bed_file(FILE *bed_file, struct chr_list chroms[], int chrom_num);

int compare_interval_by_start(const void *a, const void *b);

int compare_interval_node_by_start(const void *a, const void *b);

int compare_interval_node_by_end(const void *a, const void *b);

void free_chr_list(struct chr_list *list, int chrom_num);

int parse_bed_file_by_line(char *bed_file_name, struct bed_line **line);

int chr_list_from_bed_file(struct chr_list **list, char **chrom_names,
		int chrom_num, char *bed_file_name);

int chr_array_from_list(struct chr_list *list, struct bed_line **array, 
		int chrom_num);

int trim(struct chr_list *universe, struct chr_list *interval_set, 
		 int chrom_num);

int chr_name_to_int(char *name);
#endif
