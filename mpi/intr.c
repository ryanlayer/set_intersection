#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <dirent.h>
#include <string.h>
#include <limits.h>
#include <mpi.h>
#include "rand_model.h"

struct file_list 
{
	char *name;
	struct file_list *next;
};

struct file_list* make_list(char *dir_name) 
{
	DIR *dir = opendir(dir_name);
	struct dirent *file;

	struct file_list *head, *curr;

	head = NULL;

	while ( file = readdir(dir) ) {
		if ( (strstr(file->d_name, ".bed") != NULL) ) {
			char *n = (char *) malloc( (strlen(file->d_name) + 1) *
					sizeof(char) );
			strcpy(n, file->d_name);
			curr = (struct file_list *) malloc( sizeof(struct file_list) );
			curr->name = n;
			curr->next = head;
			head = curr;
		}
	}
	closedir(dir);
	
	return head;
}

int main(int argc, char** argv)
{
	MPI_Status status;
	int rank, size, buf, max_char = 500;
	char pair[max_char];

	MPI_Init (&argc, &argv);    /* starts MPI */
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);  /* get current process id */
	MPI_Comm_size (MPI_COMM_WORLD, &size);  /* get number of processes */


	if (argc < 4) {
		printf("usage:\t cintr <source dir> <target dir> <universe> <iters>\n");
		return 1;
	}

	char *s_dir = argv[1], *t_dir = argv[2];
	char *cmd = argv[3];

	if (rank == 0) {

		struct file_list *s_head = make_list(s_dir), 
						 *t_head = make_list(t_dir), 
						 *s_curr, *t_curr;

		s_curr = s_head;
		while (s_curr != NULL) {
			t_curr = t_head;
			while (t_curr != NULL) {
				sprintf(pair,"%s/%s,%s/%s\0", 
						s_dir, s_curr->name, 
						t_dir, t_curr->name);
				MPI_Recv(&buf, 1, MPI_INTEGER, MPI_ANY_SOURCE, 0,
						MPI_COMM_WORLD, &status);
				MPI_Send(pair, strlen(pair) + 1, MPI_CHARACTER, status.MPI_SOURCE,
						0, MPI_COMM_WORLD);
				t_curr = t_curr->next;
			}
			s_curr = s_curr->next;
		}

		pair[0] = '\0';

		int left = size - 1;
		while (left > 0) {
			MPI_Recv(&buf, 1, MPI_INTEGER, MPI_ANY_SOURCE, 0,
					MPI_COMM_WORLD, &status);
			MPI_Send(pair, strlen(pair) + 1, MPI_CHARACTER, status.MPI_SOURCE,
					0, MPI_COMM_WORLD);
			left--;
		}
	} else {
		pair[1] = '\0';
		char full_cmd[500];
		while ( 1 ) {
			MPI_Send(&buf, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD);
			MPI_Recv(&pair, max_char, MPI_CHARACTER, MPI_ANY_SOURCE, 0,
				MPI_COMM_WORLD, &status);
			if ( strlen(pair) == 0 )
				break;

			char *file_1 = strtok(pair, ",");
			char *file_2 = strtok(NULL, ",");

			//sprintf(full_cmd,"%s %s %s", cmd, file_1, file_2);
			//system(full_cmd);
			double *test = test_intersection(argv[1], argv[2], argv[3],
					atoi(argv[4]));
			printf("obs:%f mean:%f p:%f\n", test[0], test[1], test[2]);

			//printf("%d\t%s\n", rank, full_cmd);

		}
	}

	MPI_Finalize();
	return 0;
}
