#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <dirent.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include <mpi.h>

#include "rand_model.h"

#define XRAND_MAX (RAND_MAX*(RAND_MAX + 2))

unsigned int xrand ()
{
	unsigned int a = rand();
	unsigned int b = rand();
	unsigned int one = 1;
	return a * (RAND_MAX + one) + b;
}

void print_usage (char * prog) 
{
	fprintf(stderr,
		"\nUSAGE\n"
		"%s [options]\n\n"

		"DESCRIPTION\n"
		"  Simple test to have threads generate random numbers and send\n"
		"  them back to thread 0.\n\n"

		"OPTIONS\n"
		"  -h       Print this help message\n"
		"  -s  seed Random seed, each thread will add id (default time())\n"
		, prog
	);
}

int main(int argc, char** argv)
{
	MPI_Status status;
	int rank, size;
	unsigned buf;

	MPI_Init (&argc, &argv);    /* starts MPI */
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);  /* get current process id */
	MPI_Comm_size (MPI_COMM_WORLD, &size);  /* get number of processes */

	/*
	if (argc < 4) {
		printf("usage:\t cintr <source dir> <target dir> <universe> <iters>\n");
		return 1;
	}
	*/

	int c;
	unsigned int seed = time(0);

	while ((c = getopt (argc, argv, "hs:")) != -1)
		switch (c) {
			case 'h':
				print_usage(argv[0]);
				MPI_Finalize();
				return 1;
			case 's':
				seed = atoi(optarg);
				break;
			default:
				//fprintf (stderr, "Unknown option `-%c'.\n", optopt);
				print_usage(argv[0]);
				MPI_Finalize();
				return 1;
		}


	if (rank == 0) {
		int seen = 0;
		while (seen < (size - 1)) {
			MPI_Recv(&buf, 1, MPI_UNSIGNED, MPI_ANY_SOURCE, 0,
					MPI_COMM_WORLD, &status);
			printf("%d\tgot:%u\tfrom:%d\n", rank, buf, status.MPI_SOURCE);
			++seen;
		}
	} else {
		srand(seed + rank);
		buf = xrand();
		MPI_Send(&buf, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD);
	}

	MPI_Finalize();
	return 0;
}
