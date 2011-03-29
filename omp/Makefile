CC=gcc
AR=ar
CFLAGS= -Wall -g -O3 
#CFLAGS= -Wall -Wno-deprecated -g
LDFLAGS= -lgsl -lgslcblas -lm 

EXECS= num_sim_bsearch_omp \
	  enumerate_bsearch_omp \
	  count_bsearch_omp

BED_FILE_OBJS=../lib/bed.o
SET_FILE_OBJS=../lib/set_intersect.o
SET_OMP_FILE_OBJS=../lib/set_intersect_omp.o
TIMER_OBJS=../lib/timer.o
MT_OBJS=../lib/mt.o

OBJS=$(BED_FILE_OBJS) $(SET_FILE_OBJS) $(SET_OMP_FILE_OBJS) $(TIMER_OBJS) \
	 $(MT_OBJS)

all: $(EXECS)

num_sim_bsearch_omp: num_sim_bsearch_omp.o $(OBJS) 
	$(CC) -fopenmp -o $@  $< $(OBJS)

count_bsearch_omp: count_bsearch_omp.o $(OBJS) 
	$(CC) -fopenmp -o $@  $< $(OBJS)

enumerate_bsearch_omp: enumerate_bsearch_omp.o $(OBJS) 
	$(CC) -o $@  $< $(OBJS)

../lib/set_intersect_omp.o: ../lib/set_intersect_omp.c
	$(CC) -fopenmp -c $(CFLAGS) $< -o $@

%.o: %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f $(EXECS) *.a *.o $(LIBS) $(OBJS)
