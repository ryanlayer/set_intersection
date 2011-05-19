CC=mpicc
AR=ar
CFLAGS= -O3 
#CFLAGS= -Wall -Wno-deprecated -g
LDFLAGS= -lm 

EXECS=num_sim_bsearch_mpi \
	  big_count_bsearch_mpi

BED_FILE_OBJS=../lib/bed.o
SET_FILE_OBJS=../lib/set_intersect.o
TIMER_OBJS=../lib/timer.o
MT_OBJS=../lib/mt.o

OBJS=$(BED_FILE_OBJS) $(SET_FILE_OBJS) $(TIMER_OBJS) $(MT_OBJS)

all: $(EXECS)

#order_mpi: order_mpi.o $(BED_FILE_OBJS) $(SET_FILE_OBJS)
	#$(CC) -o $@  $< $(BED_FILE_OBJS) $(SET_FILE_OBJS)

num_sim_bsearch_mpi: num_sim_bsearch_mpi.o $(OBJS)
	$(CC) -o $@  $< $(OBJS) $(LDFLAGS)

big_count_bsearch_mpi: big_count_bsearch_mpi.o $(OBJS)
	$(CC) -o $@  $< $(OBJS) $(LDFLAGS)

%.o: %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o *.a a.out $(EXECS) $(OBJS)

