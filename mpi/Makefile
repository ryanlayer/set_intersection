CC=mpicc
AR=ar
CFLAGS= -Wall -O3 
#CFLAGS= -Wall -Wno-deprecated -g
LDFLAGS= -lm 

EXECS=order_mpi num_sim_bsearch_mpi

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

%.o: %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o *.a a.out $(EXECS) $(OBJS)

