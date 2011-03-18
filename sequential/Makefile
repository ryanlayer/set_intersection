CC=gcc
AR=ar
CFLAGS= -Wall -g -O3 
#CFLAGS= -Wall -Wno-deprecated -g
LDFLAGS= -lgsl -lgslcblas -lm 

EXECS=order_seq bsearch_seq get_rank_pairs get_maped_intervals \
	  scan_seq brute_force

BED_FILE_OBJS=../lib/bed.o
SET_FILE_OBJS=../lib/set_intersect.o
TIMER_OBJS=../lib/timer.o
MT_OBJS=../lib/mt.o

OBJS=$(BED_FILE_OBJS) $(SET_FILE_OBJS) $(TIMER_OBJS) $(MT_OBJS)

all: $(EXECS)

bsearch_seq: bsearch_seq.o $(BED_FILE_OBJS) $(SET_FILE_OBJS) $(TIMER_OBJS)
	$(CC) -o $@  $< $(BED_FILE_OBJS) $(SET_FILE_OBJS) $(TIMER_OBJS)

scan_seq: scan_seq.o $(BED_FILE_OBJS) $(SET_FILE_OBJS)
	$(CC) -o $@  $< $(BED_FILE_OBJS) $(SET_FILE_OBJS)

num_sim_scan_seq: num_sim_scan_seq.o $(OBJS) 
	$(CC) -o $@  $< $(OBJS)

brute_force: brute_force.o $(BED_FILE_OBJS) $(SET_FILE_OBJS)
	$(CC) -o $@  $< $(BED_FILE_OBJS) $(SET_FILE_OBJS)

order_seq: order_seq.o $(BED_FILE_OBJS) $(SET_FILE_OBJS)
	$(CC) -o $@  $< $(BED_FILE_OBJS) $(SET_FILE_OBJS)

get_rank_pairs: get_rank_pairs.o $(BED_FILE_OBJS) $(SET_FILE_OBJS)
	$(CC) -o $@  $< $(BED_FILE_OBJS) $(SET_FILE_OBJS)

get_maped_intervals: get_maped_intervals.o $(BED_FILE_OBJS) $(SET_FILE_OBJS)
	$(CC) -o $@  $< $(BED_FILE_OBJS) $(SET_FILE_OBJS)

%.o: %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f $(EXECS) *.a *.o $(LIBS) $(OBJS)
