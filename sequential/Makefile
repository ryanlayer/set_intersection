CC=gcc
AR=ar
#CFLAGS= -Wall -g -O3 
CFLAGS= -Wall -g
#CFLAGS= -Wall -Wno-deprecated -g
LDFLAGS= -lgsl -lgslcblas -lm 

EXECS=order_seq bsearch_seq get_rank_pairs get_maped_intervals \
	  scan_seq count_brute_force_seq num_sim_bsearch_seq \
	  enumerate_bsearch_seq \
	  count_bsearch_seq \
	  count_scan_seq \
	  count_sweep_seq \
	  num_sim_scan_seq \
	  s_e_list \
	  nc_list 

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

num_sim_bsearch_seq: num_sim_bsearch_seq.o $(OBJS) 
	$(CC) -o $@  $< $(OBJS)

count_bsearch_seq: count_bsearch_seq.o $(OBJS) 
	$(CC) -o $@  $< $(OBJS)

count_sweep_seq: count_sweep_seq.o $(OBJS) 
	$(CC) -o $@  $< $(OBJS)

count_scan_seq: count_scan_seq.o $(OBJS) 
	$(CC) -o $@  $< $(OBJS)

enumerate_bsearch_seq: enumerate_bsearch_seq.o $(OBJS) 
	$(CC) -o $@  $< $(OBJS)

count_brute_force_seq: count_brute_force_seq.o $(OBJS)
	$(CC) -o $@  $< $(OBJS)

s_e_list: s_e_list.o $(BED_FILE_OBJS) $(SET_FILE_OBJS)
	$(CC) -o $@  $< $(BED_FILE_OBJS) $(SET_FILE_OBJS)

nc_list: nc_list.o $(BED_FILE_OBJS) $(SET_FILE_OBJS)
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
