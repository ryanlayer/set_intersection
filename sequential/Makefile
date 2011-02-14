CC=gcc
AR=ar
CFLAGS= -Wall -Wno-deprecated -O3 
#CFLAGS= -Wall -Wno-deprecated -g
LDFLAGS= -lgsl -lgslcblas -lm 

EXECS=order_seq bsearch_seq get_rank_pairs

BED_FILE_OBJS=../lib/bed.o
SET_FILE_OBJS=../lib/set_intersect.o

OBJS=$(BED_FILE_OBJS) $(SET_FILE_OBJS)

all: $(EXECS)

bsearch_seq: bsearch_seq.o $(BED_FILE_OBJS) $(SET_FILE_OBJS)
	$(CC) -o $@  $< $(BED_FILE_OBJS) $(SET_FILE_OBJS)

order_seq: order_seq.o $(BED_FILE_OBJS) $(SET_FILE_OBJS)
	$(CC) -o $@  $< $(BED_FILE_OBJS) $(SET_FILE_OBJS)

get_rank_pairs: get_rank_pairs.o $(BED_FILE_OBJS) $(SET_FILE_OBJS)
	$(CC) -o $@  $< $(BED_FILE_OBJS) $(SET_FILE_OBJS)

%.o: %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f $(EXECS) *.a *.o $(LIBS) $(OBJS)
