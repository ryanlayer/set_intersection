CC=gcc
AR=ar
CFLAGS= -Wall -Wno-deprecated -O3 
#CFLAGS= -Wall -Wno-deprecated -g
LDFLAGS= -lgsl -lgslcblas -lm 

EXECS=order_seq

BED_FILE_OBJS=bed.o

OBJS=$(BED_FILE_OBJS) 

all: $(EXECS)

order_seq: order_seq.o $(BED_FILE_OBJS)
	$(CC) -o $@  $< $(BED_FILE_OBJS)

%.o: %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f $(EXECS) *.a *.o $(LIBS) $(OBJS)
