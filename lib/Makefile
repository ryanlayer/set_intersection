CC=gcc
CFLAGS= -Wall -Wno-deprecated -O3 
#CFLAGS= -Wall -Wno-deprecated -g

%.o: %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f $(EXECS) *.a *.o $(LIBS) $(OBJS)
