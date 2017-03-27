SHELL = /bin/sh
OBJS = main.o paralleloptimized.o
CFLAG = -Wall -g -O3
CC = g++
LIBS = -lm -fopenmp
#-msse3
all:exec

main.o : main.cpp paralleloptimized.h
	$(CC) -g -c main.cpp -fopenmp

paralleloptimized.o : paralleloptimized.cpp paralleloptimized.h
	$(CC) -g -c paralleloptimized.cpp -fopenmp -O3


exec: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LIBS)

.cpp.o: $(CC) $(CFLAGS)  -c $<

clean:
	-rm -f *.o core *.core
