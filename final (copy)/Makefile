SHELL = /bin/sh
OBJS = parallel.o paralleloptimized.o
CFLAG = -Wall -g -O3
CC = gcc
INCLUDE =
LIBS = -lm -fopenmp

parallel.o: parallel.cpp paralleloptimized.h
	$(CC) -c parallel.cpp -fopenmp

paralleloptimized.o: paralleloptimized.cpp paralleloptimized.h
		$(CC) -c paralleloptimized.cpp -fopenmp -O3

exec:
	${OBJS}	${CC} ${CFLAGS} ${INCLUDES} -o $@ ${OBJS} ${LIBS}
clean:
	-rm -f *.o core *.core
.cpp.o:
	${CC} ${CFLAGS} ${INCLUDES} -c $<
