all: blastest eigs_test

OBJECTS = ../src/matrices/blasglue.o ../src/matrices/sqmatrix.o ../src/matrices/eigensolver.o   ../src/util/debug_malloc.o  ../src/matrices/evectmatrix.o   ../src/util/mpi_utils.o  ../src/matrices/matrices.o

CC = cc
CFLAGS = -g -I../src/matrices -I../src/util

LIBS = -lessl -lxlf -lxlf90 -lblas -llapack -lm

blastest: blastest.o $(OBJECTS)
	$(CC) $(CFLAGS) blastest.o $(OBJECTS) $(LIBS) -o $@

eigs_test: eigs_test.o $(OBJECTS)
	$(CC) $(CFLAGS) eigs_test.o $(OBJECTS) $(LIBS) -o $@

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f blastest eigs_test core a.out blastest.o eigs_test.o $(OBJECTS)
