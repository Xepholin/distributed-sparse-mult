CC=mpicc
CFLAGS=-g -Wall -Wextra -march=native

OFLAGS=-O3

LFLAGS=-lm -fopenmp -lmpi

FILES=dgemv_mpi.c

NPROCESS=2

all: smult

run: smult
	mpirun -n $(NPROCESS) ./smult data/bcsstk03.mtx 10

smult: $(FILES)
	$(CC) $(CFLAGS) $(OFLAGS) $(FILES) -o $@ $(LFLAGS)

clean:
	@rm -Rf smult

