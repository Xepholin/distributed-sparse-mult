CC=mpicc
CFLAGS=-g -Wall -Wextra -march=native

OFLAGS=-O3

LFLAGS=-lmpi

FILES=main.c

NPROCESS=4

all: smult

run: smult
	mpirun -n $(NPROCESS) ./smult data/bcsstk03.mtx 100

smult: $(FILES)
	$(CC) $(CFLAGS) $(OFLAGS) $(FILES) -o $@ $(LFLAGS)

clean:
	@rm -Rf smult

