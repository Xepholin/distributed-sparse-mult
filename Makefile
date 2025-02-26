CC = mpicc
CFLAGS = -lmpi -lm
OFLAGS = -g -Wall -Wextra -O3 -march=native -funroll-loops
LFLAGS = -Iinclude

SRC_DIR = src
INCLUDE_DIR = include
OBJ_DIR = obj

FILES = main.c csr.c kernel.c utils.c

NPROCESS = 4

OBJS = $(FILES:%.c=$(OBJ_DIR)/%.o)

all: smult

smult: $(OBJS)
	$(CC) $(CFLAGS) $(OFLAGS) $(OBJS) -o $@ $(LFLAGS)

run: smult
	mpirun -n $(NPROCESS) ./smult mtx/twotone.mtx 100

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CFLAGS) $(OFLAGS) -I$(INCLUDE_DIR) -c $< -o $@

clean:
	rm -rf smult $(OBJ_DIR)