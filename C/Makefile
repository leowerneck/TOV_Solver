CC       = gcc
CFLAGS   = -Wall -march=native -O2 -std=c99
LD_FLAGS = -lm

HDF5_DIR = /usr/local/opt/hdf5
HDF5_INC = -I$(HDF5_DIR)/include
HDF5_LIB = -L$(HDF5_DIR)/lib -lhdf5

SRC = eos_table_reader.c
OBJ = $(SRC:.c=.o)
EXE = tov_solver

$(EXE): $(OBJ)
	$(CC) $(CFLAGS) $(HDF5_INC) $(OBJ) -o $(EXE) $(HDF5_LIB) $(LD_FLAGS)

$(OBJ): %.o :%.c $(INC)
	$(CC) $(CFLAGS) $(HDF5_INC) -c $< -o $@

clean:
	rm -f $(EXE) $(OBJ)
