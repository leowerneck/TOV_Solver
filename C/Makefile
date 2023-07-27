CC       = gcc
CFLAGS   = -Wall -march=native -O2 -std=c99 -I./include
LD_FLAGS = -lm

HDF5_DIR = /usr/local/opt/hdf5
HDF5_INC = -I$(HDF5_DIR)/include
HDF5_LIB = -L$(HDF5_DIR)/lib -lhdf5

SRC = tov_solver.c eos_poly_init_simple.c eos_poly_init_from_name.c eos_poly_funcs.c eos_table_reader.c eos_table_interp1d.c eos_table_interp3d.c eos_free.c
INC = tov_solver.h eos.h
EXE = tov_solver

OBJ := $(addprefix build/,   $(SRC:.c=.o))
SRC := $(addprefix src/,     $(SRC))
INC := $(addprefix include/, $(INC))

$(EXE): build $(OBJ)
	$(CC) $(CFLAGS) $(HDF5_INC) $(OBJ) -o $(EXE) $(HDF5_LIB) $(LD_FLAGS)

build:
	mkdir -p build

$(OBJ): build/%.o : src/%.c $(INC)
	$(CC) $(CFLAGS) $(HDF5_INC) -c $< -o $@

clean:
	rm -f $(EXE) $(OBJ)
