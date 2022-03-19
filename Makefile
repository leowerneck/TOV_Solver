CXX      = g++
CXXFLAGS = -Wall -march=native -O2

SRC = polytropic_eos__struct_initialization.cc \
      polytropic_eos__lowlevel_functions.cc \
      rk4.cc \
      tov_rhss.cc \
      tov_solver.cc \
      output_routines.cc

OBJ = $(SRC:.cc=.o)
EXE = tov_solver

all: $(EXE)

$(EXE): $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $(EXE)

$(OBJ): %.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(EXE)

veryclean: clean
	rm -f *.dat *.png
