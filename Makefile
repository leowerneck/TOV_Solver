CXX      = g++
CXXFLAGS = -Wall -march=native -O2

SRC = Polytropic_EOS__struct_initialization.C \
      Polytropic_EOS__lowlevel_functions.C \
      RK4.C \
      TOV_RHSs.C \
      TOV_Solver.C \
      Output_routines.C

OBJ = $(SRC:.C=.o)
EXE = TOV_solver

all: $(EXE)

$(EXE): $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $(EXE)

$(OBJ): %.o: %.C
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(EXE)

veryclean: clean
	rm -f *.dat *.png
