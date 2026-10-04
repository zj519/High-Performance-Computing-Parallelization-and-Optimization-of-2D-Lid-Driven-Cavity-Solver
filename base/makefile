CC=g++
CXX = g++

CFLAGS = -std=c++11 -Wall -O2 -DOMPI_SKIP_MPICXX -I/usr/lib/x86_64-linux-gnu/openmpi/lib
LDLIBS = -llapack -lblas -lboost_program_options -lboost_unit_test_framework

TARGET = solver
HDR_LDC=LidDrivenCavity.h SolverCG.h
HDR_CG=SolverCG.h
HDR_LDCS=LidDrivenCavity.h

LID_TEST_TARGET = lid_test
LID_TEST_SOURCES = Lid_test.cpp

SOLVER_TEST_TARGET=solvercg_test
SOLVER_TEST_SOURCES= SolverCG_test.cpp
DOXY      = doxygen
DOXYFLAGS = Doxyfile
default: all

all: $(TARGET) doc

LidDrivenCavitySolver.o:LidDrivenCavitySolver.cpp $(HDR_LDCS)
	$(CXX) $(CFLAGS) -c LidDrivenCavitySolver.cpp

LidDrivenCavity.o:LidDrivenCavity.cpp $(HDR_LDC)
	$(CXX) $(CFLAGS) -c LidDrivenCavity.cpp

SolverCG.o:SolverCG.cpp $(HDR_CG)
	$(CXX) $(CFLAGS) -c SolverCG.cpp
$(TARGET): LidDrivenCavitySolver.o LidDrivenCavity.o SolverCG.o
	$(CXX) $(CFLAGS) $^ -o $@ $(LDLIBS)

$(LID_TEST_TARGET): $(LID_TEST_SOURCES) SolverCG.o LidDrivenCavity.o
	$(CXX) $(CFLAGS) $^ -o $@ $(LDLIBS) 


$(SOLVER_TEST_TARGET): $(SOLVER_TEST_SOURCES) SolverCG.o
	$(CXX) $(CFLAGS) $^ -o $@ $(LDLIBS)

test: $(LID_TEST_TARGET) $(SOLVER_TEST_TARGET)
	./$(LID_TEST_TARGET)
	./$(SOLVER_TEST_TARGET)
.PHONY: clean test
doc:
	$(DOXY) $(DOXYFLAGS)
clean: 
	rm -f *.o $(TARGET) $(SOLVER_TEST_TARGET) $(LID_TEST_TARGET)