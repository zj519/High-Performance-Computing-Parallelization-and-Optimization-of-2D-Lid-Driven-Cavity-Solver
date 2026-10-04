CC = mpicxx
CXX = mpicxx
CXXFLAGS = -std=c++17 -Wall -O3 -DOMPI_SKIP_MPICXX -I/usr/lib/x86_64-linux-gnu/openmpi/lib
LDLIBS = -llapack -lblas -lboost_program_options #-lboost_timer

TARGET = solver

HDR_LDC=LidDrivenCavity.h SolverCG.h Communicator.h 
HDR_CG=SolverCG.h Communicator.h
HDR_LDCS=LidDrivenCavity.h 
HDR_CM=Communicator.h
TARGET_OBJECT= LidDrivenCavitySolver.o LidDrivenCavity.o SolverCG.o Communicator.o
DOXY      = doxygen
DOXYFLAGS = Doxyfile
default: all doc
all: $(TARGET)

LidDrivenCavitySolver.o:LidDrivenCavitySolver.cpp $(HDR_LDCS)
	$(CXX) $(CFLAGS) -c LidDrivenCavitySolver.cpp

LidDrivenCavity.o:LidDrivenCavity.cpp $(HDR_LDC)
	$(CXX) $(CFLAGS) -c LidDrivenCavity.cpp

SolverCG.o:SolverCG.cpp $(HDR_CG)
	$(CXX) $(CFLAGS) -c SolverCG.cpp*

Communicator.o:Communicator.cpp $(HDR_CM)
	$(CXX) $(CFLAGS) -c Communicator.cpp

$(TARGET): $(TARGET_OBJECT)
	$(CXX) $(CFLAGS) $^ -o $@ $(LDLIBS)

.PHONY: clean

doc:
	$(DOXY) $(DOXYFLAGS)

clean:
	-rm -f *.o $(TARGET) 