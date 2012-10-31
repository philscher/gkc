#
HDF5_DIR=/opt/app/hdf5/1.8.8/intel-12.1/
SO =  main.o MHD.o DataOutput.o Input.o FunctionParser/fparser.o FunctionParser/fpoptimizer.cc 
#CXX = g++ 
CXX = icc
#CXXFLAGS = -fp-model fast=2 -fast -openmp -I$(HDF5_DIR)/include
CXXFLAGS = -fast -openmp -I$(HDF5_DIR)/include
LDFLAGS = -L$(HDF5_DIR)/lib -lhdf5 -lhdf5_hl -liomp5 -lm

EXE     = rmhd2 
#
.o:
	$(CXX) -c $<
$(EXE) : $(SO)
	$(CXX) -lm -o $(EXE) $(SO) $(LDFLAGS) 
	# other wise linkage error $(CXX) -o $(EXE) $(SO) $(LDFLAGS) 
main.o : FunctionParser/fparser.cc FunctionParser/fpoptimizer.cc DataOutput.cpp Input.cpp 
clean:
	rm -f *.o rmhd2



DataOutput.o  : DataOutput.cpp DataOutput.h
	               $(CXX) $(CXXFLAGS) -c DataOutput.cpp
$(EXE) : $(SO)
Input.o       : Input.cpp Input.h
	               $(CXX) $(CXXFLAGS) -c Input.cpp
FunctionParser/fpoptimizer.o : FunctionParser/fpoptimizer.cc
	               $(CXX) $(CXXFLAGS) -c FunctionParser/fpoptimizer.cpp
fparser.o     : FunctionParser/fparser.cc 

MHD.o     : MHD.cpp MHD.h 
