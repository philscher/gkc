bin_PROGRAMS = gkc

gkc_SOURCES = main.cpp GKC.cpp FileIO.cpp Fields/Fields.cpp Fields/FieldsFFT.cpp Setup.cpp Grid.cpp  \
                 Parallel/Parallel.cpp Vlasov/Vlasov.cpp Analysis/Event.cpp Vlasov/Vlasov_Island.cpp\
                 Collisions/Collisions.cpp  Vlasov/Vlasov_Cilk.cpp Vlasov/Vlasov_Aux.cpp Vlasov/Vlasov_Optim.cpp Analysis/Analysis.cpp  Init.cpp Control.cpp\
                 Plasma.cpp Special/LA.cpp Timing.cpp Analysis/TestParticle.cpp TimeIntegration/TimeIntegration.cpp \
                FFTSolver/FFTSolver.cpp Visualization/Visualization_Data.cpp Analysis/Benchmark.cpp
 
AUTOMAKE_OPTIONS = subdir-objects
# compile FFT
gkc_SOURCES += external/FunctionParser/fparser.cc external/SFL/SFL.cpp

AM_CPPFLAGS     = -ggdb -O3 -Wall -Wno-unknown-pragmas  -I$(DIR_HDF5)/include/ -I./external  -I./FFTSolver -I./Geometry -I./Special -I./TimeIntegration -I./Fields
gkc_LDFLAGS  = -L$(DIR_HDF5)/lib/
gkc_LDADD    = $(DIR_HDF5)/lib/libhdf5_hl.a $(DIR_HDF5)/lib/libhdf5.a  -lz 
 
EXTRA_DIST      =  GKC.h FileIO.h Vlasov/Vlasov_Cilk.h Vlasov/Vlasov_Aux.h Vlasov/Vlasov_Optim.h Vlasov/Vlasov_Island.h Vlasov.h Fields.h Fields/FieldsFFT.h\
                    Setup.h Grid.h Parallel.h Init.h \
              Global.h Setup.h Analysis.h  Collisions.h Geometry  FFTSolver.h Fourier.h FunctionParser \
              gkc.doxygen  Plasma.h FFTSolver/ Timing.h Integration.h Special/GaussWeights.h Control.h Geometry.h Analysis/Benchmark.h \
              external/SHDF5/FileAttr.h external/SHDF5/TableAttr.h TimeIntegration.h Visualization.h Special/LA.h Special/HermitePoly.h


if SHARED
gkc_LDADD +=  -lhdf5_hl -lhdf5 -lz 
endif

if STATIC
gkc_LDADD += $(DIR_HDF5)/lib/libhdf5_hl.a $(DIR_HDF5)/lib/libhdf5.a  -lz -lgfortran
endif

if PAPI
AM_CPPFLAGS += -I$(DIR_PAPI)/include
gkc_LDADD   += -L$(DIR_PAPI)/lib/ -lpapi
endif


if MPI

if MPIP
gkc_LDADD += -lmpiP -lunwind
endif

#gkc_LDADD += -lmpi -lmpi_cxx -lmpi_f90 -lmpi_f77 
#AM_CPPFLAGS  += -I/usr/include/mpi  -std=c++0x 

endif

if OPENMP
AM_CPPFLAGS    += -fopenmp -DBZ_THREADSAFE
endif


################################# 3rd party libs #################

if PETSC                 
AM_CPPFLAGS    +=  -I$(DIR_PETSC)/include
gkc_SOURCES += TimeIntegration/TimeIntegration_PETSc.cpp Special/PETScMatrixVector.cpp Fields/FieldsHermite.cpp
#Fields/FieldsPETSc.cpp \
gkc_LDADD   += -L$(DIR_PETSC)/lib -lpetsc 
#gkc_LDADD   += -L$(DIR_PETSC)/lib -lpetsc -llapack  -lblas

if SLEPC
################# Eigenvalue calculations #######
AM_CPPFLAGS    +=  -I$(DIR_SLEPC)/include
gkc_SOURCES +=  Eigenvalue/Eigenvalue_SLEPc.cpp
gkc_LDADD   += -L$(DIR_SLEPC)/lib  -lslepc 
endif

endif
############################### Set Various FFT Libraries #########

#
# Definitions for the parallel version ...
#
#

# Intel stuff clean up
#gkc_LDADD += -L/home/philscher/Science/Packages/intel/lib/intel64  -liomp5 -lsvml -lirc -lifport -lifcore -lintlc -limf



##################### FFTW3 setup ################
if FFTW3

AM_CPPFLAGS    += -I$(DIR_FFTW3)/include
gkc_SOURCES +=  FFTSolver/FFTSolver_fftw3.cpp 
gkc_LDADD   += -L$(DIR_FFTW3)/lib -lfftw3
gkc_LDADD   += -lfftw3_omp

if MPI
gkc_LDADD   += -L$(DIR_FFTW3MPI)/lib -lfftw3_mpi
endif

if OPENMP
gkc_LDADD   += -lfftw3_omp
endif


endif

