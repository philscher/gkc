bin_PROGRAMS = gkc

gkc_SOURCES = main.cpp GKC.cpp FileIO.cpp Fields/Fields.cpp Fields/FieldsFFT.cpp Setup.cpp Grid.cpp  \
              Parallel/Parallel.cpp Vlasov/Vlasov.cpp Analysis/Event.cpp Vlasov/Vlasov_Island.cpp\
              Vlasov/Vlasov_Cilk.cpp Vlasov/Vlasov_Aux.cpp Vlasov/Vlasov_Optim.cpp \
				  Analysis/Analysis.cpp  Init.cpp Control.cpp \
              Plasma.cpp Special/LA.cpp Timing.cpp Analysis/TestParticle.cpp TimeIntegration/TimeIntegration.cpp \
              FFTSolver/FFTSolver.cpp Visualization/Visualization_Data.cpp Analysis/Benchmark.cpp \
				  Collisions/LenardBernstein.cpp
 
AUTOMAKE_OPTIONS = subdir-objects

gkc_SOURCES += external/FunctionParser/fparser.cc external/SFL/SFL.cpp

AM_CPPFLAGS  = -ggdb -O3 -Wall -Wno-unknown-pragmas -I./external  -I./FFTSolver -I./Geometry -I./Special -I./TimeIntegration -I./Fields
 
EXTRA_DIST =  GKC.h FileIO.h Vlasov/Vlasov_Cilk.h Vlasov/Vlasov_Aux.h Vlasov/Vlasov_Optim.h Vlasov/Vlasov_Island.h Vlasov.h Fields.h Fields/FieldsFFT.h\
              Setup.h Grid.h Parallel.h Init.h \
              Global.h Setup.h Analysis.h  Collisions.h Geometry  FFTSolver.h Fourier.h FunctionParser \
              gkc.doxygen  Plasma.h FFTSolver/ Timing.h Integration.h Special/GaussWeights.h Control.h Geometry.h Analysis/Benchmark.h \
              external/SHDF5/FileAttr.h external/SHDF5/TableAttr.h TimeIntegration.h Visualization.h Special/LA.h Special/HermitePoly.h


# Add for HDF-5 file
gkc_LDFLAGS  =   -L$(DIR_HDF5)/lib/
gkc_LDADD    =   $(DIR_HDF5)/lib/libhdf5_hl.a $(DIR_HDF5)/lib/libhdf5.a  -lz 
AM_CPPFLAGS += -I$(DIR_HDF5)/include/

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

endif

if OPENMP
AM_CPPFLAGS += -fopenmp
endif

##################### FFTW3 setup ################
if FFTW3
AM_CPPFLAGS  += -I$(DIR_FFTW3)/include
gkc_SOURCES  +=  FFTSolver/FFTSolver_fftw3.cpp 
gkc_LDADD    += -L$(DIR_FFTW3)/lib -lfftw3
gkc_LDADD    += -lfftw3_omp
if MPI
gkc_LDADD    += -L$(DIR_FFTW3MPI)/lib -lfftw3_mpi
endif
if OPENMP
gkc_LDADD    += -lfftw3_omp
endif
endif


################# SELPc (for Eigenvalue calculations) #######

if SLEPC
AM_CPPFLAGS    +=  -I$(DIR_SLEPC)/include
gkc_SOURCES +=  Eigenvalue/Eigenvalue_SLEPc.cpp
#gkc_LDADD   += -L$(DIR_SLEPC)/lib  -lslepc 
gkc_LDADD +=/opt/app/slepc/3.2_complex/intel-12.1/lib/libslepc.a
endif




################################# 3rd party libs #################
if PETSC                 
AM_CPPFLAGS +=  -I$(DIR_PETSC)/include
gkc_SOURCES += TimeIntegration/TimeIntegration_PETSc.cpp Special/PETScMatrixVector.cpp Fields/FieldsHermite.cpp
#gkc_LDADD   += -L$(DIR_PETSC)/lib -lpetsc 

# Use static linkage
gkc_LDADD +=/opt/app/petsc/3.2_complex/intel-12.1/lib/libpetsc.a

endif

# Atlas packages
# Link MKL after fftw-3 otherwise will be overwritten (MKL crashes for our transforms)
# And crashes for SLEPc (wtf)
#if MKL
AM_CPPFLAGS += -DMKL_ILP64 -I$(MKLROOT)/include
#gkc_LDADD   += -L$(MKLROOT)/lib/intel64 -lmkl_cdft_core  -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_ilp64 -lpthread -lm -lX11
gkc_LDFLAGS   += -Wl,-rpath,/home/opt/systemBC/intel/composer_xe_2011_sp1.9.293/mkl/lib/intel64 
gkc_LDADD   += -L/home/opt/systemBC/intel/composer_xe_2011_sp1.9.293/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lX11
#gkc_LDADD   +=   -L$(MKLROOT)/lib/intel64 -lmkl_rt -lpthread -lm
#gkc_LDFLAGS += $(MKLROOT)/lib/intel64/libmkl_cdft_core.a $(MKLROOT)/lib/intel64/libmkl_intel_ilp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_ilp64.a -lpthread -lm
#endif
#gkc_LDADD   += -L/usr/lib64/atlas/  -llapack -lblas





#ifdef PAPI
#AM_CPPFLAGS += -DMKL_ILP64 -I$(MKLROOT)/include
gkc_LDADD   += -L/home/b/b30683/papi-5.0.0/lib -lpapi
#endif


