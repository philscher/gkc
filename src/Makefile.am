bin_PROGRAMS = gkc

## Include source files

gkc_SOURCES = main.cpp GKC.cpp FileIO.cpp Fields/Fields.cpp Fields/FieldsFFT.cpp Setup.cpp Grid.cpp  \
              Parallel/Parallel.cpp Vlasov/Vlasov.cpp Analysis/Event.cpp Vlasov/Vlasov_Island.cpp\
              Vlasov/Vlasov_Cilk.cpp Vlasov/Vlasov_Aux.cpp Vlasov/Vlasov_Optim.cpp \
				  Analysis/Analysis.cpp  Init.cpp Control.cpp \
              Plasma.cpp Special/LA.cpp TimeIntegration/Timing.cpp Analysis/TestParticle.cpp TimeIntegration/TimeIntegration.cpp \
              FFTSolver/FFTSolver.cpp Visualization/Visualization_Data.cpp Analysis/Benchmark.cpp \
				  Collisions/LenardBernstein.cpp



## Include corresponding header files

gkc_SOURCES += GKC.h FileIO.h Fields/Fields.h Fields/FieldsFFT.h Setup.h Grid.h  \
               Parallel/Parallel.h Vlasov/Vlasov.h Analysis/Event.h Vlasov/Vlasov_Island.h\
               Vlasov/Vlasov_Cilk.h Vlasov/Vlasov_Aux.h Vlasov/Vlasov_Optim.h \
		       Analysis/Analysis.h  Init.h Control.h \
               Plasma.h Special/LA.h TimeIntegration/Timing.h Analysis/TestParticle.h TimeIntegration/TimeIntegration.h \
               FFTSolver/FFTSolver.h \
			   Visualization/Visualization.h Visualization/Visualization_Data.h \
			   Analysis/Benchmark.h Global.h \
			   Collisions/Collisions.h Collisions/LenardBernstein.h \
			   Geometry/Geometry.h Geometry/Geometry2D.h Geometry/GeometryShear.h \
               Geometry/GeometrySlab.h Geometry/GeometrySA.h Geometry/GeometryCHEASE.h\
			   System.h Special/SpecialMath.h Special/HermitePoly.h Tools.h Special/Vector3D.h 
	
# Integration sub-module
gkc_SOURCES += Special/Integrate/GaussChebychevFirstWeights.h Special/Integrate/GaussChebychevSecondWeights.h \
				   Special/Integrate/GaussHermiteWeights.h        Special/Integrate/GaussLaguerInitializer.h \
					Special/Integrate/GaussLaguerreWeights.h       Special/Integrate/GaussLegendreWeights.h \
					Special/Integrate/GaussRadauWeights.h          Special/Integrate/GaussWeights.h \
					Special/Integrate/GaussWeights.py              Special/Integrate/Integrate.h


# Interpolation sub-module
gkc_SOURCES += Special/HermiteInterpolation.h 

#gkc_SOURCES += external/SHDF5/FileAttr.h external/SHDF5/TableAttr.h
#gkc_SOURCES += external/allocate.h external/Array.h 

gkc_SOURCES += Matrix/Matrix.h Matrix/MatrixPETSc.h Matrix/MatrixSolver.h \
					Matrix/SpecialMatrix.h Matrix/MatrixPETSc.cpp

EXTRA_DIST =   gkc.doxygen external

gkc_SOURCES += external/FunctionParser/fparser.cc external/SFL/SFL.cpp 

#gkc_SOURCES += external/FunctionParser/fparser.cc external/FunctionParser/fp_identifier_parser.inc 
#					external/FunctionParser/fparser.hh external/SFL/SFL.cpp external/SFL/SFL.h 
#					external/FunctionParser/fpconfig.hh	external/FunctionParser/fptypes.hh




AUTOMAKE_OPTIONS = subdir-objects


AM_CPPFLAGS  = -ggdb -O3 -Wall -Wno-unknown-pragmas -I./external  -I./FFTSolver -I./Geometry -I./Special -I./TimeIntegration -I./Fields
 

#Vlasov/Vlasov_Cilk.h Vlasov/Vlasov_Aux.h Vlasov/Vlasov_Optim.h Vlasov/Vlasov_Island.h Vlasov/Vlasov.

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
gkc_LDADD   += -L$(DIR_PAPI)/lib/ -lpapi -lrt
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
gkc_SOURCES  +=  FFTSolver/FFTSolver_fftw3.h FFTSolver/FFTSolver_fftw3.cpp 
gkc_LDADD    += -L$(DIR_FFTW3)/lib -lfftw3 -lfftw3_omp

if MPI
gkc_LDADD    += -L$(DIR_FFTW3MPI)/lib -lfftw3_mpi
endif
if OPENMP
gkc_LDADD    += -lfftw3_omp
endif
endif


################# SELPc (for Eigenvalue calculations) #######

if SLEPC
AM_CPPFLAGS +=  -I$(DIR_SLEPC)/include
gkc_SOURCES +=  Eigenvalue/Eigenvalue_SLEPc.cpp Eigenvalue/Eigenvalue_SLEPc.h Eigenvalue/Eigenvalue.h
gkc_LDADD   += -L$(DIR_SLEPC)/lib  -lslepc 
endif




################################# 3rd party libs #################
if PETSC                 
AM_CPPFLAGS +=  -I$(DIR_PETSC)/include
gkc_SOURCES += TimeIntegration/TimeIntegration_PETSc.h TimeIntegration/TimeIntegration_PETSc.cpp\
					Special/PETScMatrixVector.h Special/PETScMatrixVector.cpp \
					Fields/FieldsHermite.h Fields/FieldsHermite.cpp 
gkc_LDADD   += -L$(DIR_PETSC)/lib -lpetsc 

endif

# Atlas packages
# Link MKL after fftw-3 otherwise will be overwritten (MKL crashes for our transforms)
# And crashes for SLEPc (wtf)
#if MKL
#AM_CPPFLAGS += -DMKL_ILP64 -I$(MKLROOT)/include
#gkc_LDADD   += -L$(MKLROOT)/lib/intel64 -lmkl_cdft_core  -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_ilp64 -lpthread -lm -lX11
#gkc_LDFLAGS   += -Wl,-rpath,/home/opt/systemBC/intel/composer_xe_2011_sp1.9.293/mkl/lib/intel64 
#gkc_LDADD   += -L/home/opt/systemBC/intel/composer_xe_2011_sp1.9.293/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lX11
#gkc_LDADD   +=   -L$(MKLROOT)/lib/intel64 -lmkl_rt -lpthread -lm
#gkc_LDFLAGS += $(MKLROOT)/lib/intel64/libmkl_cdft_core.a $(MKLROOT)/lib/intel64/libmkl_intel_ilp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_ilp64.a -lpthread -lm
#gkc_LDADD   += -L/usr/lib64/atlas/  -llapack -lblas
#endif
#







