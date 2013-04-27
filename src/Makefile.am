bin_PROGRAMS = gkc

## Include source files

gkc_SOURCES = \
   main.cpp GKC.cpp FileIO.cpp Fields/Fields.cpp Fields/FieldsFFT.cpp Setup.cpp Grid.cpp  \
   Parallel/Parallel.cpp Vlasov/Vlasov.cpp Analysis/Event.cpp Vlasov/Vlasov_Island.cpp\
   Vlasov/Vlasov_Cilk.cpp Vlasov/Vlasov_Aux.cpp Vlasov/Vlasov_Optim.cpp \
   Analysis/Diagnostics.cpp  Init.cpp Control.cpp Analysis/Moments.cpp \
   Plasma.cpp Special/LA.cpp TimeIntegration/Timing.cpp Analysis/TestParticle.cpp TimeIntegration/TimeIntegration.cpp \
   FFTSolver/FFTSolver.cpp Visualization/Visualization_Data.cpp Benchmark/Benchmark_PAPI.cpp \
   Collisions/LenardBernstein.cpp Tools/TermColor.cpp TimeIntegration/ScanLinearModes.cpp \
   TimeIntegration/ScanPoloidalEigen.cpp Collisions/PitchAngle.cpp Benchmark/Benchmark_PMPI.cpp

## Include corresponding header files

gkc_SOURCES += \
   GKC.h FileIO.h Fields/Fields.h Fields/FieldsFFT.h Setup.h Grid.h  \
   Parallel/Parallel.h Vlasov/Vlasov.h Analysis/Event.h Vlasov/Vlasov_Island.h\
   Vlasov/Vlasov_Cilk.h Vlasov/Vlasov_Aux.h Vlasov/Vlasov_Optim.h \
   Analysis/Diagnostics.h Analysis/Moments.h Init.h Control.h \
   Plasma.h Special/LA.h TimeIntegration/Timing.h Analysis/TestParticle.h TimeIntegration/TimeIntegration.h \
   FFTSolver/FFTSolver.h \
   Visualization/Visualization.h Visualization/Visualization_Data.h \
   Analysis/Benchmark.h Global.h \
   Collisions/Collisions.h Collisions/LenardBernstein.h \
   Geometry/Geometry.h Geometry/Geometry2D.h Geometry/GeometryShear.h \
   Geometry/GeometrySlab.h Geometry/GeometrySA.h Geometry/GeometryCHEASE.h\
   Tools/System.h Tools/TermColor.h Special/SpecialMath.h Special/HermitePoly.h Tools/Tools.h Special/Vector3D.h \
   TimeIntegration/ScanLinearModes.h TimeIntegration/ScanPoloidalEigen.h Collisions/PitchAngle.h
 
# Integration sub-module
gkc_SOURCES += \
   Special/Integrate/GaussChebychevFirstWeights.h Special/Integrate/GaussChebychevSecondWeights.h \
   Special/Integrate/GaussHermiteWeights.h        Special/Integrate/RectangleWeights.h \
   Special/Integrate/GaussLaguerreWeights.h       Special/Integrate/GaussLegendreWeights.h \
   Special/Integrate/GaussRadauWeights.h      \
   Special/Integrate/GaussWeights.py              Special/Integrate/Integrate.h\
   Special/Interpolate/LinearInterpolate.h Special/RootFinding.h

# Interpolation sub-module
gkc_SOURCES += Special/Interpolate/HermiteInterpolation.h 

# Matrix module
gkc_SOURCES += Matrix/Matrix.h Matrix/MatrixPETSc.h Matrix/MatrixSolver.h \
               Matrix/SpecialMatrix.h Matrix/MatrixPETSc.cpp

EXTRA_DIST =   external

gkc_SOURCES += external/FunctionParser/fparser.cc external/SFL/SFL.cpp 

# Object files are created in subdirectories
AUTOMAKE_OPTIONS = subdir-objects

AM_CPPFLAGS  = -ggdb -fno-alias -align -fast -Wall -Wno-unknown-pragmas -I./external \
					-I./FFTSolver -I./Geometry -I./Special -I./TimeIntegration -I./Fields
 

#### Link 3rd party libraries

##gkc_LDADD    = -lstdc++

# HDF-5

# Add for HDF-5 file
gkc_LDADD = $(DIR_HDF5)/lib/libhdf5_hl.a $(DIR_HDF5)/lib/libhdf5.a  -lz 
gkc_LDFLAGS  = -L$(DIR_HDF5)/lib/
AM_CPPFLAGS += -I$(DIR_HDF5)/include/

if SHARED
gkc_LDADD += -lhdf5_hl -lhdf5 -lz 
endif

if STATIC
gkc_LDADD += $(DIR_HDF5)/lib/libhdf5_hl.a $(DIR_HDF5)/lib/libhdf5.a  -lz -lgfortran
endif

if PAPI
AM_CPPFLAGS += -I$(DIR_PAPI)/include
gkc_LDADD   += -L$(DIR_PAPI)/lib/ -lpapi -lrt -lpfm
endif

if MPI

if MPIP
gkc_LDADD += -lmpiP -lunwind
endif

endif

if OPENMP
AM_CPPFLAGS += -fopenmp
gkc_LDADD   += -liomp5
endif

##################### FFTW3 setup ################
if FFTW3

AM_CPPFLAGS  += -I$(DIR_FFTW3)/include
gkc_SOURCES  +=  FFTSolver/FFTSolver_fftw3.h FFTSolver/FFTSolver_fftw3.cpp 


#if STATIC
if MPI
   gkc_LDADD += $(DIR_FFTW3)/lib/libfftw3_mpi.a
endif
if OPENMP
   gkc_LDADD += $(DIR_FFTW3)/lib/libfftw3_omp.a
endif
   gkc_LDADD += $(DIR_FFTW3)/lib/libfftw3.a

#else
#   gkc_LDADD += -L$(DIR_FFTW3)/lib/ -lfftw3
#if MPI
#   gkc_LDADD += -lfftw3_mpi
#endif
#if OPENMP
#   gkc_LDADD += -lfftw3_omp
#endif

#endif
endif

################# SELPc (for Eigenvalue calculations) #######

if SLEPC
AM_CPPFLAGS += -I$(DIR_SLEPC)/include
gkc_SOURCES += Eigenvalue/Eigenvalue_SLEPc.cpp Eigenvalue/Eigenvalue_SLEPc.h Eigenvalue/Eigenvalue.h
gkc_LDADD   += -L$(DIR_SLEPC)/lib  -lslepc 
endif

if PETSC                 
AM_CPPFLAGS += -I$(DIR_PETSC)/include
gkc_SOURCES += TimeIntegration/TimeIntegration_PETSc.h TimeIntegration/TimeIntegration_PETSc.cpp\
     Special/PETScMatrixVector.h Special/PETScMatrixVector.cpp \
     Fields/FieldsHermite.h Fields/FieldsHermite.cpp 

gkc_LDADD   += -L$(DIR_PETSC)/lib -lpetsc 

endif

# Link MKL after fftw-3 otherwise will be overwritten (MKL does not include many transforms and thus returns null planner)
if MKL
AM_CPPFLAGS += -DMKL_ILP64 -I$(DIR_MKL)/include
gkc_LDFLAGS += -Wl,-rpath,$(DIR_MKL)/lib/intel64 
gkc_LDADD   += -L$(DIR_MKL)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lX11
endif

# Atlas packages
#gkc_LDADD   += -L/usr/lib64/atlas/  -llapack -lblas


