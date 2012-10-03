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
AM_CPPFLAGS  += -I/usr/include/mpi  -std=c++0x 

endif

if OPENMP
AM_CPPFLAGS    += -fopenmp -DBZ_THREADSAFE
endif


################################# 3rd party libs #################

if PETSC                 
AM_CPPFLAGS    +=  -I$(PETSC_DIR)/include
gkc_SOURCES += TimeIntegration/TimeIntegration_PETSc.cpp Special/PETScMatrixVector.cpp Fields/FieldsHermite.cpp
#Fields/FieldsPETSc.cpp \
gkc_LDADD   += -L$(PETSC_DIR)/lib -lpetsc 
#gkc_LDADD   += -L$(PETSC_DIR)/lib -lpetsc -llapack  -lblas

if SLEPC
################# Eigenvalue calculations #######
AM_CPPFLAGS    +=  -I$(SLEPC_DIR)/include
gkc_SOURCES +=  Eigenvalue/Eigenvalue_SLEPc.cpp
gkc_LDADD   += -L$(SLEPC_DIR)/lib  -lslepc 
endif

endif
############################### Set Various FFT Libraries #########

#
# Definitions for the parallel version ...
#
#

# Intel stuff clean up
gkc_LDADD += -L/home/philscher/Science/Packages/intel/lib/intel64  -liomp5 -lsvml -lirc -lifport -lifcore -lintlc -limf



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

##################### FFTW2 setup ################

#if FFTW2
#AM_CPPFLAGS  +=  -I$(DIR_FFTW2)/include
#gkc_SOURCES += FFT/FFTSolver_fftw2.cpp
#gkc_LDADD   += -L$(DIR_FFTW2)/lib  -lrfftw -lfftw
#
#if MPI
#AM_CPPFLAGS  +=  -I$(DIR_FFTW2MPI)/include
#gkc_SOURCES += FFT/FFTSolver_fftw2_mpi.cpp
#gkc_LDADD   += -L$(DIR_FFTW2MPI)/lib -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
#endif
#
#endif

############################### Set Extra Libraries  Libraries #########
if HYPRE
gkc_SOURCES += Fields/FieldsHypre.cpp
AM_CPPFLAGS  += -I$(DIR_HYPRE)/include  
gkc_LDADD +=  -L$(DIR_HYPRE) -L$(DIR_HYPRE)/lib/ -lHYPRE
endif

if LIBMESH
AM_CPPFLAGS  += -I-I$(DIR_FEM)/
gkc_SOURCES += Fields/FieldsFEM.cpp FEM/FEMSolver.cpp
include $(DIR_FEM)/Make.common
AM_CPPFLAGS += -I$(DIR_FEM)/include/base \
                   -I$(DIR_FEM)/include/enums \
                   -I$(DIR_FEM)/include/error_estimation \
                   -I$(DIR_FEM)/include/fe \
                   -I$(DIR_FEM)/include/systems \
                   -I$(DIR_FEM)/include/geom \
                   -I$(DIR_FEM)/include/mesh \
                   -I$(DIR_FEM)/include/numerics \
                   -I$(DIR_FEM)/include/parallel \
                   -I$(DIR_FEM)/include/partitioning \
                   -I$(DIR_FEM)/include/quadrature \
                   -I$(DIR_FEM)/include/solvers \
                   -I$(DIR_FEM)/include/utils

AM_CPPFLAGS += -I/home/philscher/Code/libmesh-0.7.0.3/include/base -I/home/philscher/Code/libmesh-0.7.0.3/include/enums -I/home/philscher/Code/libmesh-0.7.0.3/include/error_estimation -I/home/philscher/Code/libmesh-0.7.0.3/include/fe -I/home/philscher/Code/libmesh-0.7.0.3/include/geom -I/home/philscher/Code/libmesh-0.7.0.3/include/mesh -I/home/philscher/Code/libmesh-0.7.0.3/include/numerics -I/home/philscher/Code/libmesh-0.7.0.3/include/parallel -I/home/philscher/Code/libmesh-0.7.0.3/include/partitioning -I/home/philscher/Code/libmesh-0.7.0.3/include/quadrature -I/home/philscher/Code/libmesh-0.7.0.3/include/solvers -I/home/philscher/Code/libmesh-0.7.0.3/include/systems -I/home/philscher/Code/libmesh-0.7.0.3/include/utils -I/home/philscher/Code/libmesh-0.7.0.3/contrib/boost -I/home/philscher/Code/libmesh-0.7.0.3/contrib/laspack  -I/home/philscher/Code/libmesh-0.7.0.3/contrib/metis/Lib -I/home/philscher/Code/libmesh-0.7.0.3/contrib/sfcurves -I/home/philscher/Code/libmesh-0.7.0.3/contrib/gzstream -I/home/philscher/Code/libmesh-0.7.0.3/contrib/tecplot/include -I/home/philscher/Code/libmesh-0.7.0.3/contrib/gmv -I/usr/include/vtk-5.2 -I/home/philscher/Code/libmesh-0.7.0.3/contrib/tetgen -I/home/philscher/Code/libmesh-0.7.0.3/contrib/triangle -I/home/philscher/Code/libmesh-0.7.0.3/contrib/exodusii/Lib/include -I/home/philscher/Code/libmesh-0.7.0.3/contrib/netcdf/Lib -I/home/philscher/Code/libmesh-0.7.0.3/contrib/nemesis/Lib -I/home/philscher/Code/libmesh-0.7.0.3/contrib/libHilbert/include -I/usr/lib/slepc/include -I/usr/lib/petsc/include -I/usr/lib/petsc/bmake/linux-gnu-c-opt -I/usr/lib/petsc/linux-gnu-c-opt/include -I/usr/lib/openmpi/include -I/usr/lib/openmpi/lib -I/usr/include -I/usr/lib/openmpi/include -I/usr/lib/openmpi/lib -I/usr/include/spooles -I/usr/include  -I/usr/include/suitesparse -I/usr/include/scotch -I/usr/lib/openmpi/include -I/usr/lib/openmpi/lib -I/home/philscher/Code/hdf5mpi/include/ -std=c++0x -MT 

#gkc_LDADD +=  -lz -lgfortran -lmpi -lmpi_cxx -lmpi_f90 -lmpi_f77 /home/philscher/Code/p3dfft//build/libp3dfft.a -lgfortran -lfftw3_mpi -lfftw3 /home/philscher/Code/libmesh-0.7.0.3/lib/x86_64-unknown-linux-gnu_opt/libmesh.so /home/philscher/Code/libmesh-0.7.0.3/contrib/lib/x86_64-unknown-linux-gnu_opt/liblaspack.so   /home/philscher/Code/libmesh-0.7.0.3/contrib/lib/x86_64-unknown-linux-gnu_opt/libsfcurves.so /home/philscher/Code/libmesh-0.7.0.3/contrib/lib/x86_64-unknown-linux-gnu_opt/libgzstream.so -lz /home/philscher/Code/libmesh-0.7.0.3/contrib/tecplot/lib/x86_64-unknown-linux-gnu/tecio.a /home/philscher/Code/libmesh-0.7.0.3/contrib/lib/x86_64-unknown-linux-gnu_opt/libgmv.so -lvtkIO -lvtkCommon -lvtkFiltering /home/philscher/Code/libmesh-0.7.0.3/contrib/lib/x86_64-unknown-linux-gnu_opt/libtetgen.so /home/philscher/Code/libmesh-0.7.0.3/contrib/lib/x86_64-unknown-linux-gnu_opt/libtriangle.so /home/philscher/Code/libmesh-0.7.0.3/contrib/lib/x86_64-unknown-linux-gnu_opt/libexodusii.so /home/philscher/Code/libmesh-0.7.0.3/contrib/lib/x86_64-unknown-linux-gnu_opt/libnetcdf.so /home/philscher/Code/libmesh-0.7.0.3/contrib/lib/x86_64-unknown-linux-gnu_opt/libnemesis.so /home/philscher/Code/libmesh-0.7.0.3/contrib/lib/x86_64-unknown-linux-gnu_opt/libHilbert.so -L/usr/lib/slepc/linux-gnu-c-opt/lib -L/usr/lib/slepc/linux-gnu-c-opt/lib -lslepc -L/usr/lib -lparpack -larpack -L/usr/lib/petsc/linux-gnu-c-opt/lib -L/usr/lib/petsc/linux-gnu-c-opt/lib -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc -lX11 -L/usr/lib -lspooles -lHYPRE -lmpi_cxx -lstdc++  -lumfpack -lamd -llapack -lblas -L/usr/lib/openmpi/lib  -lmpi -lopen-rte -lopen-pal -lnsl -lutil -lgcc_s -lpthread -lmpi_f90 -lmpi_f77 -lgfortranbegin -lgfortran -lm -lm -L/usr/lib/gcc/x86_64-linux-gnu -lm -lm -ldl -lmpi -lopen-rte -lopen-pal -lnsl -lutil -lgcc_s -lpthread -lmpi_f90 -lmpi_f77 -lgfortranbegin -lgfortran -lm -lm -lm -lm -lscotch -L/usr/lib/openmpi/lib -L/usr/lib/gcc/x86_64-linux-gnu/4.4.0 -L/lib -L/usr/lib/x86_64-linux-gnu  -lmpi -lopen-rte -lopen-pal -lnsl -lutil -lgcc_s -lpthread -lmpi_f90 -lmpi_f77 -lgfortranbegin -lgfortran -lm -lm -L/usr/lib/gcc/x86_64-linux-gnu -lm -lm -lmpi_cxx -lstdc++ -lmpi_cxx -lstdc++  -lmpi -lopen-rte -lopen-pal -lnsl -lutil -lgcc_s -lpthread -L$(DIR_FEM)/lib/x86_64-unknown-linux-gnu_opt -lmesh 

endif


#AM_F90FLAGS = -I.. -I. -O3 
