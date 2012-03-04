bin_PROGRAMS = helios

gkc_SOURCES = main.cpp Helios.cpp FileIO.cpp Fields.cpp Fields/FieldsFFT.cpp Setup.cpp Grid.cpp  \
                 Parallel.cpp Vlasov.cpp Event.cpp\
                 Collisions.cpp Fourier.cpp  Vlasov/Vlasov_Cilk.cpp Analysis.cpp  Init.cpp Control.cpp\
      		 Plasma.cpp Special/LA.cpp Timing.cpp TestParticle.cpp TimeIntegration.cpp \
 	         FFTSolver.cpp
 
# compile FFT
gkc_SOURCES += FunctionParser/fparser.cc

AM_CPPFLAGS     = -ggdb -O3 -Wall -Wno-unknown-pragmas  -I$(DIR_HDF5)/include/ -I$(DIR_BLITZ)  -I./FFT -I./Geometry -I./Special -I./TimeIntegration -I./Fields
gkc_LDFLAGS  = -L$(DIR_HDF5)/lib/ -L$(DIR_BLITZ)/lib
gkc_LDADD    =   -lblitz $(DIR_HDF5)/lib/libhdf5_hl.a $(DIR_HDF5)/lib/libhdf5.a  -lz 
 
EXTRA_DIST      =  Helios.h FileIO.h Vlasov/Vlasov_Blitz.h Vlasov/Vlasov_Cilk.h Vlasov.h Fields.h Fields/FieldsFFT.h\
                   Fields/FieldsFEM.h Setup.h Grid.h Parallel.h Init.h \
	           Global.h Setup.h Analysis.h  Collisions.h Geometry  FFTSolver.h FEMSolver.h Fourier.h FunctionParser \
	           helios.doxygen  Plasma.h FFT/ Timing.h Integration.h GaussWeights.h FieldsHypre.h Control.h Geometry.h Benchmark.h FileAttr.h \
	           FileAttr.h TimeIntegration.h Visualization.h LA.h

if SHARED
gkc_LDADD +=  -lblitz  -lhdf5_hl -lhdf5 -lz 
endif

if STATIC
gkc_LDADD += $(DIR_BLITZ)/lib/liblitz.a   $(DIR_HDF5)/lib/libhdf5_hl.a $(DIR_HDF5)/lib/libhdf5.a  -lz -lgfortran
endif

if PETSC                 
AM_CPPFLAGS  +=  -I$(PETSC_DIR)/include
gkc_SOURCES += TimeIntegration/TimeIntegration_PETSc.cpp Special/PETScMatrixVector.cpp Fields/FieldsPETSc.cpp
gkc_LDADD   += -L$(PETSC_DIR)/lib -lpetsc 

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

if OPENMP
gkc_LDADD += -lfftw3_threads
endif
gkc_LDADD += -liomp5

if MPI

if MPIP
gkc_LDADD += -lmpiP -lunwind
endif

gkc_LDADD += -lmpi -lmpi_cxx -lmpi_f90 -lmpi_f77 

AM_CPPFLAGS  += -I/usr/include/mpi 
endif





if P3DFFT
gkc_SOURCES += FFT/FFTSolver_p3dfft.cpp
gkc_LDFLAGS +=  $(DIR_P3DFFT)/build/libp3dfft.a 
AM_CPPFLAGS  += -I$(DIR_P3DFFT)/include 
gkc_LDADD +=  $(DIR_P3DFFT)/build/libp3dfft.a -lgfortran
AM_CPPFLAGS  +=   -I/home/b/b30683/fftw-3.3alpha1/include
endif


if FFTW3
AM_CPPFLAGS  +=  -I$(DIR_FFTW3)/include
gkc_SOURCES += FFT/FFTSolver_fftw3.cpp 
gkc_LDADD   += -L$(DIR_FFTW3)/lib -lfftw3

if FFTW3MPI
gkc_LDADD   += -L$(DIR_FFTW3MPI)/lib -lfftw3_mpi -lfftw3_omp

endif
endif


if FFTW2
AM_CPPFLAGS  +=  -I$(DIR_FFTW2)/include
gkc_SOURCES += FFT/FFTSolver_fftw2.cpp
gkc_LDADD   += -L$(DIR_FFTW2)/lib  -lrfftw -lfftw
endif

if FFTW2MPI
AM_CPPFLAGS  +=  -I$(DIR_FFTW2MPI)/include
gkc_SOURCES += FFT/FFTSolver_fftw2_mpi.cpp
gkc_LDADD   += -L$(DIR_FFTW2MPI)/lib -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
endif

############################### Set Extra Libraries  Libraries #########
if HYPRE
gkc_SOURCES += Fields/FieldsHypre.cpp
AM_CPPFLAGS  += -I$(DIR_HYPRE)/include  
gkc_LDADD +=  -L$(DIR_HYPRE) -L$(DIR_HYPRE)/lib/ -lHYPRE
endif

if LIBMESH
AM_CPPFLAGS  += -I-I$(DIR_FEM)/
gkc_SOURCES += Fields/FieldsFEM.cpp Fields/FEMSolver.cpp
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

AM_CPPFLAGS += -I/home/philscher/Code/libmesh-0.7.0.3/include/base -I/home/philscher/Code/libmesh-0.7.0.3/include/enums -I/home/philscher/Code/libmesh-0.7.0.3/include/error_estimation -I/home/philscher/Code/libmesh-0.7.0.3/include/fe -I/home/philscher/Code/libmesh-0.7.0.3/include/geom -I/home/philscher/Code/libmesh-0.7.0.3/include/mesh -I/home/philscher/Code/libmesh-0.7.0.3/include/numerics -I/home/philscher/Code/libmesh-0.7.0.3/include/parallel -I/home/philscher/Code/libmesh-0.7.0.3/include/partitioning -I/home/philscher/Code/libmesh-0.7.0.3/include/quadrature -I/home/philscher/Code/libmesh-0.7.0.3/include/solvers -I/home/philscher/Code/libmesh-0.7.0.3/include/systems -I/home/philscher/Code/libmesh-0.7.0.3/include/utils -I/home/philscher/Code/libmesh-0.7.0.3/contrib/boost -I/home/philscher/Code/libmesh-0.7.0.3/contrib/laspack  -I/home/philscher/Code/libmesh-0.7.0.3/contrib/metis/Lib -I/home/philscher/Code/libmesh-0.7.0.3/contrib/sfcurves -I/home/philscher/Code/libmesh-0.7.0.3/contrib/gzstream -I/home/philscher/Code/libmesh-0.7.0.3/contrib/tecplot/include -I/home/philscher/Code/libmesh-0.7.0.3/contrib/gmv -I/usr/include/vtk-5.2 -I/home/philscher/Code/libmesh-0.7.0.3/contrib/tetgen -I/home/philscher/Code/libmesh-0.7.0.3/contrib/triangle -I/home/philscher/Code/libmesh-0.7.0.3/contrib/exodusii/Lib/include -I/home/philscher/Code/libmesh-0.7.0.3/contrib/netcdf/Lib -I/home/philscher/Code/libmesh-0.7.0.3/contrib/nemesis/Lib -I/home/philscher/Code/libmesh-0.7.0.3/contrib/libHilbert/include -I/usr/lib/slepc/include -I/usr/lib/petsc/include -I/usr/lib/petsc/bmake/linux-gnu-c-opt -I/usr/lib/petsc/linux-gnu-c-opt/include -I/usr/lib/openmpi/include -I/usr/lib/openmpi/lib -I/usr/include -I/usr/lib/openmpi/include -I/usr/lib/openmpi/lib -I/usr/include/spooles -I/usr/include -I/usr/include/superlu -I/usr/include/suitesparse -I/usr/include/scotch -I/usr/lib/openmpi/include -I/usr/lib/openmpi/lib -I/home/philscher/Code/hdf5mpi/include/ -std=c++0x -MT 

gkc_LDADD +=  -lz -lgfortran -lmpi -lmpi_cxx -lmpi_f90 -lmpi_f77 /home/philscher/Code/p3dfft//build/libp3dfft.a -lgfortran -lfftw3_mpi -lfftw3 /home/philscher/Code/libmesh-0.7.0.3/lib/x86_64-unknown-linux-gnu_opt/libmesh.so /home/philscher/Code/libmesh-0.7.0.3/contrib/lib/x86_64-unknown-linux-gnu_opt/liblaspack.so   /home/philscher/Code/libmesh-0.7.0.3/contrib/lib/x86_64-unknown-linux-gnu_opt/libsfcurves.so /home/philscher/Code/libmesh-0.7.0.3/contrib/lib/x86_64-unknown-linux-gnu_opt/libgzstream.so -lz /home/philscher/Code/libmesh-0.7.0.3/contrib/tecplot/lib/x86_64-unknown-linux-gnu/tecio.a /home/philscher/Code/libmesh-0.7.0.3/contrib/lib/x86_64-unknown-linux-gnu_opt/libgmv.so -lvtkIO -lvtkCommon -lvtkFiltering /home/philscher/Code/libmesh-0.7.0.3/contrib/lib/x86_64-unknown-linux-gnu_opt/libtetgen.so /home/philscher/Code/libmesh-0.7.0.3/contrib/lib/x86_64-unknown-linux-gnu_opt/libtriangle.so /home/philscher/Code/libmesh-0.7.0.3/contrib/lib/x86_64-unknown-linux-gnu_opt/libexodusii.so /home/philscher/Code/libmesh-0.7.0.3/contrib/lib/x86_64-unknown-linux-gnu_opt/libnetcdf.so /home/philscher/Code/libmesh-0.7.0.3/contrib/lib/x86_64-unknown-linux-gnu_opt/libnemesis.so /home/philscher/Code/libmesh-0.7.0.3/contrib/lib/x86_64-unknown-linux-gnu_opt/libHilbert.so -L/usr/lib/slepc/linux-gnu-c-opt/lib -L/usr/lib/slepc/linux-gnu-c-opt/lib -lslepc -L/usr/lib -lparpack -larpack -L/usr/lib/petsc/linux-gnu-c-opt/lib -L/usr/lib/petsc/linux-gnu-c-opt/lib -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc -lX11 -L/usr/lib -lspooles -lHYPRE -lmpi_cxx -lstdc++ -lsuperlu -lumfpack -lamd -llapack -lblas -L/usr/lib/openmpi/lib  -lmpi -lopen-rte -lopen-pal -lnsl -lutil -lgcc_s -lpthread -lmpi_f90 -lmpi_f77 -lgfortranbegin -lgfortran -lm -lm -L/usr/lib/gcc/x86_64-linux-gnu -lm -lm -ldl -lmpi -lopen-rte -lopen-pal -lnsl -lutil -lgcc_s -lpthread -lmpi_f90 -lmpi_f77 -lgfortranbegin -lgfortran -lm -lm -lm -lm -lscotch -L/usr/lib/openmpi/lib -L/usr/lib/gcc/x86_64-linux-gnu/4.4.0 -L/lib -L/usr/lib/x86_64-linux-gnu  -lmpi -lopen-rte -lopen-pal -lnsl -lutil -lgcc_s -lpthread -lmpi_f90 -lmpi_f77 -lgfortranbegin -lgfortran -lm -lm -L/usr/lib/gcc/x86_64-linux-gnu -lm -lm -lmpi_cxx -lstdc++ -lmpi_cxx -lstdc++  -lmpi -lopen-rte -lopen-pal -lnsl -lutil -lgcc_s -lpthread -L$(DIR_FEM)/lib/x86_64-unknown-linux-gnu_opt -lmesh 

endif


#if STATIC
#gkc_LDADD += $(DIR_FFTW)/fftw3_mpi.a
#endif

#if SHARED
#gkc_LDADD +=   -lfftw3_mpi 
#endif


#gkc_LDADD += -lfftw3 

AM_F90FLAGS = -I.. -I. -O3 
