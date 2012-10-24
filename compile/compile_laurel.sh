# Ignore SEEK Problem
#PETSc -3.3 installation --prefix=/home/b/b30683/Packages/PETSc-3.3/ --with-scalar-type=complex --with-shared-libraries --with-mpi=yes --download-superlu=yes --download-scalapack --download-blacs --with-c++-support=1 --download-superlu-dist=yes --with-blas-lapack-dir=/opt/app/intel/mkl/lib/intel64/ --with-memalign=64
# Had to hand compile all fortran modules in src directory (roughly 7 ) using comilation command as given by executeing the makefile, and then
# rename the makefile, so that it is not compiled again, all modules have to be moved to include directory so they can be accesses
# ignore compilation error with SLEPc and directly copy linux=archdata
#./configure --enable-openmp --with-hdf5dir=/LARGE0/gr10140/Packages/hdf5-1.8.8/ --with-blitzdir=/LARGE0/gr10140/Packages/blitz-0.10/ --enable-mpi --enable-fftw3 --with-fftw3dir=/LARGE0/gr10140/Packages/fftw-3.3.2/  --enable-petsc --with-petscdir=/LARGE0/gr10140/Packages/PETSc-3.3/  --enable-slepc --with-slepcdir=/LARGE0/gr10140/Packages/SLEPc-3.3/ CXXFLAGS=-DMPICH_IGNORE_CXX_SEEK LDFLAGS=-L/opt/app/intel/impi/4.0.3.008/lib64/ && make clean && make
#./configure --enable-openmp --with-hdf5dir=/LARGE0/gr10140/Packages/hdf5-1.8.8/ --enable-mpi --enable-fftw3 --with-fftw3dir=/LARGE0/gr10140/Packages/fftw-3.3.2/  --enable-petsc --with-petscdir=/LARGE0/gr10140/Packages/PETSc-3.3/  --enable-slepc --with-slepcdir=/LARGE0/gr10140/Packages/SLEPc-3.3/ CXXFLAGS=-DMPICH_IGNORE_CXX_SEEK && make clean && make
#./configure --enable-openmp --with-hdf5dir=/opt/app/hdf5-parallel/1.8.8/intel-12.1/ --enable-mpi --enable-fftw3 --with-fftw3dir=/opt/app/fftw/3.3.1/intel-12.1/  --enable-petsc --with-petscdir=/LARGE0/gr10140/Packages/PETSc-3.3/  --enable-slepc --with-slepcdir=/LARGE0/gr10140/Packages/SLEPc-3.3/ CXXFLAGS=-DMPICH_IGNORE_CXX_SEEK && make clean && make

# Need to load the corresponding modules
module load intel
#module load intel/13.0.0
export CXX=mpiicc
export CXXFLAGS='-std=c++0x -O3 -xHOST -fno-alias -DMPICH_IGNORE_CXX_SEEK'
./configure --enable-openmp --enable-mpi                                           \
            --with-hdf5dir=/opt/app/hdf5-parallel/1.8.8/intel-12.1/                \
            --enable-fftw3 --with-fftw3dir=/LARGE0/gr10140/Packages/fftw-3.3.2/    \
            --enable-petsc --with-petscdir=/opt/app/petsc/3.2_complex/intel-12.1/  \
            --enable-slepc --with-slepcdir=/opt/app/slepc/3.2_complex/intel-12.1/  \
            && make clean && make 
