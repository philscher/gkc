# Need to load the corresponding modules
module load intel
module load intel/mkl
#module load intel/13.0.0 (compiler crashes during comilation with internal error)
module load petsc-complex/3.2_intel-12.1
module load slepc-complex/3.2_intel-12.1

export CXX=mpiicc
export CXXFLAGS='-std=c++0x -O3 -xHOST -fno-alias -DMPICH_IGNORE_CXX_SEEK'

# Configure with proper parameters
./configure --enable-mpi --enable-openmp                                           \
            --with-hdf5dir=/opt/app/hdf5-parallel/1.8.8/intel-12.1/                \
            --enable-fftw3 --with-fftw3dir=/LARGE0/gr10140/Packages/fftw-3.3.2/    \
            --enable-petsc --with-petscdir=$PETSC_DIR                              \
            --enable-slepc --with-slepcdir=$SLEPC_DIR                              \
            --enable-mkl   --with-mkldir=$MKLROOT                                  \
            && make clean && make 
