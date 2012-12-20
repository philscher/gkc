# Need to load the corresponding modules
module load intel
module load intel/mkl
#module load intel/13.0.0 (compiler crashes during comilation with internal error)

export CXX=mpiicc
export CXXFLAGS='-std=c++0x -O3 -xHOST -fno-alias -DMPICH_IGNORE_CXX_SEEK'

# Configure with proper parameters
./configure --enable-mpi --enable-openmp                                           \
            --with-hdf5dir=/opt/app/hdf5-parallel/1.8.8/intel-12.1/                \
            --enable-fftw3 --with-fftw3dir=/LARGE0/gr10140/Packages/fftw-3.3.2/    \
            --enable-petsc --with-petscdir=/opt/app/petsc/3.2_complex/intel-12.1/  \
            --enable-slepc --with-slepcdir=/opt/app/slepc/3.2_complex/intel-12.1/  \
            --enable-mkl   --with-mkldir=$MKLROOT                                  \
            && make clean && make 
