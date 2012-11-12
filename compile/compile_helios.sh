module add bullxmpi
module load intel
export CXX=mpic++
export CXXFLAGS='-std=c++0x -O3 -xHOST -fno-alias -DMPICH_IGNORE_CXX_SEEK'
export CXXFLAGS='-std=c++0x -O3 -xHOST -fno-alias -DMPICH_IGNORE_CXX_SEEK'

./configure --enable-openmp --enable-mpi                                    \
            --with-hdf5dir=/csc/home3/philscher/hdf5-1.8.8/                 \
            --enable-fftw3 --with-fftw3dir=/csc/home3/philscher/fftw-3.3.1/ \
            --enable-petsc --with-petscdir=/csc/home3/philscher/PETSc-3.2/  \
            --enable-slepc --with-slepcdir=/csc/home3/philscher/SLEPc-3.2/  \
            --enable-papi  --with-papidir=/csc/home3/philscher/papi-4.4.0/  \
            && make clean && make

