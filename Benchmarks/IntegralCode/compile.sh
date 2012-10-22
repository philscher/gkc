# Compile Special Mathfunction library, create library
g++  -Ofast -flto -march=native -funroll-loops -fPIC  -c SpecialMath.cpp  -o spmath.o
#g++ -shared -Wl,-soname,spmath.so -o libspmath.so  spmath.o

# Compile Cython functions
g++ -Ofast -flto -march=native -fopenmp -funroll-loops -fPIC -shared  -lrt -c iCode-func.cpp

# Create cython file, compile to python module
cython iCode.pyx

g++ -shared -fPIC -O3 -Wall -fno-strict-aliasing  -I/usr/include/python2.7 -o iCode.so iCode.c  iCode-func.o spmath.o -lgomp 

# Clean-up temporary files
rm *.o
rm iCode.c
rm -rf build

