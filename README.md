The GKC++ code
==============================

Description
-----------

GKC++ (GyroKinetic Code++) is a δf gyro-kinetic 
code. Currently features are :
      
  * Various geometries (sheared, shearless, toroidal)
  * multiple species support 
  * fully electromagnetic code 
  * Interfaces for various field and Vlasov solvers
  * Eigensolver interface using 
    [SLEPc](http://www.grycap.upv.es/slepc/)/[PETSc](http://www.mcs.anl.gov/petsc/)
  * Decomposition into eigenvectors using  [elemental](https://code.google.com/p/elemental/)
  * Preliminary support for global simulations
  * 3-stage parallelization (Vectorization, OpenMP and MPI) 
  * Data output using [HDF-5](www.hdfgroup.org/HDF5). 
  * Data analysis using 
    [ipython](http://ipython.org/)/[matplotlib](http://matplotlib.org/)/
    [scipy](http://www.scipy.org/)/[pytables](http://www.pytables.org/)
 
However, this code is still considered *experimental*
and benchmark to other codes is pending. You may want to consider using
other gyro-kinetic codes such as [GENE](http://www.ipp.mpg.de/~fsj/gene/),
[gkw](http://www.gkw.org.uk/) or [gs2](http://gs2.sourceforge.net/) for
production usage.

GKC++ is programmed in C++-11/[Cilk+](http://cilkplus.org/) 
which (for our opinion) greatly simplifies readability and handling of
multi-dimensional arrays. However, note that only a limited subset of compilers support
the Cilk+ extensions such as

* [Intel Compiler](http://software.intel.com/en-us/intel-compilers)
* [GCC-4.9] (http://gcc.gnu.org/wiki/cilkplus-merge)
* [Cilkplus/LLVM](http://cilkplus.github.io/) 

Why C++ but not Fortran
------------------------
  
Fortran is popular in numerical science communities due to
it's power of handling multi-dimensional arrays.
However, the lack of templates, function overloading, pointers
and classes (although some of these features are now supported
by the Fortran-03 standard) makes Fortran code difficult to read
and extend as well as accessing other 3rd party libraries.

To handle multiple dimension in C++, we make use of CEAN
(C/C++ Extension for Array Notation) as provided by Cilk+.
Additionally, Cilk+ allows to use C99 complex numbers within
C++, which facilitates vectorization. Although the code
is non-standard, the authors believe that Cilk+ will establish
as a new-standard such as OpenMP did (e.g. array sections in Section 2.6 OpenMP 4.0rc2).

License 
-------------------------

The code is license under the [GNU Public License Version 3](http://www.gnu.org/licenses/gpl.html)
(or at your opinion any later version). 
 
 
