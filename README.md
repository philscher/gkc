The GKC++ code
==============================

Description
-----------

GKC++ (GyroKinetic Code++) is a delta-f gyro-kinetic 
code. Currently features are 
      
  * Various geometries (sheared, shearless, toroidal)
  * multiple species support 
  * fully electromagnetic code 
  * Interfaces for various field and Vlasov solvers
  * Preliminary support for global simulations
  * 3-stage parallelization (Vectorization, OpenMP and MPI) 
 
However, this code is still considered experimental
and benchmark to other codes is still pending.

GKC++ is programmed in C++-11/Cilk+ which on greatly simplifies
readability and handling of multi-dimensional arrays.
However note that only a limited subset of compilers support
the Cilk+ extensions such as

* Intel Compiler (12.1+)
* GCC compiler (gcc-cilk plus branch)
* LLVM compiler (cilkplus branch) 

More information on Cilk+ can be found in [Cilk+ Homepage](http://cilkplus.org/)

Notes on data output 
----------------------
   
Data output is implemented using the  [HDF-5 library](www.hdfgroup.org/HDF5). 
     
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
as a new-standard such as OpenMP did.

License 
-------------------------

The code is license under the  
[GNU Public License Version 3](http://www.gnu.org/licenses/gpl.html)
(or at your opinion any later version). 
 
 
