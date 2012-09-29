@mainpage
  
 
GKC (GyroKinetic Code) is a local, delta-f gyro-kinetic 
code. Currently support is restricted to
      

* Various geometries (sheared, shearless, toroidal [untested])
* Multiple species support 
* fully electro-magnetic code (warning untested)  
* Interfaces for various field and Vlasov solvers
* Preliminary support for global simulations
* hybrid parallelization (MPI and OpenMP)
* Data output using [HDF-5](<www.hdfgroup.org/HDF5>)
* Data analysis script for python
 
Although it is not feature-rich as other gyro-kinetic codes,
it provides an object-oriented approach using C++. Thus
the emphasis is on simplicity and readability using abstraction.


Notes on Speed :
================
  
  For efficient numerical calculation, gkc++ makes use of CEAN.
  A extension to the C++ language currently supported by
  gcc-4.7 and Intel 12.1.R6 or newer. CEAN allows efficient 
  handling of multi-dimensional arrays, with speeds
  comparable to Fortran implementations.
        
License 
==================
   The code is license under the GNU Public License Version 3 ( or
   any later version). 
   

Other gyro-kinetic code
=========================

 * [gs2](http://gs2.sourceforge.net/)
 * [GENE](http://www.ipp.mpg.de/~fsj/gene/)
 * [gkw](http://www.gkw.org.uk/)
 * [GYRO](https://fusion.gat.com/theory/Gyro)


   

Profiling file: Profiling.md Parallelization file:Parallelization.md
 
