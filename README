
     GKC (GyroKinetic Code) is a local, delta-f gyro-kinetic 
     code. Currently features are 
      
     * Various geometries (sheared, shearless, toroidal)
     * multiple species support 
     * fully electromagnetic code (warning untested)  
     * Interfaces for various field and Vlasov solvers
     * Preliminary support for global simulations
     * 3-stage parallelization (Vectorization, OpenMP and MPI) 
 
     However, handle with care, as the code is still experimental
     and benchmark to other codes is still pending.


     Although it is not feature-rich as other gyro-kinetic codes,
     it provides an object-oriented approach using C++. Thus
     the emphasis is on simplicity and readability using abstraction.
 
     (A little technical notes)
  
     Notes on parallelization :
    
        The code uses OpenMP and MPI to achieve a high parallelization rate.
        The Parallel class provides an abstraction for most MPI function, 
        thus direct calls to MPI functions are not necessary. Function 
        overloading is provided to achieve a common interface.
        Current vectorization is optimized for AVX instruction set.
        
   
     Notes on data output :
   
        (Parallel) data output is implemented using the HDF-5 library
        (www.hdfgroup.org/HDF5). 
     
     Notes on Speed :
  
        Fortran is popular in numerical science communities due to
        it's power of handling multi-dimensional arrays, as well the
        many mathematical functions it supports.
  
        However, the lack of templates, function overloading, pointers
        and classes (although some of these features are now supported
        by the Fortran-03 standard) makes Fortran code difficult to read
        and extend.
        Additionally for most function calls speed is not of a concern
        e.g. in the initialization phase - or consist of a library call
        to MPI, FFT routine, PETSc, HDF-5 etc. Using Profiling,
        we concluded that about 80% of the computational spend in the
        Vlasov equation solver, and 10% in the underlying FFT solver,
        with only a minor fraction (around %5) inside other GKC classes.
       
        To handle multiple dimension in C++ this code makes use of
        use of CEAN (C/C++ Extension for Array Notation) supported
        by Intel C/C++ (12.1) and GCC (4.8) onwards as part
        of Cilk Plus. Additionally C99 complex numbers are used.
        This makes the code non C++ standard conform, and probably
        would not compile using other compilers and non x86 platforms.
        Hopefully, other compiler vendor will support it in the near
        future.

        CEAN greatly simplifies using multi-dimensional arrays in C++
        and a proper vectorization of array manipulations.

        
    Notes on License :

        The code is license under the GNU Public License Version 3 ( or
        at your opinion any later version). 
 
 
