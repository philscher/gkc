/*
 * =====================================================================================
 *
 *       Filename: GKC.h
 *
 *    Description: Main file
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __GKC_H_
#define __GKC_H_

#include "config.h"

#include <string>

#include "Global.h"

#include "Vlasov.h"
#include "Setup.h"
#include "FileIO.h"
#include "Fields.h"
#include "Init.h"
#include "Parallel.h"
#include "Collisions.h"
#include "FFTSolver.h"
#include "Control.h"
#include "Timing.h"
#include "Benchmark.h"
#include "Visualization.h"
#include "Event.h"
#include "Eigenvalue.h"
#include "TestParticle.h"
#include "TimeIntegration.h"
#include "Transport.h"
#include "Geometry.h"
#include "GeometryShear.h"
#include "GeometrySlab.h"
#include "Geometry2D.h"

/**
 *    @mainpage
 *
 *         
 * 
 *    GKC (GyroKinetic Code) is a local, delta-f gyro-kinetic 
 *    code. Currently support is restricted to
 *     
 *    * Various geometrie (sheared, shearless, toroidal [untested])
 *    * multiple species support 
 *    * fully electro-magnetic code (warning untested)  
 *    * Interfaces for various field and vlasov solvers
 *    * Prelimenary support for global simulations
 *    * hybrid parallelization (MPI and OpenMP) 
 *
 *    Alltough it is not featurich as other gyro-kinetic codes,
 *    it provides an object-oriented approach using C++. Thus
 *    the emphasis is on simplicity and readbility using abstraction.
 *
 *    (A little technical notes)
 * 
 *    Notes on parallelization :
 *   
 *       The code uses OpenMP (for poloidal direction only) and MPI
 *       to achieve a high parallelization rate. The Parallel class
 *       provides an abstraction for most MPI function, thus direct
 *       calls to MPI functions are not necessary. Function 
 *       overloading is provided to achieve a common interface.
 *       
 *       Future work : Show parallelization rate on Helios for CBC
 *        
 *       Vectorization support is assured by using proper phragmas
 *     
 *  
 *    Notes on data output :
 *  
 *       (Parallel) data output is implemented using the HDF-5 library
 *       (www.hdfgroup.org/HDF5). 
 *    
 *    Notes on Speed :
 * 
 *       Fortran is popular in numerical science comunnities due to
 *       it's power of handling multi-dimensional arrays, as well the
 *       many mathematical functions it supports.
 * 
 *       However, the lack of templates, function overloading, pointers
 *       and classes (allhough some of these features are now supported
 *       by the Fortran-03 standard) makes the code difficulat to read.
 *       Additionally for most function calls speed is not of a concern
 *       e.g. in the intialization phase - or consist of a library call
 *       to MPI, FFT routine, PETSc, HDF-5 etc. Using Profiling for a
 *       benchmar case, we concluded that about 80% is spend in the Vlasov
 *       equation solver, and 10% in the underlying FFT solver, with only
 *       a minor fraction (around %5) inside other GKC classes.
 *       
 *       To handle multiple dimension in C++ this code makes use of 
 *       blitz++ (blitz.sourceforge.net), a meta-template library, which
 *       promised to allow speed comparable to Fortran. However, we
 *       observe a speed of only 40% using blitz-0.9. Thus for the Vlasov
 *       equation, implementation in Fortra-03 and Intel-Cilk Array Notation 
 *       is provided. Overall the speed compared to a native Fortran 
 *       implementation varies between (95-105%).
 *    
 *       
 *   Notes on License :
 *       The code is license under the GNU Public License Version 3 ( or
 *       any later version). You are free to use the code for your research
 *       referencing the code, citing the code or making the authors of
 *       this code co-authors is NOT required. Alltough depending on your
 *       choice highly appreciated, and will benefit further development.
 *
 *    @brief  The main class handles all initializations and
 *            time step iterations
 *
 */
class Helios
{
private:
    /** @brief class for data input */
    FileIO        *fileIO;
    Vlasov        *vlasov;
    Grid          *grid;
    Setup        *setup;
    Parallel      *parallel;
    Analysis      *analysis;
    Fields        *fields;   
    Control       *control;
    FFTSolver     *fftsolver;
    TestParticles  *particles;
    Eigenvalue     *eigenvalue;
    Geometry<HELIOS_GEOMETRY> *geometry;
    Event         *event;
    Init          *init;
    Benchmark *bench;
    Visualization *visual;
    TimeIntegration *timeIntegration;

    /**
    * @brief Run the code, as "IVP" (initial value code)
             or "Eigenvalue" code.
    */
    std::string Helios_Type; 
  

    void printSettings();

public:
    /**
    * @param setup configuration parameters
    */
    Helios(Setup *setup);
   ~Helios();  
   
   int mainLoop();
};


#endif // __HELIOS_H_
