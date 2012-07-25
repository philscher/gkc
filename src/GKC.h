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
#include "Visualization.h"
#include "Event.h"
#include "Eigenvalue.h"
#include "TestParticle.h"
#include "TimeIntegration.h"
#include "Transport.h"
#include "Geometry.h"
#include "Benchmark.h"

#include "config.h"

#include <string>

/**
 *    @mainpage
 *    @brief The main GKC class which initializes sub-modules
 *         
 * 
 *    GKC (GyroKinetic Code) is a local, delta-f gyro-kinetic 
 *    code. Currently support is restricted to
 *     
 *    * Various geometries (sheared, shearless, toroidal [untested])
 *    * multiple species support 
 *    * fully electro-magnetic code (warning untested)  
 *    * Interfaces for various field and Vlasov solvers
 *    * Preliminary support for global simulations
 *    * hybrid parallelization (MPI and OpenMP) 
 *
 *    Although it is not feature-rich as other gyro-kinetic codes,
 *    it provides an object-oriented approach using C++. Thus
 *    the emphasis is on simplicity and readability using abstraction.
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
 *       Future work : Show parallelization rate on GKC for CBC
 *        
 *       Vectorization support is assured by using proper phragmas
 *     
 *  
 *    Notes on data output :
 *  
 *       (Parallel) data output is implemented using the HDF-5 library
 *       (<www.hdfgroup.org/HDF5>). 
 *    
 *    Notes on Speed :
 * 
 *       Fortran is popular in numerical science communities due to
 *       it's power of handling multi-dimensional arrays, as well the
 *       many mathematical functions it supports.
 * 
 *       However, the lack of templates, function overloading, pointers
 *       and classes (although some of these features are now supported
 *       by the Fortran-03 standard) makes the code difficult to read.
 *       Additionally for most function calls speed is not of a concern
 *       e.g. in the initialization phase - or consist of a library call
 *       to MPI, FFT routine, PETSc, HDF-5 etc. Using Profiling for a
 *       benchmark case, we concluded that about 80% is spend in the Vlasov
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
 *       any later version). You are free to use the code for your research,
 *       referencing the code, citing the code or making the authors of
 *       this code co-authors is NOT required. Although depending on your
 *       choice highly appreciated, and will benefit further development.
 *
 *    @brief  The main class handles all initializations and
 *            time step iterations
 *
 */
class GKC
{
private:
    FileIO        *fileIO;             ///< Used for Data Input/Output
    Vlasov        *vlasov;             ///< Vlasov equation solver
    Grid          *grid;               ///< Grid initialization and boundaries
    Setup        *setup;               ///< Reads configuration files
    Parallel      *parallel;           ///< Parallel communication functions
    Analysis      *analysis;           ///< Data Analysis
    Fields        *fields;             ///< Source calculation and field solvers
    Control       *control;            ///< Program flow control
    FFTSolver     *fftsolver;          ///< FFTSolver used
    TestParticles  *particles;         ///< Test particles
    Eigenvalue     *eigenvalue;        ///< Eigenvalue calculations
    Geometry *geometry;                ///< Geometry module
    Event         *event;              ///< Programmable events 
    Init          *init;               ///< Initialization for plasma
    Visualization *visual;             ///< Visualization
    Benchmark     *bench;              ///< Interface for profiling and benchmarking
    TimeIntegration *timeIntegration;  ///< Numerical time integration

    /**
    * @brief Run the code, as "IVP" (initial value code)
             or "Eigenvalue" code.
    */
    std::string Helios_Type; 
  

    void printSettings();

public:
    /**
    * @param setup configuration parameters
    **/
    GKC(Setup *setup);
   ~GKC();  
   
   int mainLoop();
};


#endif // __GKC_H_
