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

#include "Vlasov/Vlasov.h"
#include "Setup.h"
#include "FileIO.h"
#include "Fields/Fields.h"
#include "Init.h"
#include "Parallel/Parallel.h"
#include "Collisions/Collisions.h"
#include "Collisions/LenardBernstein.h"
#include "FFTSolver/FFTSolver.h"
#include "Control.h"
#include "Timing.h"
#include "Visualization/Visualization.h"
#include "Eigenvalue/Eigenvalue.h"
#include "TimeIntegration/TimeIntegration.h"
#include "TimeIntegration/ScanLinearModes.h"
#include "Geometry.h"
#include "Collisions/Collisions.h"
#include "Analysis/Diagnostics.h"
#include "Analysis/Benchmark.h"
#include "Analysis/Event.h"
#include "Analysis/TestParticle.h"

#include "config.h"

#include <string>

/**
 *    @brief The main GKC class which initializes sub-modules
 *         
 * 
 *    The main class handles all initializations and
 *    time step iterations
 *
 */
class GKC
{
private:
    FileIO          *fileIO;             ///< Used for Data Input/Output
    Vlasov          *vlasov;             ///< Vlasov equation solver
    Collisions      *collisions;         ///< Collisional operator
    Grid            *grid;               ///< Grid initialization and boundaries
    Setup           *setup;              ///< Reads configuration files
    Parallel        *parallel;           ///< Parallel communication functions
    Diagnostics     *diagnostics;        ///< Data Diagnostic
    Fields          *fields;             ///< Source calculation and field solvers
    Control         *control;            ///< Program flow control
    FFTSolver       *fftsolver;          ///< FFTSolver used
    TestParticles   *particles;          ///< Test particles
    Eigenvalue      *eigenvalue;         ///< Eigenvalue calculations
    Geometry        *geometry;           ///< Geometry module
    Event           *event;              ///< Programmable events 
    Init            *init;               ///< Initialization for plasma
    Visualization   *visual;             ///< Visualization
    Benchmark       *bench;              ///< Interface for profiling and benchmarking
    TimeIntegration *timeIntegration;    ///< Numerical time integration
    ScanLinearModes *scanModes;          ///< Scan over linear modes

    /**
    * @brief Run the code, as "IVP" (initial value code)
             or "Eigenvalue" code.
    */
    std::string gkc_SolType; 
  

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
