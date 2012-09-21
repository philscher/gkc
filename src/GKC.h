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
#include "Eigenvalue.h"
#include "TimeIntegration.h"
#include "Geometry.h"
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
