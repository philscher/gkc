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


//#include "MagneticIsland.h"
/// Include necessary fftw3 headers



#include "Geometry.h"
#include "GeometryShear.h"
#include "GeometrySlab.h"
#include "Geometry2D.h"

/**
 *    Helios - The main class handles all initializations and
 *             time step iterations
 */

class Helios
{
private:
    FileIO        *fileIO;
    Vlasov        *vlasov;
    Grid          *grid;
    Setup        *setup;
    Parallel      *parallel;
    Analysis      *analysis;
    Fields        *fields;   
    Control       *control;
    FFTSolver     *fftsolver;
//    FEMSolver     *femsolver;
    TestParticles  *particles;
    Eigenvalue     *eigenvalue;
    Geometry<HELIOS_GEOMETRY> *geometry;
    Event         *event;
    Init          *init;
    Benchmark *bench;
    Visualization *visual;
    TimeIntegration *timeIntegration;
    // Plugin
//    MagneticIsland *island;


    std::string Helios_Type; 
  

    void printSettings();

public:
    Helios(Setup *_setup);
   ~Helios();  
   
   int mainLoop();
};


#endif // __HELIOS_H_
