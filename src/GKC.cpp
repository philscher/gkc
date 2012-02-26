/*
 * =====================================================================================
 *
 *       Filename: GKC.cpp
 *
 *    Description: Main file
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include "GKC.h"
#include "FEMSolver.h"

#include "Setup.h"
#include "config.h"
#include<ctime>

#include "Plasma.h"
#include "FieldsFFT.h"
#include "FieldsFEM.h"
#include "FieldsHypre.h"

#include "TimeIntegration.h"
#include "TimeIntegration_PETSc.h"

//#include "Vlasov/Vlasov_Fortran.h"
//#include "Vlasov/Vlasov_Blitz.h"

//#ifdef HELIOS_CILK
#include "Vlasov/Vlasov_Cilk.h"
//#endif
#include "Eigenvalue/Eigenvalue_SLEPc.h"

Plasma *plasma;

#include "FFT/FFTSolver_fftw3.h"

#include <unistd.h>

GeneralArrayStorage<6> HeliosStorage;
GeneralArrayStorage<4> HeliosStorage4;

Helios::Helios(Setup *_setup) : setup(_setup)  {

    // define Storage for helios (used by blitz++ arrays)
    HeliosStorage.ordering() = fourthDim, firstDim, secondDim,  thirdDim,  fifthDim, sixthDim; 
    HeliosStorage.base() = NxLlB, NkyLlD, NzLlB, NvLlB, NmLlB, NsLlB;
    
    // define Storage for helios (used by blitz++ arrays)
    HeliosStorage4.ordering() = fourthDim, firstDim, secondDim,  thirdDim; 

	// *************** Load subsystems ************** //
    parallel  = new Parallel(setup);
    fileIO    = new FileIO(parallel, setup);
    grid      = new Grid(setup, parallel, fileIO);
    geometry  = new HELIOS_GEOMETRY(setup, fileIO);
    plasma    = new Plasma(setup, fileIO, geometry);


    std::string fft_solver_name = setup->get("Helios.FFTSolver", "fftw3");
    
    if(fft_solver_name == "") check(-1, DMESG("No FFT Solver Name given"));
#ifdef FFTW2MPI
//    else if(fft_solver_name == "fftw2_mpi")  fftsolver = new FFTSolver_fftw2_mpi    (setup, parallel, geometry);
#endif
#ifdef FFTW3
   else if(fft_solver_name == "fftw3") fftsolver = new FFTSolver_fftw3    (setup, parallel, geometry);
#endif
   else check(-1, DMESG("No such FFTSolver name"));
 
//    femsolver  = new FEMSolver(setup, parallel);
  
  
    std::string vlasov_type = setup->get("Vlasov.Solver", "Cilk");
   if(vlasov_type == "None" ) check(-1, DMESG("No Vlasov Solver Selected"));
//    if(vlasov_type == "Blitz"  ) vlasov     = new VlasovBlitz  (grid, parallel, setup, fileIO, geometry, fftsolver);
//#ifdef HELIOS_CILK
    else if(vlasov_type == "Cilk"  ) vlasov     = new VlasovCilk  (grid, parallel, setup, fileIO, geometry, fftsolver);
//#endif
    else   check(-1, DMESG("No such Fields Solver"));
	
    std::string psolver_type = setup->get("Helios.FieldsSolver", "DFT");
    if     (psolver_type == "DFT"  ) fields   = new FieldsFFT(setup, grid, parallel, fileIO, geometry, fftsolver);
#ifdef LIBMESH
    else if(psolver_type == "FEM"  ) fields   = new FieldsFEM(setup, grid, parallel, geometry, femsolver);
#endif
#ifdef HELIOS_HYPRE
    else if(psolver_type == "Hypre") fields   = new FieldsHypre(setup, grid, parallel, fileIO,geometry, fftsolver);
#endif
    else    check(-1, DMESG("No such Fields Solver"));
 
    analysis = new Analysis(parallel, vlasov, fields, grid, setup, fftsolver, fileIO, geometry); 
    visual   = new Visualization(grid, parallel, setup, fileIO, vlasov, fields);
    event    = new Event(setup, grid, parallel, fileIO, geometry);
    
    eigenvalue = new Eigenvalue_SLEPc(setup, grid, parallel); 
    // *******************************************   //

   
    // call Initital conditions
    init      = new  Init(grid, setup, vlasov, fields, geometry);

    control  = new Control(setup, parallel, analysis);
    bench    = new Benchmark(setup);
 
    particles = new TestParticles(fileIO, setup, parallel);
    Helios_Type = setup->get("Helios.Type", "IVP");
   

    timeIntegration = new TimeIntegration_PETSc(setup, grid, parallel, vlasov, fields, eigenvalue);

    printSettings();	
    setup->check_config();


    }
	       

int Helios::mainLoop()   {

   if (Helios_Type == "IVP") {
 
            Timing timing;
  
            for(; control->checkOK(timing, timeIntegration->maxTiming);){
        
         	    event->checkEvent(timing, vlasov, fields);
         	    analysis->writeData(timing, timeIntegration->dt);
         	    fields->writeData(timing, timeIntegration->dt);
         	    visual->writeData(timing, timeIntegration->dt);
        	
        	    timeIntegration->solveTimeStep(vlasov, fields, particles, timing);        

                // OK Time Step finished
       	        bench->benchmark(vlasov, fields, timing);
 		    }

   
        control->printLoopStopReason();

   }  
   else if(Helios_Type == "Eigenvalue") {
	eigenvalue->solve(vlasov, fields, visual, control);

   } else  check(-1, DMESG("No Such Helios.Type Solver"));


   return HELIOS_SUCCESS;
}


Helios::~Helios(){
   // Shutdown submodules : ORDER IS IMPORTANT !!
        delete fields;
        delete visual;
        delete vlasov;
        delete geometry;
     //   delete femsolver;
        delete grid;
        delete analysis;
        delete particles;
        delete fileIO;
        delete island;
        // no need to delete table ?! Anyway SLEPc/PETSc seems to crash
        delete eigenvalue;
    	parallel->barrier();       
        // some times FFT crashes, be sude that it is below fileIO (hdf-5 is closed)
        delete fftsolver;
        delete control;
        delete parallel;
        delete init;
        delete bench;
        delete timeIntegration;
}



void Helios::printSettings() {
    std::stringstream infoStream;
  
    time_t start_time = std::time(0); 
  infoStream 
            << PACKAGE_NAME << "(version " << PACKAGE_VERSION <<") : " << std::ctime(&start_time)         
            << "-------------------------------------------------------------------------------" << std::endl
            << *grid << *plasma << *fileIO << *setup << *vlasov  << *fields << *geometry << *init << *parallel << *fftsolver;
            infoStream << "Timing     |  " << timeIntegration;
            infoStream << "Control    |  File : " << control->cntrl_file_name << std::endl;
            infoStream << "-------------------------------------------------------------------------------" << std::endl << std::endl << std::flush;


   parallel->print(infoStream);

}
    
