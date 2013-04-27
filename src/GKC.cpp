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

#include "config.h"

#include "GKC.h"
#include "Setup.h"
#include "Plasma.h"

#include "FieldsFFT.h"
#include "FieldsHermite.h"

#include "FFTSolver/FFTSolver_fftw3.h"

#include "TimeIntegration_PETSc.h"

#include "Collisions/LenardBernstein.h"
#include "Collisions/PitchAngle.h"

#include "Vlasov/Vlasov_Cilk.h"
#include "Vlasov/Vlasov_Aux.h"
#include "Vlasov/Vlasov_Optim.h"
#include "Vlasov/Vlasov_Island.h"

#include "Geometry/GeometrySA.h"
#include "Geometry/Geometry2D.h"
#include "Geometry/GeometryShear.h"
#include "Geometry/GeometrySlab.h"
#include "Geometry/GeometryCHEASE.h"

#include "Tools/System.h"

#include "Eigenvalue/Eigenvalue_SLEPc.h"
#include "Visualization/Visualization_Data.h"


GKC::GKC(Setup *_setup) : setup(_setup)  
{

  // Read Setup 
  std::string fft_solver_name = setup->get("GKC.FFTSolver"       , "fftw3");
  std::string psolver_type    = setup->get("Fields.Solver"       , "DFT"  );
  std::string vlasov_type     = setup->get("Vlasov.Solver"       , "Aux"  );
  std::string timeInt_type    = setup->get("GKC.TimeIntegration" , "Explicit"  );
  std::string collision_type  = setup->get("Collisions.Solver"   , "None" );
  gkc_SolType                 = setup->get("GKC.Type"            , "IVP"  );
  std::string geometry_Type   = setup->get("GKC.Geometry"        , "2D"   );


  ///////////////////////////////   Load subsystems   /////////////////////////////

  parallel  = new Parallel(setup);

  parallel->print("Initializing GKC++\n");

  fileIO    = new FileIO(parallel, setup);
  grid      = new Grid(setup, parallel, fileIO);

  // ugly here, however in parallel constructor fileIO is not defined yet
  parallel->initData(setup, fileIO);
  
  bench     = new Benchmark(setup, parallel, fileIO); 

  if     (geometry_Type == "SA"   ) geometry  = new GeometrySA(setup, grid, fileIO);
  else if(geometry_Type == "2D"   ) geometry  = new Geometry2D(setup, grid, fileIO);
  else if(geometry_Type == "Slab" ) geometry  = new GeometrySlab(setup, grid, fileIO);
  else if(geometry_Type == "Shear") geometry  = new GeometryShear(setup, grid, fileIO);
  else check(-1, DMESG("No such Geometry"));

  plasma    = new Plasma(setup, fileIO, geometry);

  // Load fft-solver 
  if(fft_solver_name == "") check(-1, DMESG("No FFT Solver Name given"));
#ifdef FFTW3
  else if(fft_solver_name == "fftw3") fftsolver = new FFTSolver_fftw3(setup, parallel, geometry);
#endif
  else check(-1, DMESG("No such FFTSolver name"));
 
  // Load field solver
  if     (psolver_type == "DFT"    ) fields   = new FieldsFFT(setup, grid, parallel, fileIO, geometry, fftsolver);
#ifdef GKC_HYPRE
  else if(psolver_type == "Hypre"  ) fields   = new FieldsHypre(setup, grid, parallel, fileIO,geometry, fftsolver);
#endif
  else if(psolver_type == "Hermite") fields   = new FieldsHermite(setup, grid, parallel, fileIO,geometry);
  else    check(-1, DMESG("No such Fields Solver"));
  
  // Load Collisonal Operator
  if     (collision_type == "None"      ) collisions = new Collisions                (grid, parallel, setup, fileIO, geometry); 
  else if(collision_type == "LB"        ) collisions = new Collisions_LenardBernstein(grid, parallel, setup, fileIO, geometry); 
  else if(collision_type == "PitchAngle") collisions = new Collisions_PitchAngle     (grid, parallel, setup, fileIO, geometry); 
  else    check(-1, DMESG("No such Collisions Solver"));
    
  // Load Vlasov Solver
  if(vlasov_type == "None" ) check(-1, DMESG("No Vlasov Solver Selected"));
  else if(vlasov_type == "Cilk"   ) vlasov  = new VlasovCilk  (grid, parallel, setup, fileIO, geometry, fftsolver, bench, collisions);
  else if(vlasov_type == "Aux"    ) vlasov  = new VlasovAux   (grid, parallel, setup, fileIO, geometry, fftsolver, bench, collisions);
  else if(vlasov_type == "Island" ) vlasov  = new VlasovIsland(grid, parallel, setup, fileIO, geometry, fftsolver, bench, collisions);
  else if(vlasov_type == "Optim"  ) vlasov  = new VlasovOptim (grid, parallel, setup, fileIO, geometry, fftsolver, bench, collisions);
  else   check(-1, DMESG("No such Fields Solver"));
  
  // Load some other general modules
  diagnostics     = new Diagnostics(parallel, vlasov, fields, grid, setup, fftsolver, fileIO, geometry); 
  visual          = new Visualization_Data(grid, parallel, setup, fileIO, vlasov, fields);
  event           = new Event(setup, grid, parallel, fileIO, geometry);
  eigenvalue      = new Eigenvalue_SLEPc(fileIO, setup, grid, parallel); 
  init            = new Init(parallel, grid, setup, fileIO, vlasov, fields, geometry);
  control         = new Control(setup, parallel, diagnostics);
  particles       = new TestParticles(fileIO, setup, parallel);
 
  if     (timeInt_type == "Explicit") timeIntegration = new TimeIntegration      (setup, grid, parallel, vlasov, fields, particles, eigenvalue, bench);
  else if(timeInt_type == "Implicit") timeIntegration = new TimeIntegration_PETSc(setup, grid, parallel, vlasov, fields, particles, eigenvalue, bench);
  else   check(-1, DMESG("No such TimeIntegratioScheme in GKC.TimeIntegration"));
  
  //scanModes       = new ScanLinearModes  (setup, grid, parallel, vlasov, fields, fileIO);
  //scanEigen       = new ScanPoloidalEigen(setup, grid, parallel, vlasov, fields, fileIO, control, visual, diagnostics, timeIntegration);
  // Optimize values to speed up computation
  bench->bench(vlasov, fields);
  
  printSettings();   
  setup->check_config();

}
          

int GKC::mainLoop()   
{

  parallel->print("Running main loop");

  //int bench_id = benchmark->start("MainLoop", Benchmark::time)

  unsigned int start_time = System::getTime();
   
  if (gkc_SolType == "IVP") {
 
    Timing timing(0,0.) ;
      
    timeIntegration->setMaxLinearTimeStep(eigenvalue, vlasov, fields);

    bool isOK = true; 
         
    ////////////////////////  Starting OpenMP global threads   /////////////////////////
    //
    //  OpenMP threads exists throughout the main iteration loop.
    //  Synchronization occurs e.g. during implicit barriers (e.g. #pragma omp for)
    //
    //  We need to take care as every variable which is defined prior to this point
    //  (e.g. classes, arrays) are shared variables which may lead to race conditions or
    //  false sharing if not used properly. 
    //  After this point, stack variables are private, thus if a domain is decomposed
    //  reductions clauses are required.
    //
    //  @todo Inside parallel region we can decompose NkylD, NkyLuD, when set to private 
    //        in omp parallel. Test if more efficient.
    //
    //#pragma omp parallel  private(NkyLlD, NkyLuD)
    #pragma omp parallel private(threadID) 
    {
      parallel->setThreadID();
   
      do {
        bench->start("A", 1);        
        // integrate for one time-step, give current dt as output
        const double dt = timeIntegration->solveTimeStep(vlasov, fields, particles, timing);     
        bench->stop("A", 1);
        // Analysis results and output data (currently singlethreaded)
        #pragma omp master
        {
          
          vlasov->writeData(timing, dt);
          fields->writeData(timing, dt);
          visual->writeData(timing, dt);

          diagnostics->writeData(timing, dt);
          event->checkEvent(timing, vlasov, fields);

          // flush in regular intervals in order to minimize HDF-5 file 
          // corruption in case of an abnormal program termination
          fileIO->flush(timing, dt);  
             
          isOK = control->checkOK(timing, timeIntegration->maxTiming);

        }
        #pragma omp barrier
    
      } while(isOK);

    } // parallel section
   
    control->printLoopStopReason();

  }  
  else if(gkc_SolType == "Eigenvalue") {
   
    eigenvalue->solve(vlasov, fields, visual, control);


  } 
  else if(gkc_SolType == "Scan") {
   // scanModes->solve(vlasov, fields, timeIntegration, eigenvalue, init, visual); 
  }
  else if(gkc_SolType == "Eigen") {
    //scanEigen->solve(vlasov, fields, timeIntegration, eigenvalue, init, visual); 

  }
  else  check(-1, DMESG("No Such gkc.Type Solver"));
  
  bench->save("MainLoopTime", System::getTime() - start_time);

  parallel->print("Simulation finished normally ... ");

  return 0;
}


GKC::~GKC()
{
   
  // Shutdown submodules : ORDER IS IMPORTANT !!
  delete fields;
  delete visual;
  delete vlasov;
  delete collisions; 
  delete geometry;
  delete grid;
  delete diagnostics;
  delete particles;
  delete eigenvalue;
  delete event;
  //delete scanModes;
  delete bench;
//  delete scanEigen;
  delete fileIO; // once this is successful, file cannot
                 // be corrupted anymore.
  
  delete plasma;
  delete fftsolver;
  delete control;
  delete parallel;
  delete init;
  delete timeIntegration;

}


void GKC::printSettings() 
{
  
  std::stringstream infoStream;
  
  time_t start_time = std::time(0); 
  
  infoStream 
  
    << "Welcome to " << PACKAGE_NAME << " (" << PACKAGE_VERSION <<")  " << PACKAGE_BUGREPORT <<  "      Date :  " << std::ctime(&start_time)         
    << "-------------------------------------------------------------------------------" << std::endl
    << *grid << *plasma << *fileIO << *setup << *vlasov  << *fields << *geometry << *init << *parallel << *fftsolver << *timeIntegration << *collisions;
    
  infoStream << *control << *bench << std::endl;
  infoStream << "-------------------------------------------------------------------------------" << std::endl << std::endl << std::flush;

  parallel->print(infoStream.str());

}
    
