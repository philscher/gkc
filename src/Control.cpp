/*
 * =====================================================================================
 *
 *       Filename: Control.cpp
 *
 *    Description: Handles signals and other triggers
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include "Control.h"
#include "GKC.h"
#include "Tools/System.h"

#include <fenv.h>
#include <csignal>
#include <string.h>

bool force_exit = false;
int control_triggered_signal=0;

Parallel *gl_parallel; // need extern variables for signal handler


/**
*
*   @brief Signal handler for gkc catched signals
*
*   Is the usage of std::cerr allowed in the signal handler? 
*   (in C printf is not !)
*
**/ 
void signal_handler(int sig)
{
  
  switch(sig) {

    case(SIGFPE ) :
                std::cerr << "SIGFPE received. Exiting ..." << std::endl; 
                // ignore subsequent signals (otherwise program may slow down) 
                signal(sig, SIG_IGN);
                System::printStackTrace();
                control_triggered_signal |= sig;
    break;
    case(SIGSEGV) : 
                std::cerr << "SIGSEGV received. Exiting ..." << std::endl; 
                // ignore subsequent signals (otherwise program may slow down) 
                signal(sig, SIG_IGN);
                System::printStackTrace();
    abort();
                break;
    default       :
                std::cerr << strsignal(sig) << " received. Exiting ..." << std::endl; 
                control_triggered_signal |= sig;
      
 
    }
  
    if(force_exit == true) check(-1, DMESG("Signal received and \"force exit\" set"));
}

 

///////////// Control class member function definitions //////////////////


Control::Control(Setup *setup, Parallel *_parallel, Diagnostics *_diagnostics) 
: parallel(_parallel), diagnostics(_diagnostics) 
{
  
  gl_parallel = parallel;
   
  // distributed unified id (the process id of MPI master process) to all processes
  // useAControlFile to stop process 
  // [ method used by the GKW code to controll MPI processes ]
  if (setup->get("Control.useControlFile", 0)) {
        
        int master_process_id = parallel->reduce( (parallel->myRank == 0) ? System::getProcessID() : 0, Op::sum, DIR_ALL);
        cntrl_file_name   = "gkc_" + Setup::num2str(master_process_id) + ".stop";
   } 
   else cntrl_file_name = "";

   maxKineticEnergy  = setup->get("Control.MaxKineticEnergy", 1.e355);
   maxElectricEnergy = setup->get("Control.MaxElectricEnergy", 1.e355);
   maxMagneticEnergy = setup->get("Control.MaxMagneticEnergy", 1.e355);
        
   // Some Tests to check the function
   maxRunningTime = Setup::getSecondsFromTimeString(setup->get("Control.MaxRunningTime", "0s"));
   // set start Time
   startTime = time(NULL);
  
   // set signal handler to catch SIGUSR1, SIGFPE and SIGTRP
   setSignalHandler();

};
   
Control::~Control() 
{
        // clean up stop file if control file is used
        if(cntrl_file_name != "") remove(cntrl_file_name.c_str());
}


void Control::signalForceExit(bool val) 
{ 
  force_exit = val; 
}
     

void Control::setSignalHandler() 
{


     // Set Floating Point Capturing if in Debug mode
    feenableexcept(FE_DIVBYZERO | FE_INVALID);
//    feenableexcept(FE_ALL_EXCEPT);
    // Set the signal handler:
    signal(SIGFPE  , signal_handler);
    signal(SIGINT  , signal_handler);
    signal(SIGUSR1 , signal_handler);
    signal(SIGUSR2 , signal_handler);
    signal(SIGSEGV , signal_handler);
    signal(SIGBUS  , signal_handler);
// catch all signals (SIGKILL and SIGSTOP are uncatchable)
    for ( int i = 0; i < 32; i++ ) signal( i, signal_handler );
//    for ( int i = 0; i < 32; i++ ) signal( 0x1111111111111111, signal_handler );


//#endif



};


bool Control::checkOK(Timing timing, Timing maxTiming) 

{
   int isOK = true;
  // #pragma omp single 
  // {

      cntrl.check(timing <= maxTiming, "(1) : Time Limit for simulation reached");
      if(cntrl_file_name != "") cntrl.check(std::ifstream(cntrl_file_name.c_str()) == NULL, "(1) : Manual stop bu using file.stop trigger");
      
      
      cntrl.check(((time(NULL)-startTime) < maxRunningTime) || (maxRunningTime == 0), "(1) : Running Time Limit Reached");
  
      // Problems with OpenMP, due to use of FFT ?
      // check for some physical limits exceeded
//      double phiEnergy, ApEnergy, BpEnergy;
//      analysis->getFieldEnergy(phiEnergy, ApEnergy, BpEnergy);
//      
//      cntrl.check( phiEnergy <= maxElectricEnergy, "(2) Electric Energy is over Maximum"     );
      //cntrl.check( ApEnergy  <= maxApEnergy      , "(2) Magnetic (Ap) Energy is over Maximum");
      //cntrl.check( BpEnergy  <= maxMagneticEnergy, "(2) Magnetic (Bp) Energy is over Maximum");
      
      // check Signals
      cntrl.check(!(control_triggered_signal & SIGINT ), "(3) Interrupted by SIGINT" );
      cntrl.check(!(control_triggered_signal & SIGTERM), "(3) Interrupted by SIGTERM");
      cntrl.check(!(control_triggered_signal & SIGUSR1), "(3) Interrupted by SIGUSR1");
      cntrl.check(!(control_triggered_signal & SIGUSR2), "(3) Interrupted by SIGUSR2");
      cntrl.check(!(control_triggered_signal & SIGFPE ), "(3) Interrupted by SIGFPE" );
      cntrl.check(!(control_triggered_signal & SIGCONT), "(3) Interrupted by SIGCONT" );


      isOK = cntrl.isOK(); 
      // Check if stopping signal was triggered on other processes
      //cntrl.check(parallel->collect(cntrl.isOK(), Op::LAND, DIR_ALL) == true, "(4) Interupted by other processor"); 
      cntrl.check(parallel->reduce(isOK, Op::land, DIR_ALL) == true, "(4) Interupted by other processor"); 
   //}
  // #pragma omp flush
      return isOK;

}
    
void Control::printLoopStopReason() 
{

        parallel->print(std::string("\nMain Loop finished due to ") + cntrl.getMessage());

};
    
void Control::runningException(int status, char *error_message) {
        // catch secondary exceptions
        std::cerr << error_message << std::endl;
    
        parallel->barrier();

};
        

void Control::printOn(std::ostream &output) const 
{
            output << "Control    | phi^2 " << (maxElectricEnergy > 0. ? Setup::num2str(maxElectricEnergy) : "off") << std::endl;
}
