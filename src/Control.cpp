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

#include <fenv.h>
#include <csignal>

bool force_exit = false;

Parallel *gl_parallel;

Control::Control(Setup *setup, Parallel *_parallel, Analysis *_analysis) : parallel(_parallel), analysis(_analysis) {

  gl_parallel = parallel;
        // set our control file, this is same for all processes and set to mpi root process id 
   if (setup->get("Control.useControlFile", 0)) {
        cntrl_file_name   = "gkc_" + Setup::number2string(parallel->master_process_id) + ".stop";
   } else cntrl_file_name = "";
        maxKineticEnergy  = setup->get("Control.MaxKineticEnergy", 1.e355);
        maxElectricEnergy = setup->get("Control.MaxElectricEnergy", 1.e355);
        maxMagneticEnergy = setup->get("Control.MaxMagneticEnergy", 1.e355);
        


        // Some Tests to check the function
        maxRunningTime = Setup::getSecondsFromTimeString(setup->get("Control.MaxRunningTime", "0s"));

        // set start Time
        startTime = time (NULL);



     setSignalHandler();
};
   
void Control::signalForceExit(bool val) { force_exit = val; }


Control::~Control() {
        // clean up stop file
        if(cntrl_file_name != "") remove(cntrl_file_name.c_str());
 }
     
    void printLoopStopReason();

int control_triggered_signal=0;

void signal_handler(int sig)
   {
    switch(sig) {
      case(SIGFPE)  :  std::cerr << "Floating point exception occured. Exiting" << std::endl;
                       // we have to unmask the signal otherwie program will slow down
                       control_triggered_signal |= SIGFPE;
                       //delete helios;
                       abort();
                       //signal(SIGFPE, SIG_IGN);
                       // now we raise SIGUSR1 which is propagetaed by mpirun to other processes
                       //raise(SIGUSR2);
                       break;

                       // when SIGINT or SIGTERM appears, e.g. openmpi first propagates SIGTERM too all procecess, waits a couple
                       // of seconds and then exists, (we should take care that this time is enough to finish the job)
      case(SIGINT)   : std::cerr << "SIGINT received, raising SIGTERM" << std::endl;
                       control_triggered_signal |= SIGINT;
                 //      signal(SIGINT, SIG_IGN);
                 //      raise(SIGTERM); 
                       break;
      case(SIGTERM)  : //helios->runningException(GKC_EXIT); 
                       std::cout << "SIGTERM" << std::endl;
                       control_triggered_signal |= SIGTERM;
                       break;
      case(SIGUSR1) :  gl_parallel->print("SIGUSR1 received");
                       control_triggered_signal |= SIGUSR1;
                       break;
      case(SIGUSR2) :  gl_parallel->print("SIGUSR2 received");
                       control_triggered_signal |= SIGUSR2;
                       break;
      default       :  std::cerr << "Unkown signal .... Ignoring\n";
    }


    if(force_exit == true) check(-1, DMESG("Signal received and \"force exit\" set"));
}
   
#ifdef OS_DARWIN
#include <Accelerate/Accelerate.h>
#include <xmmintrin.h>
#endif


void checkSignal() {



};


void Control::setSignalHandler() {


// Set Floating Point Capturing if in Debug mode

//#ifdef OS_LINUX
    feenableexcept(FE_DIVBYZERO | FE_INVALID);
//    feenableexcept(FE_ALL_EXCEPT);
    // Set the signal handler:
    signal(SIGFPE,   signal_handler);
    signal(SIGINT ,  signal_handler);
    signal(SIGUSR1 , signal_handler);
    signal(SIGUSR2 , signal_handler);
//#endif

#ifdef OS_DARWIN 
       _mm_setcsr( _MM_MASK_MASK &~  (_MM_MASK_OVERFLOW | _MM_MASK_INVALID | _MM_MASK_DIV_ZERO) );
#endif



};



bool Control::checkOK(Timing timing, Timing maxTiming) {
        
      cntrl.check(timing <= maxTiming, "(1) : Time Limit for simulation reached");
      if(cntrl_file_name != "") cntrl.check(ifstream(cntrl_file_name.c_str()) == NULL, "(1) : Manual stop bu using file.stop trigger");
      
      
      cntrl.check(((time(NULL)-startTime) < maxRunningTime) || (maxRunningTime == 0), "(1) : Running Time Limit Reached");
     
   
      // check for some physical limits exceeded
      double phiEnergy, ApEnergy, BpEnergy;
      analysis->getFieldEnergy(phiEnergy, ApEnergy, BpEnergy);
      
      cntrl.check( phiEnergy <= maxElectricEnergy, "(2) Electric Energy is over Maximum"     );
      //cntrl.check( ApEnergy  <= maxApEnergy      , "(2) Magnetic (Ap) Energy is over Maximum");
      //cntrl.check( BpEnergy  <= maxMagneticEnergy, "(2) Magnetic (Bp) Energy is over Maximum");
      
      // check Signals
      cntrl.check(!(control_triggered_signal & SIGINT ), "(3) Interrupted by SIGINT" );
      cntrl.check(!(control_triggered_signal & SIGTERM), "(3) Interrupted by SIGTERM");
      cntrl.check(!(control_triggered_signal & SIGUSR1), "(3) Interrupted by SIGUSR1");
      cntrl.check(!(control_triggered_signal & SIGUSR2), "(3) Interrupted by SIGUSR2");
      cntrl.check(!(control_triggered_signal & SIGFPE ), "(3) Interrupted by SIGUSR2");


      // FIXME no bool datatype available 
      cntrl.check(parallel->collect((int) cntrl.isOK(), OP_BAND) > 0, "(4) Interupted by other processor"); 

      return cntrl.isOK();
    }
    
void Control::printLoopStopReason() {
        parallel->print(std::string("\nMain Loop finished due to ") + cntrl.getMessage());
    };
    
void Control::runningException(int status, char *error_message) {
        // catch secondary exceptions
        std::cerr << error_message << std::endl;
    //    if     (status == GKC_FINISH) abort_run = 1;
    //    else if(status == GKC_EXIT  ) delete fileIO;
   //     else    check(-1, DMESG("No such status"));
    #ifdef GKC_PARALLEL_MPI
        parallel->barrier();
    #endif
     
    };
        

void Control::printOn(ostream &output) const 
{
            output << "Control    | phi^2 " << (maxElectricEnergy > 0. ? Setup::number2string(maxElectricEnergy) : "off") << std::endl;
}
