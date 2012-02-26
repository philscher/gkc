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
#include "Helios.h"


extern Helios *helios;

int control_triggered_signal=0;

void signal_handler(int sig)
	{
    switch(sig) {
      case(SIGFPE)  :  std::cerr << "Floating point exception occured. Exiting" << std::endl;
                       // we have to unmask the signal otherwie program will slow down
                       delete helios;
                       abort();
                       //signal(SIGFPE, SIG_IGN);
                       // now we raise SIGUSR1 which is propagetaed by mpirun to other processes
                       raise(SIGUSR2);
                       break;

                       // when SIGINT or SIGTERM appears, e.g. openmpi first propagates SIGTERM too all procecess, waits a couple
                       // of seconds and then exists, (we should take care that this time is enough to finish the job)
      case(SIGINT)   : std::cerr << "SIGINT received, raising SIGTERM" << std::endl;
                       control_triggered_signal |= SIGINT;
                 //      raise(SIGTERM); 
                       break;
      case(SIGTERM)  : //helios->runningException(HELIOS_EXIT); 
                       std::cout << "SIGTERM" << std::endl;
                       control_triggered_signal |= SIGTERM;
                       break;
      case(SIGUSR1) :  writeMessage("SIGUSR1 received");
                       control_triggered_signal |= SIGUSR1;
                       break;
      case(SIGUSR2) :  writeMessage("SIGUSR2 received");
                       control_triggered_signal |= SIGUSR2;
                       break;
      default       :  std::cerr << "Unkown signal .... Ignoring\n";
    }
	
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
	    _mm_setcsr( _MM_MASK_MASK &~  (_MM_MASK_OVERFLOW|_MM_MASK_INVALID|_MM_MASK_DIV_ZERO) );
#endif







};



