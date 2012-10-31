/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/10/2011 06:22:08 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <unistd.h>
#include <fenv.h>
#include <csignal>
#include <iostream>
#include <sys/resource.h>


#include "MHD.h"



MHD *sim;

void set_stacksize(int size_Mbytes)
{
    const rlim_t kStackSize = size_Mbytes * 1024 * 1024;   // min stack size = 16 MB
    struct rlimit rl;

    int result = getrlimit(RLIMIT_STACK, &rl);
    
    if(rl.rlim_cur == RLIM_INFINITY) std::cout << " Stack size is : Unlimited" << std::endl;
    else                             std::cout << " Stack size is : " <<
                                             rl.rlim_cur/(1024*1024)  << " MBytes" << std::endl;
    //std::cout << " Stack size is : " << (rl.rlim_cur == RLIM_INFINITY ? "Unlimited" :
    
    if ((result == 0) && (rl.rlim_cur != RLIM_INFINITY))
    {
        if (rl.rlim_cur < kStackSize)
        {
            rl.rlim_cur = kStackSize;
            result = setrlimit(RLIMIT_STACK, &rl);
            if (result != 0)
            {
                fprintf(stderr, "setrlimit returned result = %d\n", result);
            }
        }
    
       result = getrlimit(RLIMIT_STACK, &rl);
       std::cout << " Stack size set to : " << rl.rlim_cur/(1024*1024) << " MBytes" << std::endl;
    }
    


    return;
};


void signal_handler(int sig)
	{
      switch(sig) {
      case(SIGFPE)  :  std::cerr << "Floating point exception occured. Exiting" << std::endl;
                       // we have to unmask the signal otherwie program will slow down
                       break;

                 //      raise(SIGTERM); 
                       break;
      case(SIGTERM)  : //helios->runningException(HELIOS_EXIT); 
                       std::cerr << "SIGTERM" << std::endl;
                       break;
      case(SIGUSR1) :  std::cerr << "SIGUSR1 received" << std::endl;
                       break;
      case(SIGUSR2) :  std::cerr << "SIGUSR2 received" << std::endl;
                       break;
      default       :  std::cerr << "Unkown signal .... Ignoring\n" << std::endl;
                       
    }
                       // if signal is catched by different threads, make sure
                       // that sim is deleted only once !
                       #pragma omp single
                       delete sim;
                       exit(0);
	
    };
	


template<class T> T max(T a, T b) { return a > b ? a : b; };


int main(int argc, char **argv)
{
    // Set stack size (not necesart for jobs) 
    set_stacksize(max(32, MHD::Nx*MHD::Nky/1000));

    feenableexcept(FE_DIVBYZERO | FE_INVALID);
    signal(SIGFPE,   signal_handler);
    signal(SIGINT ,  signal_handler);
    signal(SIGUSR1 , signal_handler);
    signal(SIGUSR2 , signal_handler);


   
   std::string setup_filename(""), setup_Xoptions("");

    int c;
   while ((c = getopt(argc, argv, "c:o:")) != -1) {
            
        switch(c) {
            case 'c' : setup_filename      = std::string(optarg);
                       break;
            case 'o' : setup_Xoptions      = std::string(optarg);
                       break;
        default:
            std::cout << "Unknown option : " << optarg << " please check if compiled for this mode " << std::endl;
            exit(1);
        }
    }

    // Start as ./MHD -c File -o Options
    
    /***********************************************************
   // Start Helios engine and main loop
   ***********************************************************/
    sim = new MHD(setup_filename, setup_Xoptions);
    sim->startMainLoop();



    delete sim;
}



