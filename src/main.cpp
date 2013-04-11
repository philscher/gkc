#include "Global.h"

#include <string>
#include <iostream>
#include <unistd.h>
#include <cstdlib>

#include "GKC.h"
#include "Tools/System.h"
#define MASTER_PROCESS 0

extern int process_rank;

/**
*   @brief Program starting point
*
*   The main program. Reads in program arguments and intilizes
*   the Setup class and the gkc class and handles over the
*   control to mainLoop
*
**/
int main(int argc, char **argv)
{
  // Set stack size (necesary to set for local jobs) 
  System::set_min_stacksize(256);
  
  // Some internal variables
  time_t start_sim_time, start_all_time, end_sim_time, end_all_time;
    
  // Rank can be changed by MPI Process
  process_rank = 0;
    
  time(&start_all_time);

  //////////////////////////////////////////////////////////////////
  //
  //   Check command line parameters (note, the user can also 
  //   speify PETSc parameters
  //
  std::string setup_filename(""), setup_Xoptions(""), setup_ExArgv(""), stack_size("");
  // set it to automatic ?
  std::string setup_decomposition("Auto Decomposition");
  int gkcFlags=0;
  
  /////////////////// Read In Command Line Arguments //////////////////
  int c;
  extern char *optarg;

  while ((c = getopt(argc, argv, "x:o:c:d:v;s:;f;i")) != -1) {

    switch(c) {

      case 'c' : setup_filename      = std::string(optarg);
                       break;
      case 'o' : setup_Xoptions      = std::string(optarg);
                       break;
      case 'x' : setup_ExArgv        = std::string(optarg);
                       break;
      case 'd' : setup_decomposition = std::string(optarg);
                       break;
      case 'v' : gkcFlags           |= Setup::GKC_VERBOSE;
                       break;
      case 's' : stack_size          = std::string(optarg); 
                       break;
      case 'f' : gkcFlags           |= Setup::GKC_OVERWRITE;
                       break;
      case 'i' : gkcFlags           |= Setup::GKC_READ_STDIN;
                       break;
      //case 'b' : gkcFlags         |= GKC_SILENT;
      //               break;
      default:
      
        std::cout << "Unknown option : " << optarg << " please check if compiled for this mode " << std::endl;
        exit(1);
    }
    
  }

  // Get simulation properties and setup
  Setup *setup    = new Setup(argc, argv, setup_filename, setup_decomposition, setup_Xoptions, setup_ExArgv, gkcFlags);

  //////////////////    Start gkc engine and main loop /////////////////////

  // define static, as destructor is called @ exit() (see, § 3.6.3 of C++03) 
  static GKC *gkc =  new GKC(setup);
   
  time(&start_sim_time);
    
  int GKC_status = gkc->mainLoop();
       
  time(&end_sim_time);
  
  //////////////////////////////////////////////////////////////////////////

  // Clean up and print basic informations.
  delete gkc;
  delete setup;
    
  time(&end_all_time);

  if(process_rank == MASTER_PROCESS) {

    std::cout << "Total Time : " << end_all_time - start_all_time << "s" 
              << " Loop Time : " << end_sim_time - start_sim_time << "s" <<std::endl;
  }

  return 0;
}

