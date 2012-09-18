#include "Global.h"

#include<string>
#include<iostream>
#include <unistd.h>

#include "GKC.h"


#include <cstdlib>
   
#define MASTER_PROCESS 0

int process_rank;
GKC *gkc;


int main(int argc, char **argv)
{
    // Some internal variables
    time_t start_sim_time, start_all_time, end_sim_time, end_all_time;
    

    // Rank can be changed by MPI Process
    process_rank = 0;
    
    time(&start_all_time);

   /***************************************************************************
   * Check command line parameters (note, the user can also 
   * speify PETSc parameters
   */
   int c;
   extern char *optarg;
   
   std::string setup_filename(""), setup_Xoptions(""), setup_ExArgv("");
   // set it to automatic ?
   std::string setup_decomposition("Auto Decomposition"), setup_scalingFileName("ScalingTime.txt");
   int gkcFlags=0;
   
   while ((c = getopt(argc, argv, "x:o:c:d:v;s:;f;i")) != -1) {
            
        switch(c) {
            case 'c' : setup_filename      = std::string(optarg);
                       break;
            case 'o' : setup_Xoptions      = std::string(optarg);
                       break;
            case 'x' : setup_ExArgv        = std::string(optarg);
                       break;
            case 'd' : 
#ifdef GKC_PARALLEL_MPI            
                       setup_decomposition = std::string(optarg);
#else
            std::cout <<  "ERROR :  Decompostion only available when compiled with MPI support ... exiting." << std::endl;
            exit(1);
#endif            
                       break;
            case 'v' : gkcFlags      |= Setup::GKC_VERBOSE;
                       break;
            case 's' : gkcFlags      |= Setup::GKC_STATISTICS;
                     setup_scalingFileName = std::string(optarg);
                       break;
            case 'f' : gkcFlags      |= Setup::GKC_OVERWRITE;
                       break;
            case 'i' : gkcFlags      |= Setup::GKC_READ_STDIN;
                       break;
            //case 'b' : gkcFlags      |= GKC_SILENT;
            //           break;

        default:
            std::cout << "Unknown option : " << optarg << " please check if compiled for this mode " << std::endl;
            exit(1);
        }
    }

    int process_pid = getpid();
    // Get simulation properties and setupurations
   Setup *setup    = new Setup(argc, argv, setup_filename, setup_decomposition, setup_Xoptions, setup_ExArgv, process_pid, gkcFlags);

    
   //////////////////    Start gkc engine and main loop /////////////////////
   gkc =  new GKC(setup);
   
  
    time(&start_sim_time);
    
    //#pragma omp parallel reduction(+:GKC_status)
    int GKC_status = gkc->mainLoop();
       
    time(&end_sim_time);
  
   // write simulation time to file, so we can late investigate it
    // in case for scaling we want to collect some basic statistics about the running
   
    if((gkcFlags & Setup::GKC_STATISTICS) && (process_rank == 0)) {
   int numThreads=1;
#ifdef GKC_PARALLEL_OPENMP
   #pragma omp parallel
   {
      numThreads = omp_get_num_threads();
   }   
#endif
        ofstream scalingFile(setup_scalingFileName.c_str(), ios::out | ios::app);
        scalingFile << setup_decomposition << ":" << numThreads << " " << end_sim_time - start_sim_time << std::endl;

        scalingFile.close();
    }


   // Clean up and print basic informations.
   delete gkc;
   delete setup;
    
    time(&end_all_time);
    if(process_rank == MASTER_PROCESS) std::cout << "Running Time : Total : " << end_all_time - start_all_time << "s   Loop Time :  " << end_sim_time - start_sim_time << "s" <<std::endl;





}
