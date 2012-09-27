#include "Global.h"

#include<string>
#include<iostream>
#include <unistd.h>

#include "GKC.h"


#include <cstdlib>
#include <sys/resource.h>
   
#define MASTER_PROCESS 0

int process_rank;


void set_stacksize(int size_Mbytes)
{
    const rlim_t kStackSize = size_Mbytes * 1024 * 1024;   // min stack size = 16 MB
    struct rlimit rl;

    int result = getrlimit(RLIMIT_STACK, &rl);
    
    std::cout << " Stack size is : " << rl.rlim_cur/(1024*1024) << " MBytes" << std::endl;
    
    if (result == 0)
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
    }
    
    result = getrlimit(RLIMIT_STACK, &rl);
    std::cout << " Stack size is : " << rl.rlim_cur/(1024*1024) << " MBytes" << std::endl;


    return;
};




int main(int argc, char **argv)
{
    // Set stack size (not necesart for jobs) set_stacksize(32);
    // Some internal variables
    time_t start_sim_time, start_all_time, end_sim_time, end_all_time;
    

    // Rank can be changed by MPI Process
    process_rank = 0;
    
    time(&start_all_time);

   /***************************************************************************
   * Check command line parameters (note, the user can also 
   * speify PETSc parameters
   */
   
   std::string setup_filename(""), setup_Xoptions(""), setup_ExArgv("");
   // set it to automatic ?
   std::string setup_decomposition("Auto Decomposition"), setup_scalingFileName("ScalingTime.txt");
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

    // Get simulation properties and setup
    Setup *setup    = new Setup(argc, argv, setup_filename, setup_decomposition, setup_Xoptions, setup_ExArgv, gkcFlags);

    
   //////////////////    Start gkc engine and main loop /////////////////////
   GKC *gkc =  new GKC(setup);
   
  
    time(&start_sim_time);
    
    //#pragma omp parallel reduction(+:GKC_status)
    int GKC_status = gkc->mainLoop();
       
    time(&end_sim_time);
  
   // write simulation time to file, so we can late investigate it
    // in case for scaling we want to collect some basic statistics about the running
  
    ////////////////////// get somer statistic information //////////////
    /*
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
    */
   //////////////////////////////////////////////////////////////////////////

   // Clean up and print basic informations.
   delete gkc;
   delete setup;
    
    time(&end_all_time);
    if(process_rank == MASTER_PROCESS) std::cout << "Running Time : Total : " << end_all_time - start_all_time << "s   Loop Time :  " << end_sim_time - start_sim_time << "s" <<std::endl;





}
