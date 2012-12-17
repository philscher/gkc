#include <functional>

#include "Benchmark.h"
#include "Vlasov/Vlasov.h"
#include "Fields/Fields.h"

Benchmark::Benchmark(Setup *setup, Parallel *_parallel) : parallel(_parallel)
{
      
   int EventSet = PAPI_NULL;
   int retval;
  

   BlockSize_X=1, BlockSize_V=1;

   useBenchmark = setup->get("Benchmark.usePAPI", 0);

   if(useBenchmark) {
   //int Events[2] = { PAPI_TOT_CYC, PAPI_TOT_INS };
/*  
   // Initialize the PAPI library 
   retval = PAPI_library_init(PAPI_VER_CURRENT);
   check(retval != PAPI_VER_CURRENT, DMESG("PAPI library init error!\n"));
   
   // check avaiable counters
   check((num_hwcntrs = PAPI_num_counters()) <= PAPI_OK, DMESG("Could not get number of avaiable counters"));

   // Create our event set
   check(PAPI_create_eventset(&EventSet) != PAPI_OK, DMESG("Could not create Eventset"));
   check(PAPI_add_event(EventSet, PAPI_TOT_INS) != PAPI_OK, DMESG("Could not add Event"));
*/
// reset counters PAPIF_stop_counters(NULL, array_length, check)
   }
};
   
Benchmark::~Benchmark()
   {
      if(useBenchmark) PAPI_shutdown();
   };


void Benchmark::bench(Vlasov *vlasov, Fields *fields) 
{
      if(!useBenchmark) return;
      
     // Set difference optimization options
     for(BlockSize_X=1; BlockSize_X<Nx/8; BlockSize_X+=4) {
     for(BlockSize_V=1; BlockSize_V<Nv/8; BlockSize_V+=4) {
        

         auto benchtest = [=] (std::function<void ()> func) -> double { 
        
            double gflops = 0.;
            const int niter=16;

       start("Vlasov"); 
       for(int i=0;i<niter;i++) func() ; 
       gflops = stop("Vlasov"); //((double) niter);
      
            return gflops;
    };

     auto vlasov_plain = [=] (void) { 
             const double rk[] = { 1., 2., 0.};
      vlasov->solve(vlasov->getEquationType(), fields, vlasov->fss , vlasov->f, 1.e-3 , 2, rk); 
     };
     
          auto vlasov_full = [=] (void) { 
             const double rk[] = { 1., 2., 0.};
      vlasov->solve(fields, vlasov->fss , vlasov->f , 1.e-3 , 2, rk); 
     };
          
     auto poisson = [=] (void) { 
           fields->solve(vlasov->f0,vlasov->fs);
     };

     auto full_step = [=] (void) {
             vlasov_full(); poisson(); 
    };
               
    // auto full_step = [=] (void) {
         //   double dt = timeIntegration->solveTimeStep(vlasov, fields, particles, timing);        
    //};


          std::cout << "Vlasov plain " << benchtest(vlasov_plain) << std::endl;
          std::cout << "Vlasov full  " << benchtest(vlasov_full ) << std::endl;
          std::cout << "Poisson      " << benchtest(poisson     ) << std::endl;
          double ggflops = benchtest(full_step  );
          simMaxGFLOPS =  max(simMaxGFLOPS, ggflops);
          std::cout << "One Step     " << ggflops << std::endl;
   
       }
   }
   

};

 
void Benchmark::writeData(const Timing &timing, const double dt) 
{


}
 
void Benchmark::closeData()
{


}

void Benchmark::printOn(std::ostream &output) const
 {

  double totalFLOPS = parallel->reduce(simMaxGFLOPS, Op::SUM) * parallel->numThreads;
  //output << "Benchmark  | using PAPI  " << PAPI_VERSION_MAJOR << "." << PAPI_VERSION_MINOR <<  "   Available Counters : " << num_hwcntrs << std::endl;
  output << "           | Total GFLOPS : " << std::setprecision(3) << totalFLOPS    << "   " <<
                        " Average GFLOPS/CPU : "  << totalFLOPS/(parallel->numProcesses * parallel->numThreads) << std::endl;
                          

}


  
std::string Benchmark::getPAPIErrorString(int error_val)
  {
     const int max_str_len = 128;
     char error_str[max_str_len];
     //PAPI_perror(error_val, error_str, max_str_len);
     return std::string(error_str);
  };

   
void Benchmark::start(std::string id, int type)

{ 
      if(!useBenchmark) return; 
      int ret;
      // Setup PAPI library and begin collecting data from the counters
      if((ret = PAPI_flops( &M.rtime, &M.ptime, &M.flpops, &M.mflops)) < PAPI_OK)
      ;//   check(-1, DMESG(getPAPIErrorString(ret))); 
    clock_gettime(CLOCK_REALTIME, &ts_start); 

};

double Benchmark::stop(std::string id, int type)
   {
      if(!useBenchmark) return 0.; 

      int ret;
      // Setup PAPI library and begin collecting data from the counters
      if((ret = PAPI_flops( &M.rtime, &M.ptime, &M.flpops, &M.mflops)) < PAPI_OK)
      ;//   check(-1, DMESG(getPAPIErrorString(ret)));
    
      clock_gettime(CLOCK_REALTIME, &ts_end); 
    
       double secs  = (ts_end.tv_sec - ts_start.tv_sec);
       double nsecs = (ts_end.tv_nsec - ts_start.tv_nsec);
    
       double time = secs + nsecs * 1.e-9;

       long long flop = M.flpops;//parallel->reduce(M.flpops);

       const double gflops = M.mflops/1.e3;

       return gflops;
   //  std::cout << std::endl << id << " @ GFLOP/s : " << std::setprecision(2) << std::setiosflags(ios::fixed) << setw(5) << (M.mflops/1.e3) << std::endl;
  //   std::cout << std::endl << id << "M@ GFLOP/s : " << std::setprecision(2) << std::setiosflags(ios::fixed) << setw(5) << (flop/1.e9/time) << std::endl;
   };


    
void Benchmark::initData(Setup *setup, FileIO *fileIO)
{
   // create HDF-5 table
    
   hid_t benchGroup = check(H5Gcreate(fileIO->getFileID(), "/Benchmark",H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), DMESG("Error creating group for Geometry : H5Gcreate"));

   check(H5LTset_attribute_int   (benchGroup, ".", "NumberOfCounters",  &num_hwcntrs, 1), DMESG("H5LTset_attribute"));

         
   //////////////////////// Set Table for species.
   size_t counters_offset[]     = { HOFFSET(Counters , Cycles ), HOFFSET( Counters, Instructions ) };
   size_t counters_sizes[]      = { sizeof(long long), sizeof(long long) };
   hid_t counters_type[]        = { H5T_NATIVE_LLONG, H5T_NATIVE_LLONG };
   const char *counters_names[]  = { "Cycles", "Intructions" };

   Counters counters;

   check(H5TBmake_table("Counters", fileIO->getFileID(), "Counters", (hsize_t) 2, (hsize_t) 0, 
                        sizeof(Counters), (const char**) counters_names, counters_offset, counters_type,
                         32, NULL, 0, &counters ), DMESG("H5Tmake_table Counters"));
   
   // create table for all included species
   //H5TBappend_records (fileIO->getFileID(), "Species", 1, sizeof(Species), species_offset, species_sizes, &species[s]); 
    H5Gclose(benchGroup);

};

