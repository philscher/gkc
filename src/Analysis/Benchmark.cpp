/*
 * =====================================================================================
 *
 *       Filename: Benchmark.cpp
 *
 *    Description: PAPI Benchmark class implementation.
 *
 *         Author: Paul P. Hilscher (2012-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include <functional>

#include "Benchmark.h"
#include "Vlasov/Vlasov.h"
#include "Fields/Fields.h"

Benchmark::Benchmark(Setup *setup, Parallel *_parallel, FileIO *fileIO) : parallel(_parallel)
{
  
  int retval;
  
  BlockSize_X=1, BlockSize_V=1;

  useBenchmark = setup->get("Benchmark.usePAPI", 0);

  if(useBenchmark) {

   // Initialize the PAPI library 
   if(PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT) DMESG("PAPI library init error!\n");
   
   // check avaiable counters
   check((num_hwcntrs = PAPI_num_counters()) <= PAPI_OK, DMESG("Could not get number of avaiable counters"));
  
   // Set Hardware counters here (use setup for that !)
   events[:] =  PAPI_NULL;
   events[0] =  PAPI_DP_OPS;
   events[1] =  PAPI_VEC_DP;
   events[2] =  PAPI_FP_OPS;

  }
 
  initData(setup, fileIO);
}
   
Benchmark::~Benchmark()
{
  H5Gclose(benchGroup);
 
  closeData();

  delete eventTable;
  if(useBenchmark) PAPI_shutdown();
}


void Benchmark::bench(Vlasov *vlasov, Fields *fields) 
{
  if(!useBenchmark) return;
  return;    
    // Set difference optimization options
    for(BlockSize_X = 1; BlockSize_X < Nx/8; BlockSize_X += 4) {
    for(BlockSize_V = 1; BlockSize_V < Nv/8; BlockSize_V += 4) {
        
      auto benchtest = [=] (std::function<void ()> func) -> double { 
        
        double gflops   = 0.;
        const int niter = 16;

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
               
    std::cout << "Vlasov plain " << benchtest(vlasov_plain) << std::endl;
    std::cout << "Vlasov full  " << benchtest(vlasov_full ) << std::endl;
    std::cout << "Poisson      " << benchtest(poisson     ) << std::endl;
    double ggflops = benchtest(full_step  );
    //simMaxGFLOPS   = std::max(simMaxGFLOPS, ggflops);
    std::cout << "One Step     " << ggflops << std::endl;
   
   } } // BlockSize_V, BlockSize_X
}

 
void Benchmark::writeData(const Timing &timing, const double dt) 
{


}
 
void Benchmark::closeData()
{


}


void Benchmark::printOn(std::ostream &output) const
 {

  //double totalFLOPS    = parallel->reduce(simMaxGFLOPS, Op::sum) * parallel->numThreads;
  //double flops_per_cpu = totalFLOPS/(parallel->numProcesses * parallel->numThreads) ;

  double totalFLOPS = 0.;
  output << "PAPI       | Total GFLOPS : " << std::setprecision(3) << totalFLOPS << std::endl;
  //       <<     "   Average GFLOPS/CPU : " << flops_per_cpu        << std::endl;

}

  
std::string Benchmark::getPAPIErrorString(int error_val)
{

  const int max_str_len = 128;
  char error_str[max_str_len];
     //PAPI_perror(error_val, error_str, max_str_len);
  return std::string(error_str);

}

   
void Benchmark::start(std::string id, int type)

{
  if(!useBenchmark) return;
   // measure floating point
   PAPI_start_counters(events, event_num);
   time_usec_start = PAPI_get_real_usec();
}

double Benchmark::stop(std::string id, int type)
{
  if(!useBenchmark) return;

  long long time_usec_end = PAPI_get_real_usec();
 
  // Get results and time
  if(PAPI_read_counters(event.value, event_num) != PAPI_OK) std::cout << "Error Reading counters" << std::endl;
  
  event.dtime = (time_usec_end - time_usec_start) * 1.e-6;
  
  eventTable->append(&event);

  if(PAPI_stop_counters(event.value, event_num) != PAPI_OK) std::cout << "Error Stoping counters" << std::endl;
    

  return 0;
}


void Benchmark::initData(Setup *setup, FileIO *fileIO)
{
  if(!useBenchmark) return;
   
  benchGroup = fileIO->newGroup("/Benchmark");

  check(H5LTset_attribute_int   (benchGroup, ".", "NumberOfCounters",  &num_hwcntrs, 1), DMESG("H5LTset_attribute"));


  //////////////////////// Set Table for Events
  size_t cs_offset[event_num+1]; cs_offset[0] = HOFFSET(Event, dtime); 
  for(int n=0; n < event_num; n++) cs_offset[n+1] = HOFFSET( Event, value[0]) + n * sizeof(long long);

  size_t cs_sizes[event_num+1]; cs_sizes[0] = sizeof(double);
  for(int n=0; n < event_num; n++) cs_sizes[n+1]  = sizeof(long long);

  hid_t cs_types[event_num+1];  cs_types[0] = H5T_NATIVE_DOUBLE;
  for(int n=0; n < event_num; n++) cs_types[n+1]  = H5T_NATIVE_LLONG;
 
  // Set PAPI event names as columns
  const char *cs_names[event_num+1]; cs_names[0] = std::string("dt").c_str();
  
  // that looks weird, but I do not find any information how to call event_code_to_name properly ... (HELP!) 
  char str[event_num][128]; // MAX_LEN ? Maybe 20 or so ...
  for(int n=0; n < event_num; n++) {
    PAPI_event_code_to_name(events[n], &str[n][0]);
    cs_names[n+1] = &str[n][0];
  }
  
  Event ev;
  eventTable = new TableAttr(benchGroup, "PAPICounters", event_num + 1, cs_names, cs_offset, cs_types, cs_sizes, &ev); 
}

void Benchmark::save(std::string id, int value)
{
  check(H5LTset_attribute_int(benchGroup, ".", id.c_str(),  &value, 1), DMESG("H5LTset_attribute"));
}


