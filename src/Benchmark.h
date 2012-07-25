/*
 * =====================================================================================
 *
 *       Filename:  Benchmark.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  07/25/2012 01:50:25 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef __GKC_BENCHMARK_H
#define __GKC_BENCHMARK_H

#include "Global.h"


#include <papi.h> 

/**
* @brief interface to PAPI 
*
* In case for PAPI, are hardware counters defined ?
* Check with papi_avail
*
* PAPI, the Performance Application Programming Interface
* (@ref http://icl.cs.utk.edu/papi/),
* is a C-library to faciliate the usage of hardware counters.
*
* We mainly employ it to estimate the FLOP/s (Floating point
* operations per second).
*
* @todo optimize also cache misses
* @todo save value in HDF-5 file.
* @todo check if using PAPI does have performance impact
*
*/
class Benchmark
{
   
  struct Measure {
   float rtime;  ///< Real time
   float ptime;  ///< Processor time
   long long flpops; ///< Total floating point operations
   float mflops;   ///< Million Floating Operations Per Second
   } M;


  /**
  *    @brief get error string from PAPI error
  *
  **/ 
  std::string getPAPIErrorString(int error_val)
   {
     const int max_str_len = 128;
     char error_str[max_str_len];
     PAPI_perror(error_val, error_str, max_str_len);
     return std::string(error_str);
   };


  public:
   
   /**
   *   @brief constructor which initialized PAPI
   *
   **/
   Benchmark()
   {

int retval = PAPI_library_init(PAPI_VER_CURRENT);

   };

   ~Benchmark()
   {
         PAPI_shutdown();
   };

   void start(std::string id, int type=0)
   {  
      int ret;
      // Setup PAPI library and begin collecting data from the counters
      if((ret = PAPI_flops( &M.rtime, &M.ptime, &M.flpops, &M.mflops)) < PAPI_OK)
         check(-1, DMESG(getPAPIErrorString(ret))); 
   };

   void stop(std::string id, int type=0)
   {

      int ret;
      // Setup PAPI library and begin collecting data from the counters
      if((ret = PAPI_flops( &M.rtime, &M.ptime, &M.flpops, &M.mflops)) < PAPI_OK)
         check(-1, DMESG(getPAPIErrorString(ret)));

     std::cout << std::endl << "FLOP : " << M.flpops << " GFLOP/s : " << std::setprecision(2) << (M.mflops/1000.) << std::endl;
     std::cout << std::endl << " GFLOP/s : " << std::setprecision(2) << std::setiosflags(ios::fixed) << setw(5) << (M.mflops/1000.) << std::endl;
   };


};


#endif // __GKC_BENCHMARK_H
