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

#include "Setup.h"
#include "Parallel.h"
#include "FileIO.h"
#include "Timing.h"


#include <sys/time.h>
#include <papi.h> 


class Vlasov;
class Fields;

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
**/
class Benchmark : public IfaceGKC
{
  int num_hwcntrs; ///< Number of available Hardware counters

  
   
  struct Counters
  {
   long long Cycles; ///< Total Cycles
   long long Instructions   ; ///< Total Instructions
  };


  struct Measure {
   float      rtime;  ///< Real time
   float      ptime;  ///< Processor time
   long long flpops;  ///< Total floating point operations
   float     mflops;  ///< Million Floating Operations Per Second
   } M;
    
  timespec ts_start, ts_end;
  
  bool useBenchmark;

  Parallel * parallel;
  /**
  *    @brief get error string from PAPI error
  *
  **/ 
  static std::string getPAPIErrorString(int error_val);

  double simMaxGFLOPS;
  public:
  
  int BlockSize_X, BlockSize_V;

 
   /**
   *   @brief constructor which initialized PAPI
   *
   **/
   Benchmark(Setup *setup, Parallel *_parallel);

   ~Benchmark();

   void start(std::string id, int type=0);

   double stop(std::string id, int type=0);
 
   void bench(Vlasov *vlasov, Fields *fields) ;

protected:
   virtual void writeData(Timing timing, double dt);
   void initDataOutput(Setup *setup, FileIO *fileIO);
   void closeData();
   virtual void printOn(ostream &output) const ;
};


#endif // __GKC_BENCHMARK_H
