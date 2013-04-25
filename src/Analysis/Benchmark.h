/*
 * =====================================================================================
 *
 *       Filename: Benchmark.h
 *
 *    Description: PAPI Benchmark class definition. Used for optimizing
 *                 Vlasov equation runtime.
 *
 *         Author: Paul P. Hilscher (2012-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __GKC_BENCHMARK_H__
#define __GKC_BENCHMARK_H__

#include "Global.h"

#include "Setup.h"
#include "Parallel/Parallel.h"
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
* @todo check if using PAPI does have performance impact
*
**/
class Benchmark : public IfaceGKC
{
  hid_t benchGroup; ///< HDF-5 Group id 
  
  int num_hwcntrs; ///< Number of available Hardware counters

  TableAttr *eventTable;
  
  static const int event_num = 3;
 
  struct Event {
      double    dtime;            /// Time 
      long long value[event_num]; // Values
  } event;


  long long time_usec_start;

  int events[8]      ; /// Pre-defined max number of events
  char *event_name[8]; /// Name of events

  bool useBenchmark;

  Parallel * parallel;
  /**
  *    @brief get error string from PAPI error
  *
  **/ 
  static std::string getPAPIErrorString(int error_val);

 public:
  
  int BlockSize_X, ///< Blocksize in X to use in Vlasov equation 
      BlockSize_V; ///< Blocksize in Y to use in Vlasov equation

  /**
  *   @brief constructor which initialized PAPI
  *
  **/
  Benchmark(Setup *setup, Parallel *_parallel, FileIO *fileIO);

 ~Benchmark();

  void start(std::string id, int type=0);
  double stop(std::string id, int type=0);
  
  void bench(Vlasov *vlasov, Fields *fields);

  void save(std::string id, int value);

 protected:

  virtual void writeData(const Timing &timing, const double dt);
  void initData(Setup *setup, FileIO *fileIO);
  void closeData();
  virtual void printOn(std::ostream &output) const ;

};


#endif // __GKC_BENCHMARK_H__
