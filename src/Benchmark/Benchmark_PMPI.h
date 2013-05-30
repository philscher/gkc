/*
 * =====================================================================================
 *
 *       Filename: Benchmark_PMPI.h
 *
 *    Description: Benchmarking using MPI-Profiling interface 
 *
 *         Author: Paul P. Hilscher (2013-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#ifndef __GKC_BENCHMARK_PMPI_H__
#define __GKC_BENCHMARK_PMPI_H__

#include <vector>

#include "Parallel/Parallel.h"
#include "FileIO.h"


/**
*   @brief Benchmarking using MPI-Profiling interface
*
*
**/
class Benchmark_PMPI {

  friend Parallel;
public:
  struct pmpi_vec {
     double t        ; ///< time stamp
     double dt       ; ///< time duration
     int    dest; ///< destination rank
     int    bytes    ; ///<
     char   cmd[16]  ; ///< MPI-Command
///     int    id       ; ///< 
  } _pmpi_vec;
   
  static std::vector<pmpi_vec> time_trace;

  Parallel *parallel;
  FileIO   *fileIO;

  hid_t pmpiGroup;
  static double time_start;

  int count_VlasovBoundary;

  void initData();

 public:
  Benchmark_PMPI(Setup *setup, Parallel *parallel, FileIO *fileIO);
 ~Benchmark_PMPI();
  //void start(std::string);
  //void stop ();

};

#endif // __GKC_BENCHMARK_PMPI_H__

