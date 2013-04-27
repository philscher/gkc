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

#include "Parallel/Parallel.h"
#include "FileIO.h"

/**
*   @brief Benchmarking using MPI-Profiling interface
*
*
**/
class Benchmark_PMPI {
   friend Parallel;

   Parallel *parallel;
   int count_VlasovBoundary;

 public:
   Benchmark_PMPI(Setup *setup, Parallel *parallel, FileIO *fileIO);

//   void start(std::string);
//   void stop ();

};

#endif // __GKC_BENCHMARK_PMPI_H__

