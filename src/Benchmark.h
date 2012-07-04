/*
 * =====================================================================================
 *
 *       Filename:  Benchmark.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/15/2011 06:16:28 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */
#ifndef __BENCHMARK_H__
#define __BENCHMARK_H__

#include "Global.h"

#include "Setup.h"
#include "Plasma.h"
#include "Fields.h"
#include "Vlasov.h"


class Benchmark {


  ofstream dBFile;
  bool useBenchmark;
  public:
      Benchmark(Setup *setup) {
        useBenchmark = setup->get("Benchmark.Use", 0);
        if(useBenchmark == false) return;
       std::string BFileName = setup->get("Benchmark.FileName", "default.txt");
       std::cout << "Benchmark file : " << BFileName << std::endl;
       dBFile.open(BFileName.c_str());
      };

    void benchmark(Vlasov *vlasov, Fields *fields, Timing timing) {
        if(!useBenchmark) return;

        if((plasma->nfields >= 2)) {
             dBFile << timing.time << " "<< 
               real(fields->Field0(NxLlD,NkyLlD+1,NzLlD,Field::Ap )) << "  " <<
               imag(fields->Field0(NxLlD,NkyLlD+1,NzLlD,Field::Ap )) << "  " <<
               real(fields->Field0(NxLlD,NkyLlD+1,NzLlD,Field::phi)) << "  " <<
               imag(fields->Field0(NxLlD,NkyLlD+1,NzLlD,Field::phi)) <<  std::endl;
               
    }



    };





};

#endif // _BENCHMARK_H__

