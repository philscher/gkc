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
          cmplxd dB=0;
          for(int z=NzLlD; z<= NzLuD;z++){ for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) {  for(int x=NxLlD; x<= NxLuD;x++) { 

//          const double dAp_dy     = (8. *(fields->Field(x, y+1, z  , 1,1, Field::Ap) - Field(x  , y-1, z  , 1,1,Field::Ap))  -1. *(Field(x  , y+2, z  , 1,1,Field::Ap) - Field(x  , y-2, z  , 1,1,Field::Ap)))/(12.*dy)  ;  
   //       const double dAp_dx     = (8. *(Field(x+1, y, z  , 1,1, Field::Ap) - Field(x-1  , y, z  , 1,1,Field::Ap))  -1. *(Field(x+2  , y, z  , 1,1,Field::Ap) - Field(x-2, y  , z  , 1,1,Field::Ap)))/(12.*dx)  ;  

            dB += fields->Field(x,y_k,z,1,1,Field::Ap);
  //          std::cout << - dAp_dx;
          }}}
 
             dBFile << timing.time << " "<<  dB << " " << fields->Field(5,5,5,1,1,Field::Ap) << " " <<  fields->Field(3,4,5,1,1,Field::Ap) << std::endl;
    }



    };





};

#endif // _BENCHMARK_H__

