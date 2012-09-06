/*
 * =====================================================================================
 *
 *       Filename: Vlasov_Cilk.h
 *
 *    Description: Implementation of GK Vlasov's equation using 
 *                 Intel Cilk (Array Notation)
 *
 *         Author: Paul P. Hilscher (2011-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include "Vlasov.h"


#ifndef __VLASOV_OPTIM_H
#define __VLASOV_OPTIM_H




#include "Vlasov.h"
#include "Global.h"

class Event;


typedef struct{ double re; double im; } cmplx16;


/**
*  @brief Implementation of Vlasov's equation using Cilk Array Notation
*
*  Making extensively use of Intel's Cilk Plus Array Notation to faciliate
*  array operations (especially vectorization). (http://software.intel.com/en-us/articles/intel-cilk-plus/)
*  Supported by Intel(12.1)  and GCC (svn side branch)
*
*
**/
class VlasovOptim : public Vlasov {


  friend class Event;
        
   /**
   *   
   *   @brief Solves gyro-kinetic equation in 2-d plane in sheared geometry.
   *
   *   Write about optimizations
   *
   *
   **/
   void    Vlasov_2D(
                           const cmplx16 fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplx16 fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplx16 f0 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplx16 f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplx16 ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplx16 phi[NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           cmplx16 nonLinear[NzLD][NkyLD][NxLD][NvLD],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           const double dt, const int rk_step, const double rk[3]);
  
  public:

   /**
   *    Please Document Me !
   *
   **/
   VlasovOptim(Grid *_grid, Parallel *_parallel, Setup *_setup, FileIO * fileIO, Geometry *_geo, FFTSolver *fft, Benchmark *bench); 
        
   /**
   *    Please Document Me !
   *
   **/
   int solve(std::string equation_tyoe, Fields *fields, Array6C fs, Array6C fss, double dt, int rk_step, const double rk[3]);
 
  protected :
  
};


#endif //  __VLASOV_OTIM_H

