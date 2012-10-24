/*
 * =====================================================================================
 *
 *       Filename: Vlasov_Optim.h
 *
 *    Description: Definitions of optimized algorithms for Vlasov 
 *
 *         Author: Paul P. Hilscher (2012-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */



#ifndef __VLASOV_OPTIM_H
#define __VLASOV_OPTIM_H

#include "Vlasov.h"

/**
*     @brief Optimzed implementations of various Vlasov functions
*
*     Implements optimized algorithms.
*
*     @note 
*        Here, we sacrifice readability in favor of speed.
*        Experimental and mainly broken.
*
*  Optimization used :
*     Seperat calculations of real and imaginary parts
*     Cache-blocking
*
**/
class VlasovOptim : public Vlasov {

 protected:

  typedef struct{ double re; double im; } cmplx16;

  typedef cmplx16(*A6sz)[0][0][0][0][0];
  typedef cmplx16(*A5sz)[0][0][0][0];
  typedef cmplx16(*A4sz)[0][0][0];
  typedef cmplx16(*A3sz)[0][0];


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
                           const cmplx16 Field[Nq][NsLD][NmLD][NzLB][NkyLD][NxLB+4],
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
   void solve(std::string equation_tyoe, Fields *fields, CComplex *fs, CComplex *fss, double dt, int rk_step, const double rk[3]);
 
  
};


#endif //  __VLASOV_OTIM_H

