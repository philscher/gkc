
/*
 * =====================================================================================
 *
 *       Filename: Vlasov_Aux.h
 *
 *    Description: Implementation of GK Vlasov's equation for two dimensional
 *                 slab geometries
 *
 *         Author: Paul P. Hilscher (2011-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include "Global.h"

#ifndef __VLASOV_AUX_H
#define __VLASOV_AUX_H


#include "Vlasov_Cilk.h"

#include "Geometry/Geometry2D.h"

class Event;


/**
*  @brief  Vlasov Solver Implementation for two-dimensional Geometry
*          and other special reduced Geometries.
*
*  Making extensively use of Intel's Cilk Plus Array Notation to faciliate
*  array operations (especially vectorization). (http://software.intel.com/en-us/articles/intel-cilk-plus/)
*  Supported by Intel(12.1)  and GCC (svn side branch)
*
*
**/
class VlasovAux : public VlasovCilk {



   friend class Event;
        
   /**
   *   
   *   @brief Solves gyro-kinetic equation in 2-d plane in sheared geometry. 
   *
   *   For the derivation, the gyro-kinetic equation in slab geometry 
   *    
   *    \f[ 
   *          \frac{g_1}{\partial t} = 
   *          - \frac{B_0}{B_0^\star} \left( \omega_n + \omega_T \left(v^2+\mu B_t \right) \right) 
   *            \frac{\partial \chi}{\partial y} F_{0\sigma} - \alpha v_\parallel \frac{\partial G}{\partial z} 
   *    \f]
   *
   *   Sheared geometry is expressed as \f$ B_0 = \left( 0, -x/L_s, 1 \right) \f$ with L_s the shearing length. The parallel component
   *   along the magnetic field line is thus calculated according to 
   *   \f[ 
   *        k_\parallel = B_0 * k = \left( 0, -x/L_s, 1 \right) \cdot (k_x, k_y, k_z) = k_z - x \hat{s} k_y 
   *   \f]
   *   where in our normalization \f$ \hat{s} \f$  is defined as \f$ \hat{s} = 1/L_s \f$ . 
   *   As \f$ k_z \ll k_y \f$ we assumtion \f$ k_\parallel = x \hat{s} k_y \f$ is valid. So in the Vlasov equation the z-derivative
   *   is replaced by \f$ \partial_z = \hat{s} x \partial_y \f$.
   *
   *   Note : Thus is not field lign average, thus k_z corresponds to z-direction which is NOT along the magnetic field line
   *
   *  References :
   *
   *          Wang, Z. X.; Li, J. Q.; Kishimoto, Y.; Dong, J. Q.;  Magnetic-island-induced ion temperature gradient mode ; PoP (2009)
   *          Dong, J. Q.; Guzdar, P. N.; Lee, Y. C.            ;  Finite beta effects on ion temperature gradient driven modes ; Phys.of Fluid (1987)
   *
   *
   **/
   void    Vlasov_ES(
                           const CComplex  fs       [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex fss       [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex  f0 [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex  f1 [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex ft        [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Coll[NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Fields[Nq][NsLD][NmLD][NzLB][Nky][NxLB+4],
                           CComplex nonLinear                  [Nky][NxLD  ][NvLD],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           const double dt, const int rk_step, const double rk[3]);
   /**
   *    @brief two-dimensional sheared slab goemetry (electro-magnetic)
   *
   **/
   void Vlasov_EM(
                           CComplex fs       [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex fss      [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f0 [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f1 [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex ft       [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Coll      [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Fields[Nq][NsLD][NmLD][NzLB][Nky][NxLB+4]      ,
                           CComplex    nonLinear               [Nky][NxLD  ][NvLD],
                           CComplex Xi       [NzLB][Nky][NxLB][NvLB],
                           CComplex G        [NzLB][Nky][NxLB][NvLB],
                           const double dt, const int rk_step, const double rk[3]);

  
   /**
   *    Please Document Me !
   *
   **/
   void  Landau_Damping(
                           CComplex fs       [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex fss      [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f0 [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f1 [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex ft       [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Field[Nq][NsLD][NmLD][NzLB][Nky][NxLB+4],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           const double dt, const int rk_step, const double rk[3]);


   /**
   *   @brief Calculates the parallel non-linearity term for two-dimensional geometry
   *
   *   \f[
   *      \mathcal{N}{}_{v_\parallel} = 
   *      \alpha_\sigma \frac{q_\sigma}{m_\sigma} 
   *      k_\parallel \left<\phi\right> \frac{\partial g_{1\sigma}}{\partial v_\parallel} 
   *     
   *   \f]
   *
   *   @note only electrostatic calculation is supported
   *
   **/
   void calculateParallelNonLinearity(
                                const CComplex f          [NsLD][NmLD][NzLB][Nky][NxLB   ][NvLB],
                                const CComplex Fields [Nq][NsLD][NmLD ][NzLB][Nky][NxLB+4], // in case of e-s
                                const int z, const int m, const int s                     ,
                                CComplex NonLinearTerm[Nky][NxLD][NvLD]);
   
   void calculateParallelNonLinearity2(
                                const CComplex f          [NsLD][NmLD][NzLB][Nky][NxLB   ][NvLB],
                                const CComplex Fields [Nq][NsLD][NmLD ][NzLB][Nky][NxLB+4], // in case of e-s
                                const int z, const int m, const int s                     ,
                                CComplex NonLinearTerm[Nky][NxLD][NvLD]);


  public:

   /**
   *    Please Document Me !
   *
   **/
   VlasovAux(Grid *_grid, Parallel *_parallel, Setup *_setup, FileIO * fileIO, Geometry *_geo, FFTSolver *fft, Benchmark *bench, Collisions *coll); 
        
   /**
   *    Please Document Me !
   *
   **/
   void solve(std::string equation_tyoe, Fields *fields, CComplex *fs, CComplex *fss, 
              double dt, int rk_step, const double rk[3]);
 
  protected :
 
 
   Geometry2D *geo;


   /**
   *    Please Document Me !
   *
   **/
   //void printOn(std::ostream &output) const;

};


#endif //  __VLASOV_AUX_H


