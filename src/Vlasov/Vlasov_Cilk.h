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

//#ifdef GKC_CILK

#include "Vlasov.h"


#ifndef __VLASOV_CILK_H
#define __VLASOV_CILK_H



#include "Vlasov.h"
#include "Global.h"

#include<iostream>
#include<fstream>

class Event;


/**
*  @brief Implementation of Vlasov's equation using Cilk Array Notation
*
*  Making extensively use of Intel's Cilk Plus Array Notation to faciliate
*  array operations (especially vectorization). (http://software.intel.com/en-us/articles/intel-cilk-plus/)
*  Supported by Intel(12.1)  and GCC (svn side branch)
*
*
**/
class VlasovCilk : public Vlasov {

  friend class Event;
        
   /**
   *    Please Document Me !
   *
   **/
   void setBoundaryXY(Array3d A, int dir=DIR_XY) ;
        
	Array3d transform2Real(Array5z A, const int m, const int s);

   /**
   *    Please Document Me !
   *
   **/
   Array3d transform2Real(Array6z A, const int v, const int m, const int s);
   
   // Some temporary arrays (not all are necesserally intialized !)
   Array3z dphi_dy, dphi_dx, dAp_dx, dAp_dy;
  
   /**
   *    Please Document Me !
   *
   **/
   Array4z k2p_phi;
  
   /**
   *    Please Document Me !
   *
   **/
   Array4z calculatePoissonBracket(Array5z A, Array6z B, const int m, const int s);

   /**
        needed for nonLinear transforms
   *
   **/
   Array4z nonLinearTerms; 
  
   /**
   *    Please Document Me !
   *
   **/
   Array3d xy_dphi_dx, xy_dphi_dy, xy_df1_dy, xy_df1_dx, xy_phi, xy_f1;

   /**
   *    Please Document Me !
   *
   **/
   Array3d  SendXuR, SendXlR, RecvXuR, RecvXlR;

   /**
   *   
   *   @brief Solves gyro-kinetic equation in 2-d plane in sheared geometry. 
   *
   *   For the derivation, the gyro-kinetic equation in slab geometry 
   *    
   *    \f[ 
   *          \frac{g_1}{\partial t} = 
   *          - \frac{B_0}{B_0^\star} \left( \omega_n + \omega_T \left(v^2+\mu B_t \right) \right) 
   *            \frac{\partial \Xi}{\partial y} F_{0\sigma} - \alpha v_\parallel \frac{\partial G}{\partial z} 
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
   void    Vlasov_2D(
                           cmplxd fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd vf0[NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd phi[NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           cmplxd k2_phi[plasma->nfields][NzLD][NkyLD][NxLD],
                           cmplxd nonLinear[NzLD][NkyLD][NxLD][NvLD],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           Fields *fields,
                           const double dt, const int rk_step, Array6z _fs);
   /**
   *   @brief perform full-f simulations
   *
   **/
   void    Vlasov_2D_Fullf(
                           cmplxd fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd vf0[NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd phi[NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           cmplxd k2_phi[plasma->nfields][NzLD][NkyLD][NxLD],
                           Fields *fields,
                           const double dt, const int rk_step, Array6z _fs);


   /**
   *    Please Document Me !
   *
   **/
   void  Vlasov_2D_Island(
                           cmplxd fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd vf0[NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd phi[NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           cmplxd k2_phi   [plasma->nfields][NzLD][NkyLD][NxLD],
                           cmplxd nonLinear[NzLD][NkyLD][NxLD][NvLD],
                           cmplxd dphi_dx  [NzLB][NkyLD][NxLB],
                           const double    X[NxGB], const double V[NvGB], const double M[NmGB],
                           Fields *fields,
                           const double dt, const int rk_step, Array6z _fs);

   /**
   *    Please Document Me !
   *
   **/
   void  Landau_Damping(
                           cmplxd fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd vf0[NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd phi[NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           //cmplxd k2_phi[NzLD][NkyLD][NxLD],
                           cmplxd k2_phi[plasma->nfields][NzLD][NkyLD][NxLD],
                           cmplxd dphi_dx[NzLB][NkyLD][NxLB],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           Fields *fields,
                           const double dt, const int rk_step);

   /**
   *
   *   @brief Calculated the gyro-average modified potential
   *
   *   \f[ 
   *        \bar{\Xi} = \bar{\phi}_1 - \frac{v_\parallel}{c}\bar{A}_{1\parallel} 
             *         + \frac{\mu}{q_\sigma} \bar{B}_{1\parallel}
   *   \f]
   *
   *  as well as the Gamma abbrevation functions (Eq. 2.50 below)
   *
   *  \f[ 
   *       \Gamma_{\sigma,\nu} = \partial_\nu F_{1\sigma} + 
   *       \frac{F_{0\sigma}}{T_{0\sigma}} \partial_\nu \left( q_\sigma \bar{\phi}_1 + \mu \bar{B}_{1\parallel} \right)
   *  \f]
   *
   *  and F1 which needs to be reconstructed from g and is needed for the magnetic mirror term (Eq. 2.50)
   *
   *  \f[ 
   *     F_1 = g_{1\sigma} - \frac{q_sigma}{c}\bar{A}_{1\parallel}\frac{v_\parallel}{T_{0\sigma}} F_{0\sigma} 
   *  \f]
   *
   *   Reference : @cite Note from Goerles PhD Thesis (Eq. 2.32)
   *
   **/
   void setupXiAndG(
                           const cmplxd g       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd f0       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd phi       [NsLD][NmLD][NzLB][NkyLD][NxLB+4  ],
                           const cmplxd Ap       [NsLD][NmLD][NzLB][NkyLD][NxLB+4  ],
                           const cmplxd Bp       [NsLD][NmLD][NzLB][NkyLD][NxLB+4  ],
                           cmplxd Xi       [NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd G        [NzLB][NkyLD][NxLB  ][NvLB],
                           const double V[NvGB], const double M[NmGB],
                           const int m, const int s);


   /**
   *    Please Document Me !
   *
   **/
   void Vlasov_EM(
                           cmplxd fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd vf0[NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd phi[NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           const cmplxd Ap [NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           const cmplxd Bp [NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           //cmplxd k2_phi[NzLD][NkyLD][NxLD],
                           cmplxd k2_phi[plasma->nfields][NzLD][NkyLD][NxLD],
                           cmplxd dphi_dx[NzLB][NkyLD][NxLB],
                           cmplxd Xi       [NzLD][NkyLD][NxLD  ][NvLD],
                           cmplxd G        [NzLD][NkyLD][NxLD  ][NvLD],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           Fields *fields,
                           const double dt, const int rk_step);


   /**
   *    Please Document Me !
   *
   **/
   void Vlasov_EM_2D(
                           cmplxd fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd vf0[NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd phi[NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           const cmplxd Ap [NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           const cmplxd Bp [NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           //cmplxd k2_phi[NzLD][NkyLD][NxLD],
                           cmplxd k2_phi[plasma->nfields][NzLD][NkyLD][NxLD],
                           cmplxd dphi_dx[NzLB][NkyLD][NxLB],
                           cmplxd Xi       [NzLD][NkyLD][NxLD  ][NvLD],
                           cmplxd G        [NzLD][NkyLD][NxLD  ][NvLD],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           Fields *fields,
                           const double dt, const int rk_step);



  public:

   /**
   *    Please Document Me !
   *
   **/
   VlasovCilk(Grid *_grid, Parallel *_parallel, Setup *_setup, FileIO * fileIO, Geometry<GKC_GEOMETRY> *_geo, FFTSolver *fft); 
        
   /**
   *    Please Document Me !
   *
   **/
   int solve(std::string equation_tyoe, Fields *fields, Array6z fs, Array6z fss, double dt, int rk_step, int user_boundary_type=BOUNDARY_CLEAN);
 
  protected :
  
   /**
   *    Please Document Me !
   *
   **/
   void printOn(ostream &output) const;

   /**
   *    Please Document Me !
   *
   **/
   void initDataOutput(FileIO *fileIO);

};


#endif //  __VLASOV_BLITZ_H


//#endif // GKC_CILK
