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
#include "Global.h"


#ifndef __VLASOV_CILK_H
#define __VLASOV_CILK_H

#include "Vlasov.h"


class Event;

/**
*  @brief Implementation of Vlasov's equation using Cilk Array Notation
*
*
*
**/
class VlasovCilk : public Vlasov {


  friend class Event;

  protected:     
   
  
   /**
   *    @brief Calculates the non-linear term using Arakawa type scheme
   *
   *    The non-linear term correspond to the E x B drift.
   *    It is caculated as
   *    
   *    \f[
   *        \left[< Xi >, f_{1sigma} \right] = frac{partial Xi}
   *    \f]
   *
   *    Currently implemented as Arakawa (Morinishi) type scheme.
   *    
   *    @todo Describe and add reference
   *
   *
   *   Updated the CFL (Courant-Friedrich-Levy number). For explicit time-stepening
   *   the CFL value has to be always < 1 to ensure stability of the system
   *   (practically < 0.4).
   *   
   *   Note : Stability is still not guranteed. As the system is unstable. Thus the
   *          time-steppening scheme needs to allows imaginary values e.g.
   *          (RK-3, RK-4, Heun method).
   *
   *   Calculated using ....
   *
   *   This needs only to be caluclated in the non-linear terms
   *
   *  @note get linking error if defined inline. Check Performance !
   *
   *  @depracated This is not done direclty in the non-linear term
   *
   *
   **/
   void calculatePoissonBracket(const CComplex  G              [NzLB][NkyLD][NxLB  ][NvLB],   // in case of e-m
                                const CComplex Xi              [NzLB][NkyLD][NxLB+4][NvLB],   // in case of e-m
                                const CComplex  f [NsLD][NmLD ][NzLB][NkyLD][NxLB  ][NvLB],   // in case of e-s
                                const CComplex Fields [Nq][NsLD][NmLD ][NzLB][NkyLD][NxLB+4], // in case of e-s
                                const int z, const int m, const int s                     ,
                                CComplex ExB[NkyLD][NxLD][NvLD], double Xi_max[3], const bool electroMagnetic);


   /**
   *
   *   @brief Calculated the gyro-average modified potential
   *
   *   \f[ 
   *        \bar{\Xi} = \bar{\phi}_1 - \frac{v_\parallel}{c}\bar{A}_{1\parallel} 
   *                + \frac{\mu}{q_\sigma} \bar{B}_{1\parallel}
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
                           const CComplex g           [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex f0         [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex Fields [Nq][NsLD][NmLD ][NzLB][NkyLD][NxLB+4],
                           CComplex Xi                         [NzLB][NkyLD][NxLB  ][NvLB],
                           CComplex G                          [NzLB][NkyLD][NxLB  ][NvLB],
                           const double V[NvGB], const double M[NmGB],
                           const int m, const int s);


   /**
   *    Please Document Me !
   *
   **/
   void Vlasov_EM(
                           const CComplex fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           CComplex fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex f0 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           CComplex ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           CComplex Coll     [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex Fields [Nq][NsLD][NmLD ][NzLB][NkyLD][NxLB+4],
                           CComplex Xi       [NzLD][NkyLD][NxLD  ][NvLD],
                           CComplex G        [NzLD][NkyLD][NxLD  ][NvLD],
                           CComplex ExB            [NkyLD][NxLB  ][NvLB],
                           const double Kx[NzLD][NxLD], const double Ky[NzLD][NxLD], const double dB_dz[NzLD][NxLD], // Geometry stuff
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           const double dt, const int rk_step, const double rk[3]);

   /**
   *
   *  Set the Krook operator 
   *  \f[
   *     \frac{\partial g_{1sigma}}{\partial t} = \dots - \nu(x) g_{1\sigma}
   *  \f]
   * Is used to damp oscillations close to the simulation boundary.
   *
   *  Note : 
   *          * Is this the no-slip boundary condition ?
   *          * Violates conservation of particles, energy and momentum and
   *            needs to be fixed by modifing the fields. See Lapillone.
   *
   *
   **/
   double krook_nu;

  public:

   /**
   *    Please Document Me !
   *
   **/
   VlasovCilk(Grid *_grid, Parallel *_parallel, Setup *_setup, FileIO * fileIO, Geometry *_geo, FFTSolver *fft, Benchmark *bench, Collisions *coll); 
        
   /**
   *    Please Document Me !
   *
   **/
   virtual void solve(std::string equation_tyoe, Fields *fields, CComplex *fs, CComplex *fss,
                      double dt, int rk_step, const double rk[3]);
   
 
  protected :
  
   /**
   *    Please Document Me !
   *
   **/
   void printOn(std::ostream &output) const;

   /**
   *    Please Document Me !
   *
   **/
   void initData(Setup *setup, FileIO *fileIO);

};


#endif //  __VLASOV_BLITZ_H


//#endif // GKC_CILK
