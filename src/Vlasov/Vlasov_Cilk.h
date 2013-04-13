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
   *    It is calculated as
   *    
   *    \f[
   *      \mathcal{N}{}_{\left< \chi \right> \times g_{1\sigma}} = 
   *        \frac{\partial \left< \chi \right>}{\partial x} \frac{\partial f_{1\sigma}}{\partial y}
   *      - \frac{\partial \left< \chi \right>}{\partial y} \frac{\partial f_{1\sigma}}{\partial x}
   *    \f]
   *
   *    Currently implemented using Arakawa (Morinishi) type scheme.
   *    
   *    @todo Describe and add reference
   *
   *
   *   Updated the CFL (Courant-Friedrich-Levy number). For explicit time-stepping
   *   the CFL value has to be always < 1 to ensure stability of the system
   *   (practically < 0.4).
   *   
   *   Note : Stability is still not guaranteed. As the system is unstable. Thus the
   *          time-stepping scheme needs to allows imaginary values e.g.
   *          (RK-3, RK-4, Heun method).
   *
   *   Calculated using ....
   *
   *   This needs only to be calculated in the non-linear terms
   *
   *  @note get linking error if defined inline. Check Performance !
   *
   *  @depreciated This is not done directly in the non-linear term
   *
   *  Conserved the L2 norm (L1 ? ), thus guarantees energy conservation.
   *  [ Probably particle conservation too .... !]
   *
   *
   **/
   void calculatePoissonBracket(const CComplex  G              [NzLB][Nky][NxLB  ][NvLB],   // in case of e-m
                                const CComplex Xi              [NzLB][Nky][NxLB+4][NvLB],   // in case of e-m
                                const CComplex  f [NsLD][NmLD ][NzLB][Nky][NxLB  ][NvLB],   // in case of e-s
                                const CComplex Fields [Nq][NsLD][NmLD ][NzLB][Nky][NxLB+4], // in case of e-s
                                const int z, const int m, const int s                     ,
                                CComplex ExB[Nky][NxLD][NvLD], double Xi_max[3], const bool electroMagnetic);
   /** 
   *
   *   @brief Calculate the parallel non-linearity as give by Goerler, PhD Thesis, Eq.(2.52)
   *
   *   @note implementation broken
   *
   *   \f[
   *      \mathcal{N}{}_{v_\parallel} = 
   *       \left\{ v_\parallel \hat{b}_0 \left( q_\sigma \nabla \left< \phi \right> 
   *     + \frac{q_\sigma}{c} \left< \dot{A}_{1\parallel} \right> \hat{b}_0
   *     + \mu \nabla \left< B_{1\parallel} \right> \right)
   *     + \frac{B_0}{B_{0\parallel}^\star} \left( v_\zeta + v_{\nabla B} + v_c \right) \cdot
   *          \left( q_\sigma \nabla \left< \phi \right> 
   *     + \mu \nabla \left( B_0 + \left< B_{1\parallel}\right> \right) \right) \right\}
   *          \frac{1}{m_\sigma v_\parallel} \frac{\partial f_{1\sigma}}{\partial v_\parallel}
   *   \f]
   *
   **/
   virtual void calculateParallelNonLinearity(
                                const CComplex f          [NsLD][NmLD][NzLB][Nky][NxLB   ][NvLB],
                                const CComplex Fields [Nq][NsLD][NmLD ][NzLB][Nky][NxLB+4], // in case of e-s
                                const int z, const int m, const int s                     ,
                                CComplex NonLinearTerm[Nky][NxLD][NvLD]);

   /**
   *
   *   @brief Calculated the gyro-average modified potential
   *
   *   \f[ 
   *        \left<\chi \right> = \left<\phi_1 \right> - \frac{v_\parallel}{c}\eft<{A}_{1\parallel} \right>
   *                + \frac{\mu}{q_\sigma} \left< {B}_{1\parallel} \right>
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
   *     f_{1\sigma} = g_{1\sigma} - \frac{q_\sigma}{c}\bar{A}_{1\parallel}\frac{v_\parallel}{T_{0\sigma}} F_{0\sigma} 
   *  \f]
   *
   *   Reference : @cite Note from Goerler's PhD Thesis (Eq. 2.32)
   *
   *  Note : In [ Xi, G] in electro-static simulations reduces to [ phi, f1 + sigma phi f_0 ],
   *        equals [phi, f1] + [phi, sigma phi f_0], as f_0 does not depend on y or x 
   *        [ Here, we do not use local assumption with partial_x f_0 = eta ... 
   *          but indeed use only numerically value ]. This changes of course once we go 
   *          to global version 
   *
   **/
   virtual void setupXiAndG(
                           const CComplex g          [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f0         [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Fields [Nq][NsLD][NmLD][NzLB][Nky][NxLB+4],
                                 CComplex Xi                     [NzLB][Nky][NxLB+4][NvLB],
                                 CComplex G                      [NzLB][Nky][NxLB  ][NvLB],
                           const int m, const int s);


   /**
   *    Please Document Me !
   *
   **/
   void Vlasov_EM(
                           const CComplex fs         [NsLD][NmLD][NzLB][Nky][NxLB   ][NvLB],
                           CComplex fss              [NsLD][NmLD][NzLB][Nky][NxLB   ][NvLB],
                           const CComplex f0         [NsLD][NmLD][NzLB][Nky][NxLB   ][NvLB],
                           const CComplex f1         [NsLD][NmLD][NzLB][Nky][NxLB   ][NvLB],
                           CComplex ft               [NsLD][NmLD][NzLB][Nky][NxLB   ][NvLB],
                           CComplex Coll             [NsLD][NmLD][NzLB][Nky][NxLB   ][NvLB],
                           const CComplex Fields [Nq][NsLD][NmLD ][NzLB][Nky][NxLB+4]      ,
                           CComplex Xi                            [NzLD][Nky][NxLB+4][NvLD],
                           CComplex G                             [NzLD][Nky][NxLD  ][NvLD],
                           CComplex ExB                                 [Nky][NxLB  ][NvLB],
                           const double Kx[NzLD][NxLD], const double Ky[NzLD][NxLD], 
                           const double dB_dz[NzLD][NxLD], // Geometry stuff
                           const double dt, const int rk_step, const double rk[3]);

   /**
   *
   *  Set the Krook operator 
   *  \f[
   *     \frac{\partial g_{1\sigma}}{\partial t} = \dots - \nu(x) g_{1\sigma}
   *  \f]
   * 
   *   Is used to damp oscillations close to the simulation boundary.
   *
   *  @note   * Is this the no-slip boundary condition ?
   *          * Violates conservation of particles, energy and momentum and
   *            needs to be fixed by modifying the fields. See Lapillone.
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

