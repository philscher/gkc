/*
 * =====================================================================================
 *
 *       Filename: Vlasov_Cilk.h
 *
 *    Description: Implementation of GK Vlasov equation using 
 *                 Intel Cilk (Array Notation)
 *
 *         Author: Paul P. Hilscher (2011-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __VLASOV_CILK_H
#define __VLASOV_CILK_H

#include "Vlasov.h"

class Event;

/**
*  @brief Implementation of the Vlasov equation 
*
*  @section Calculating the non-linearity (Non-Linearity)
*
*  The Poisson brackets are calculated using the Mornishi scheme 
*  (see e.g. \cite{Morinishi_1998:FDScheme}, \cite{Arakawa_1966:DesignLong},
*  \cite{Idomura_2007:NewConservative}).
*
*  In the discretized equation system using CD-4 discretization, e.g. $\chi_{x;x,y}$ is given by
*  
*  \f[
*      \left[\frac{\partial \chi}{\partial x}\right]_{x,y} = \chi_{y;x,y} 
*     = \frac{8\left(\chi_{x+1,y} - \chi_{x-1,y}\right) - \left(\chi_{x+2,y} - \chi_{x-2,y} \right)}{12 \Delta x} \quad,
*  \f]
*   
*  with the usual definition of an equidistant grid discretization, e.g. $\Delta x = L_x / N_x$.
*  The Poisson bracket for \f$\mathcal{N}_{\chi \times g} \f$ is calculated by 
*   
*   \f{align}{
*    \Xi_{yx} =8
*    &\left( \left( \chi_{y;x,y} + \chi_{y;x+1,y} \right) g_{y,x+1} - 
*          \left( \chi_{y;x,y} + \chi_{y;x-1,y} \right) g_{y,x-1} \right) \\
*    - &
*          \left( \left( \chi_{y;x,y} + \chi_{y;x+2,y} \right) g_{y,x+2} - 
*          \left( \chi_{y;x,y}        + \chi_{y;x-2,y} \right) g_{y,x-2} \right) 
*    \f}
*
*   \f{align}{
*    \Xi_{xy} =8
*      &\left( \left( \chi_{y;x,y} + \chi_{y+1;x,y} \right) g_{y+1,x} - 
*       \left( \chi_{y;x,y} + \chi_{y-1;x,y} \right) g_{y-1,x} \right) \\
*     - &
*       \left( \left( \chi_{y;x,y} + \chi_{y+2;x,y} \right) g_{y+2,x} - 
*       \left( \chi_{y;x,y}        + \chi_{y-2;x,y} \right) g_{y-2,x} \right) 
*  \f}
*  This results
*  \f[
*    \mathcal{N}_{\chi \times g} = \frac{\Xi_{xy}}{24 \Delta x} + \frac{\Xi_{yz}}{24 \Delta y} \quad.
*  \f]
*
**/
class VlasovCilk : public Vlasov 
{

  friend class Event;

 protected:     
   
  
  /**
  *  @brief Calculates the \f$ \left[\chi, \g \right] \f$ non-linear term using Arakawa type scheme
  *
  *  The non-linear term correspond to the E x B drift. It is calculated as
  *    
  *  \f[
  *     \mathcal{N}{}_{\left< \chi \right> \times g_{1\sigma}} = 
  *     \frac{\partial \left< \chi \right>}{\partial x} \frac{\partial f_{1\sigma}}{\partial y}
  *   - \frac{\partial \left< \chi \right>}{\partial y} \frac{\partial f_{1\sigma}}{\partial x}
  *  \f]
  *
  *  Currently implemented using Arakawa (Morinishi) type scheme.
  *    
  *  @note Updates \f$ max(chi)_{x,yz} \f$ the for calculating the CFL for non-linear variable time stepping
  *
  **/
  void calculateExBNonLinearity(const CComplex  G              [NzLB][Nky][NxLB  ][NvLB],   // in case of e-m
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
  *        \left<\chi \right> = \left<\phi_1 \right> - \frac{v_\parallel}{c}\left<{A}_{1\parallel}\right>
  *                + \frac{\mu}{q_\sigma} \left< {B}_{1\parallel} \right>
  *   \f]
  *
  *  as well as the Gamma abbreviation functions (Eq. 2.50 below)
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
  *  Note : In \f$ \left[\Xi, G \right]\f$ in electro-static simulations reduces to 
  *        \f$ \left[ \phi, f_1 + \sigma \phi f_0 \right] \f$, equals 
  *        \f$ \left[\phi, f_1 \right] + \left[\phi, \sigma \phi f_0 \right] \f$, as \f$ f_0 \f$ does not depend 
  *         on \f$ y \f$ or \f$ x \f$ 
  *        [ Here, we do not use local assumption with \f$ partial_x f_0 = \eta ... $\f 
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

#endif //  __VLASOV_CILK_H

