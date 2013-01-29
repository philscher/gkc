/*
 * =====================================================================================
 *
 *       Filename: Moments.h
 *
 *    Description: Calculation of the Moments of the Phase-space function
 *
 *         Author: Paul P. Hilscher (2012-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __MOMENTS_H__
#define __MOMENTS_H__

#include<cmath>

#include "Global.h"
#include "Setup.h"

#include "Grid.h"
#include "Vlasov/Vlasov.h"
#include "Fields/Fields.h"


/**
*    @brief Calculate Moments of the Vlasov equations
*
*    Calculates various moments of the Vlasov equation. 
*
*    Based on Görler, PhD Thesis, Eq.(2.73), \cite{GoerlerPhD}
*
*    Generally the moments of the Vlasov equation can be calculated
*    using 
*
*    \f{multline}{
*        M_{ab,\sigma}(x) = n_{ref} \hat{n}_{0\sigma}(x_0) 
*                          c_{ref}^{a+b} v_{T\sigma}^{a+b}(x_0) \frac{\rho_{ref}}{L_{ref}}
*           \left\{
*              \pi B_0^{b/2} \int \int B_{0\parallel}^\star 
*                   \left< f_{1\sigma} \right> v_\parallel^a \mu^{b/2} dv_\parallel d\mu
*              -   \frac{n_{p\sigma}}{T_{0\sigma}} T_{p\sigma}^{(a+b)/2}
*                  \left[ \Upsilon(a) 
*                  + \beta_{ref} \frac{T_{0\sigma}}{B_0^2} 
*                  \frac{j_{0\parallel}}{q_\sigma v_{T\sigma}} \Upsilon(a+1)\right]i
*                    \right. \\ \quad \quad \quad \left.
*                  \left((b/2)! q_\sigma \phi - \left( \frac{B_0}{T_{p\sigma}}\right)^{b/2+1}
*                  \int \left(  q_\sigma  \left< \left< \phi \right> \right> 
*                  + T_{0\sigma}(x_0) \mu \left< \left< B_{1\parallel}\right> \right> \right)
*                    e^\frac{- \mu B_0}{T_{\sigma}} \mu^{b/2} d\mu
*                  \right)
*           \right\}
*    \f}
*
*    From Go\"rler, PhD, Eq.(2.72).
*   
*    \f[
*         \Upsilon(a) = \frac{1}{\sqrt{\pi}} \int_\infty^\infty x^a e^{-x^2} dx
*
*         = \left\{ \begin{array}{lll} 0, & & a \textrm{odd}  \\
*                                     1,  & & a \textrm{even} \\
*                                     \frac{1 \cdot 3 \dots (a-1)}{\sqrt{2}^a}
*                  \end{array}  \right.
*    \f]
*
*    where for \f$ k_\perp^2 > 1 \f$, additional correction terms are included.
*
*   The perturbed parallel temperature is defined as, Goerler, PhD Thesis, Eq. (3.81)
*   \f[
*        \hat{T}_{parallel1,sigma} = 
*            \frac{2 M_{20} - T_{p\sigma} M_{00}}{\hat{n}_{p\sigma}}
*   \f]
*  
*  and the perpendicular temperature as, Eq.(3.82)
*   
*   \f[
*        \hat{T}_{\perp,1\sigma} = 
*            \frac{M_{02} - T_{p\sigma} M_{00}}{\hat{n}_{p\sigma}}
*   \f]
*
*   where, the second term on the RHS arises, as we substract the constribution from the
*   equilibrium temperature.
*
**/
class Moments 
{

  Vlasov   *vlasov  ;
  Fields   *fields  ;
  Grid     *grid    ;
  Parallel *parallel;
   
  bool doFieldCorrections; ///< Include gyro-kinetic field corrections

  /**
  *    @brief 
  *  
  *    From T. Go\"rler, PhD thesis, Eq.(2.72),
  *   
  *    \f[
  *         \Upsilon(a) = \frac{1}{\sqrt{\pi}} \int_\infty^\infty x^a e^{-x^2} dx
  * 
  *         = \left\{ \begin{array}{lll} 0, & & a \textrm{odd}  \\
  *                                     1,  & & a \textrm{even} \\
  *                                     \frac{1 \cdot 3 \dots (a-1)}{\sqrt{2}^a}
  *                  \end{array}  \right.
  *    \f]
  *
  *    which can be simplified to, X. Lapillone (PhD thesis) Eq.(2.78),
  *
  *    \f[
  *         \Upsilon(a) =
  *         = \left\{ \begin{array}{lll} 0, & & a \textrm{odd}  \\
  *              \frac{\Gamma \left(\frac{m+1}{2} \right)}{\sqrt{2}} 
  *                                         & & a \textrm{even} 
  *                  \end{array}  \right.
  *    \f]
  *
  **/
  double Y(int a)  { return (a % 2) ? 0.
                                    : std::tgamma((a+1.)/2.)/sqrt(M_PI); };

 public:  

  /**
  *    @brief Class for calculating moments
  *
  *    setup : includeFieldCorrections
  *
  *
  **/
  Moments(Setup *setup, Vlasov *vlasov, Fields *fields, Grid *grid, Parallel *parallel); 

  /**
  *   @brief Calculates specific moments of the distribution function
  *
  **/
  void getMoment(const CComplex     f    [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                 const CComplex Field0[Nq][NzLD][Nky][NxLD],
                 CComplex Mom[8][NsLD][NzLD][Nky][NxLD], 
                 const int a, const int b, const int idx);

  /**
  *   @brief Calculates moments of the phase-space function
  *
  *   @param s species index
  *
  **/
  void getMoments(const CComplex     f    [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                  const CComplex   Field0[Nq][NzLD][Nky][NxLD],
                        CComplex Mom[8][NsLD][NzLD][Nky][NxLD]); 

};

#endif // __MOMENTS_H__
