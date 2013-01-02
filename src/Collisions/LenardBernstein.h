/*
 * =====================================================================================
 *
 *       Filename: Collision_LenardBernstein.h
 *
 *    Description: Lenard-Bernstein collisional operator including
 *                 energy & momentum conservational terms
 *
 *         Author: Paul P. Hilscher (2012-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __GKC_COLLISION_LENARD_BERNSTEIN_H__
#define __GKC_COLLISION_LENARD_BERNSTEIN_H__

#include "Collisions/Collisions.h"


/**
*   @brief Lenard-Bernstein collision term
*
*   The Lenard-Bernstein operator was descirbed in \cite{LenardBernstein_1958:DiffusionVP}.
*
*   The Lenard-Bernstein operator has the form of
*   \f[ 
*       \mathcal{C}_{LB} =  \beta \frac{\partial}{\partial v} \left(  v F_{1}
*                                 + v^2_0 \frac{\partial F_1}{\partial v} \right)
*   \f]
*
*   with \f$ v_0 \f$, the root mean square speed corresponding to the equilibrium distribution \f[ F_0 \f]
*  Quare of Root mean square defined as \f[ v_{rms}^2 = \int_0^{\infty} v^2 F_{0\sigma}  = 1/2\f] and
*  Folowing Lenard and Bernstein, \f[ \beta \f] can be roughly estimated by comapring with the true
*  Fokker-Planck equation. One get's approximately \f[ \beta \approx \frac{4 \ pi e^4}{m^2 v_0^3} \f].
*
*   Here, we include additionally correction terms to ensure particle, momentum and energy
*   conservation as described by \cite{Satake_2008:DevTranpCode}. The correction
*   terms are described by
*     
*
*
*  
*
*
**/
class Collisions_LenardBernstein : public Collisions {
    
 protected:   
    
  double beta;         ///< Collisionality 
  bool   consvMoment;  ///< Set if 0-2 Moments are conserved

  double *nu,   ///<            \f$ \nu = \left( v_\parallel^2 + 2 \mu \right) / v_{th}^2       \f$
         *a ,   ///< Pre-factor \f$ a   = 1 - \frac{\pi}{2}\left( erf - derf \right) \nu^{-1/2} \f$
         *b ,   ///< Pre-factor \f$ b   = v_\parallel x\nu^{-3/2} erf(x)                        \f$
         *c ;   ///< Pre-factor \f$ c   = \nu^{-1/2} \left(  erf(\nu) - erf'(\nu) \right)       \f$

  CComplex *dn, ///< \f$ \int f_{1\sigma} dv_\parallel d\mu             \f$
           *dP, ///< \f$ \int f_{1\sigma} v_\parallel dv_\parallel d\mu \f$
           *dE; ///< \f$ \int f_{1\sigma} \nu dv_\parallel d\mu         \f$
   
  nct::allocate ArrayPreFactors    , ///< Array class for dn, dP, dE
                ArrayCorrectionTerm; ///< Array class for nu, a, b, c

  /**
  *   @brief calculates the pre-terms
  *
  *   \f{align}{
  *         a &= 1 - 3 \sqrt{\pi}{2} \left( erf(x) - erf'(nu)  \right) * v^{-1/2} \\
  *         b &= v_\parallel x^{-3/2} erf(x)                                      \\
  *         c &= x^{-1/2} \left( erf(x) - erf'(x) \right)
  *   \f}
  *   
  *
  **/ 
  void calculatePreTerms(double a[NsLD][NmLD][NvLD], double  b[NsLD][NmLD][NvLD], 
                         double c[NsLD][NmLD][NvLD], double nu[NsLD][NmLD][NvLD]);

 public:

  
  /**
  *
  *   @brief constructor
  *
  *   accepts following setup parameters
  *
  *
  **/
  Collisions_LenardBernstein(Grid *grid, Parallel *parallel, Setup *setup, FileIO *fileIO, Geometry *geo); 
       
 
  /**
  *
  *
  **/ 
 ~Collisions_LenardBernstein();
 
  /**
  *   Set Data output parameters
  *
  *
  **/
  virtual void initData(hid_t fileID); 

  
  /**
  *   Calculate Collisional corrections
  *
  *
  **/
  void solve(Fields *fields, const CComplex  *f, const CComplex *f0, CComplex *Coll, double dt, int rk_step); 

 protected:

  
  /**
  * Program output
  *
  *
  **/
  virtual void printOn(std::ostream &output) const;


};




#endif // __GKC_COLLISION_LENARD_BERNSTEIN_H__
