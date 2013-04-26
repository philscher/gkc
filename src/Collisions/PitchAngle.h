/*
 * =====================================================================================
 *
 *       Filename: PitchAngle.h
 *
 *    Description: Implementation of the Pitch angle scattering 
 *
 *         Author: Paul P. Hilscher (2012), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef COLLISIONS__PITCHANGLE_H
#define COLLISIONS__PITCHANGLE_H

#include "Collisions/Collisions.h"

/* 
*
*
*   Helper function for collisional operator as defined in 
*
*   Reference:
*
*     PhD Thesis of Merz
*
*
*/
class Collisions_PitchAngle : Collision {
    
    double nu; ///< Collisional frequency 
public:   
    
  Collisions_PitchAngle(Setup *setup,  Grid *grid);

  /**
  *    \f[ \left< C_{\sigma j}^{\perp} \right> = - \frac{1}{v^5} \left( 2 v^2 F_1 + 3 B_0 \mu F_2  \right) \nabla_\perp f_\sigma \f]
  *
  **/
  inline double C_xy(const double df_dp, const int x, const int y, const int z, const int v, const int m, const int s) {
        
  /**
  *    \f[ \left< C_{\sigma j}^{\perp} \right> = - \frac{1}{v^5} \left( 2 v^2 F_1 + 3 B_0 \mu F_2  \right) \nabla_\perp f_\sigma \f]
  *
  **/
  inline double C_vm(const double df_dv, const double df_dm, const double f, const int x, const int y, const int z, const int v, const int m, const int s) {
  
  /**
  *   Calculate Collisional corrections
  *
  *
  **/
  void solve(Fields *fields, const CComplex  *f, const CComplex *f0, CComplex *Coll, double dt, int rk_step); 

  
  void initData(hid_t fileID);
     
 protected:

   virtual void printOn(std::ostream &output) const;
};

#endif // COLLISIONS__PITCHANGLE_H
