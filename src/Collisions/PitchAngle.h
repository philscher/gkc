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
*  @brief Pitch-Angle scattering 
*
*  Warning this is package is broken and does nothing yet.
*
*   Reference:
*
*     PhD Thesis of Merz
*
*
*/
class Collisions_PitchAngle : public Collisions {
    
    double nu; ///< Collisional frequency 

public:   
    
  Collisions_PitchAngle(Grid *grid, Parallel *parallel, Setup *setup, FileIO *fileIO, Geometry *geo);

  /**
  *    \f[ \left< C_{\sigma j}^{\perp} \right> = - \frac{1}{v^5} \left( 2 v^2 F_1 + 3 B_0 \mu F_2  \right) \nabla_\perp f_\sigma \f]
  *
  **/
  void C_xy(const double df_dp, const int x, const int y, const int z, const int v, const int m, const int s);
        
  /**
  *    \f[ \left< C_{\sigma j}^{\perp} \right> = - \frac{1}{v^5} \left( 2 v^2 F_1 + 3 B_0 \mu F_2  \right) \nabla_\perp f_\sigma \f]
  *
  **/
  void C_vm(const double df_dv, const double df_dm, const double f, const int x, const int y, const int z, const int v, const int m, const int s);
  
  /**
  *   Calculate Collisional corrections
  *
  *
  **/
  void solve(Fields *fields, const CComplex  *f, const CComplex *f0, CComplex *Coll, double dt, int rk_step); 

  
  virtual void initData(hid_t fileID);
     
 protected:

  virtual void printOn(std::ostream &output) const;
};

#endif // COLLISIONS__PITCHANGLE_H
