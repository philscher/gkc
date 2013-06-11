/*
 * =====================================================================================
 *
 *       Filename: Collision_HyperDiffusion.h
 *
 *    Description: Hyper-Diffusion collisional operator to test collapse
 *                 of CvK modes
 *
 *         Author: Paul P. Hilscher (2013-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __GKC_COLLISION_LENARD_DIFFUSION_H__
#define __GKC_COLLISION_LENARD_DIFFUSION_H__

#include "Collisions/Collisions.h"

/**
* 
*  @brief Hyper-diffusion collisional operator
*
*  The collisional operator are given by
*
*  \f[
*     \mathcal{C}_{\textrm{hyp}, n} = \beta_c\; \partial_{v_\parallel}^n f_{1} \quad,
*  \f]
*  where $\beta_c$ is the collisionality and $n$ is the order of hyper-diffusivity
*  ( here up to 8-th order is implemented) using either second order or fourth order
*  finite difference derivatives.
* 
*
*
**/
class Collisions_HyperDiffusion : public Collisions {
    
 protected:   
    
  double beta[SPECIES_MAX+1];         ///< Collisionality 

  int hyp_order;
 public:

  /**
  *
  *   @brief constructor
  *
  *   accepts following setup parameters
  *
  *
  **/
  Collisions_HyperDiffusion(Grid *grid, Parallel *parallel, Setup *setup, FileIO *fileIO, Geometry *geo); 
       
  /**
  *
  *
  **/ 
 ~Collisions_HyperDiffusion();
 
  /**
  *   Calculate Collisional corrections
  *
  *
  **/
  void solve(Fields *fields, const CComplex  *f, const CComplex *f0, CComplex *Coll, double dt, int rk_step); 

 protected:

  /**
  *   Set Data output parameters
  *
  *
  **/
  virtual void initData(Setup *setup, FileIO *fileIO);
  
  /**
  * Program output
  *
  *
  **/
  virtual void printOn(std::ostream &output) const;

};

#endif // __GKC_COLLISION_HYPER_DIFFUSION_H__
