/*
 * =====================================================================================
 *
 *       Filename: Collisions.h
 *
 *    Description: Defines helper functions for Collisionality
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __GKC_COLLISIONS_H__
#define __GKC_COLLISIONS_H__

#include "Global.h"

#include "Setup.h"
#include "Geometry.h"
#include "Grid.h"
#include "Fields.h"

class Vlasov;

/**
 * 
 *  @brief Base class for collisions
 *
 *  Thus functions does merly define some popular functions. 
 *  As well as the collisionless operator.
 *
 *
 */
class Collisions : public IfaceGKC {
  
protected:

  Grid *grid;
  Parallel *parallel;

public:   

  /**
  *  @brief Constructor
  *
  **/
  Collisions(Grid *_grid, Parallel *_parallel, Setup *_setup, FileIO *fileIO, Geometry *_geo)
     : grid(_grid), parallel(_parallel)
  {
   
  };

  /**
  *
  *
  **/
  virtual void solve(Fields *fields, const CComplex  *fs, const CComplex *f0, CComplex *Coll, double dt, int rk_step) 
  {
            // we have collisionless system
            return;

  };

  /**
  *  @brief Derivative of the error function
  *  @image html Deriv_ErrorFunction.png
  *  \f[ 
  *      Derf(x) = \frac{d erf(x)}{dx} = \frac{2}{\sqrt{\pi}} \exp\left(-x^2 \right)  
  *  \f]
  **/
  __declspec(vector) static inline double Derf(const double x) { return 2./sqrt(M_PI)*exp(-pow2(x)); };


  /**
  *  @brief Calculates helper function F1
  *  @image html Collision_F1.png
  *
  *  \f[ 
  *     F_1(x) = x \frac{d erf(x)}{dx} + \left( 2 x^2 - \right) erf(x) 
  *  \f]
  *
  **/
  static inline double F1(const double x) { return x*Derf(x) + (2. * pow2(x) - 1.) * erf(x); };
   
  /**
  *  @brief Calculates helper function F2
  *  @image html Collision_F2.png
  *
  *  \f[ 
  *     F_2(x) = \left( 1 - \frac{2}{3} x^2 \right) erf(x) - x \frac{d erf}{dx} 
  *  \f]
  *
  **/
  static inline double F2(const double x) { return (1.-2./3.*pow2(x)) * erf(x) - x* Derf(x); };

  /**
  *  @brief Calculates helper function F3
  *  @image html Collision_F3.png
  *
  *  \f[ 
  *      F_2(x) = F_1(x) + 3 F_2(x) 
  *  \f]
  **/
  static inline double F3(const double x) { return F1(x)+3.*F2(x); };

    
  /**
  *   @brief The Chandrasekhar function
  *   @image html Chandrasekhar.png
  *
  *   \f[ 
  *       Chandra(x) = \frac{erf(x) - x erf'(x)}{2 x^2} 
  *   \f]
  *
  *
  **/
  static inline double Chandra(const double x) {
   
     return (erf(x)- x* Derf(x))/(2.*pow2(x));

  }

  /**
  *     Document me please
  *
  **/ 
  virtual void printOn(std::ostream &output) const {};
  
  /**
  *     Document me please
  *
  **/ 
  virtual void initData(Setup *setup, FileIO *fileIO) {
  
    hid_t collisionGroup = fileIO->newGroup("Collisions");
     
    check(H5LTset_attribute_string(collisionGroup, ".", "Model", "Collisionless"), DMESG("H5LTset_attribute"));
            
    H5Gclose(collisionGroup);


};


#endif // __GKC_COLLISIONS_H__
