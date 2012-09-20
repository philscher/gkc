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

#ifndef COLLISIONS_H
#define COLLISIONS_H

#include "Global.h"
#include "Setup.h"
#include "Geometry.h"
#include "Grid.h"

/**
 * 
 *  @brief implements various function which are defined in 
 *         PhD Thesis of Merz
 *
 */
class Collisions {
    
public:   

   /**
   *  @brief Constructor
   **/
   Collisions(Grid *_grid, Parallel *_parallel, Setup *_setup, FileIO *fileIO, Geometry *_geo, FFTSolver *(_fft)) 
   {
   
   };

   /**
    *
    *
   **/
   virtual int solve(std::string equation_type, Fields *fields, Array6C fs, Array6C fss, double dt, int rk_step) 
   {
            // we have collisionless system
            return 0;

   };

   /**
   *  @brief Derivative of the error function
   *  @image html Deriv_ErrorFunction.png
   *  \f[ 
   *      Derf(x) = \frac{d erf(x)}{dx} = \frac{2}{\sqrt{\pi}} \exp\left(-x^2 \right)  
   *  \f]
   **/
   __declspec(vector) inline double Derf(const double x) { return 2./sqrt(M_PI)*exp(-pow2(x)); };


   /**
   *  @brief Calculates helper function F1
   *  @image html Collision_F1.png
   *  \f[ F_1(x) = x \frac{d erf(x)}{dx} + \left( 2 x^2 - \right) erf(x) \f]
   **/
   inline double F1(const double x) { return x*Derf(x) + (2. * pow2(x) - 1.) * erf(x); };
   
   /**
   *  @brief Calculates helper function F2
   *  @image html Collision_F2.png
   *  \f[ F_2(x) = \left( 1 - \frac{2}{3} x^2 \right) erf(x) - x \frac{d erf}{dx} \f]
   **/
   inline double F2(const double x) { return (1.-2./3.*pow2(x)) * erf(x) - x* Derf(x); };

   /**
   *  @brief Calculates helper function F3
   *  @image html Collision_F3.png
   *  \f[ F_2(x) = F_1(x) + 3 F_2(x) \f]
   **/
   inline double F3(const double x) { return F1(x)+3.*F2(x); };

    
   /**
   *   @brief The Chandrasekhar function
   *
   *   \f[ Chandra(x) = \frac{erf(x) - x erf'(x)}{2 x^2} \f]
   *   @image html Chandrasekhar.png
   *
   **/
   inline double Chandra(const double x) {
   
     return (erf(x)- x* Derf(x))/(2.*pow2(x));

   }
   
};




#endif // COLLISIONS_H
