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

	/* Helper function for collisional operator as defined in 
	*  PhD Thesis of Merz
	*/
class Collisions {
    
    /*
 Collisional frequency 
*/
public:   

     Collisions(Grid *_grid, Parallel *_parallel, Setup *_setup, FileIO *fileIO, Geometry<HELIOS_GEOMETRY> *_geo, FFTSolver *(_fft)) 
     {


     };

     virtual int solve(std::string equation_type, Fields *fields, Array6z fs, Array6z fss, double dt, int rk_step) 
     {
            // we have collisionless system
            return 0;

     };

    /**
	 * Derivative of the error function 
     * By definition
     *  \f[ Derf(x) = \frac{d erf(x)}{dx} = \frac{2}{\sqrt{\pi}} \exp\left(-x^2 \right)  \f]
     * */
	inline double Derf(const double x) { return 2./sqrt(M_PI)*exp(-pow2(x)); };


    /**
     *  Calculates helper functions
     *  \f[ F_1(x) = x \frac{d erf(x)}{dx} + \left( 2 x^2 - \right) erf(x) \f]
     * */
	inline double F1(const double x) { return x*Derf(x) + (2. * pow2(x) - 1.) * erf(x); };
    /**
     *  Calculates helper functions
     *  \f[ F_2(x) = \left( 1 - \frac{2}{3} x^2 \right) erf(x) - x \frac{d erf}{dx} \f]
     * */
	inline double F2(const double x) { return (1.-2./3.*pow2(x)) * erf(x) - x* Derf(x); };
    /**
     *  Calculates helper functions
     *  \f[ F_2(x) = F_1(x) + 3 F_2(x) \f]
     * */
	inline double F3(const double x) { return F1(x)+3.*F2(x); };

    
/**
 *   The Chandrasekhar function
 *
 *   \f[ Chandra(x) = \frac{erf(x) - x erf'(x)}{2 x^2} \f]
 *
 * */
      inline double Chandra(const double x) {
            return (erf(x)- x* Derf(x))/(2.*pow2(x));


      }
   
};




#endif // COLLISIONS_H
