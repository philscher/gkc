/*
 * =====================================================================================
 *
 *       Filename: HermiteInterpolation.cpp
 *
 *    Description: Inlcudes definitions for Hermite Interpolations
 *
 *         Author: Paul P. Hilscher (2012), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#include "config.h"
#include "Global.h"

#ifndef HERMITE_INTERPOLATION_H
#define HERMITE_INTERPOLATION_H

    // needs to be defined between 0 and 1
/**
 *
 *
 *
 *   HN_nm - N is interpolation order, n is radial function,
 *   m it's derivative [ 0, 1, 2, ..] -> [ no, first, second, ... ]
 *
 *
 *
 **/
class HermiteInterpolation
{
 //   static double sign(double x) { return ( x >= 0.) ? 1. : -1 ; }; 

  public:

   

    // Basis functions for 3rd-order interpolation
    static double H3_00( const double x) { return  ((x >= 0.) && (x <= 1.)) ?  2. * pow3(x) - 3. * pow2(x) + 1. : 0. ; } ;  
    static double H3_10( const double x) { return  ((x >= 0.) && (x <= 1.)) ? -2. * pow3(x) + 3. * pow2(x)      : 0. ; } ; 
    static double H3_01( const double x) { return  ((x >= 0.) && (x <= 1.)) ?       pow3(x) - 2. * pow2(x) + x  : 0. ; } ; 
    static double H3_11( const double x) { return  ((x >= 0.) && (x <= 1.)) ?       pow3(x) -      pow2(x)      : 0. ; } ; 
    
    // Basis functions for 5th-order interpolation
    static double H5_00( const double x) { return ((x >= 0.) && (x <= 1.)) ? - 6.  * pow5(x) + 15.  * pow4(x) - 10.  * pow3(x) + 1.            : 0. ; };
    static double H5_10( const double x) { return ((x >= 0.) && (x <= 1.)) ?   6.  * pow5(x) - 15.  * pow4(x) + 10.  * pow3(x)                 : 0. ; };
    static double H5_01( const double x) { return ((x >= 0.) && (x <= 1.)) ? - 3.  * pow5(x) +  8.  * pow4(x) -  6.  * pow3(x) + x             : 0. ; };
    static double H5_11( const double x) { return ((x >= 0.) && (x <= 1.)) ? - 3.  * pow5(x) +  7.  * pow4(x) -  4.  * pow3(x)                 : 0. ; };
    static double H5_02( const double x) { return ((x >= 0.) && (x <= 1.)) ? - 0.5 * pow5(x) +  1.5 * pow4(x) -  1.5 * pow3(x) + 0.5 * pow2(x) : 0. ; };
    static double H5_12( const double x) { return ((x >= 0.) && (x <= 1.)) ?   0.5 * pow5(x) -        pow4(x) +  0.5 * pow3(x)                 : 0. ; };

    // more ?
};


#endif // HERMITE_INTERPOLATION_H
