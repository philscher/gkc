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

/**
*
*   @brief Hermite interpolation function
*
*   HN_nm - N is interpolation order, n is radial function,
*   m it's derivative [ 0, 1, 2, ..] -> [ no, first, second, ... ]
*   
*   Hermite polynomial interpolation functions of n-th order $H^n$, where
*   for x is
*   \f[
*       x \set )0,1( \eq 0
*   \f]
*
*   @todo add 9th order
*   @todo add reference
*   @bug  H_02 is not defined correctly 
*   
**/
class HermiteInterpolation
{
  
  public:

   /**
   *  @brief Basis functions for 3rd-order interpolation
   *  @image html  HermiteInterpolation_3rdOrderFunctions.png
   **/
   /// 3rd order basis function \f$ H^3_{00}(x) = 2 x^3 - 3 x^2 + 1 \f$
   static double H3_00( const double x) { return  ((x >= 0.) && (x <= 1.)) ?  2. * pow3(x) - 3. * pow2(x) + 1. : 0. ; } ;  
   /// 3rd order basis function \f$ H^3_{10}(x) = - 2 x^3 + 3 x^2     \f$
   static double H3_10( const double x) { return  ((x >= 0.) && (x <= 1.)) ? -2. * pow3(x) + 3. * pow2(x)      : 0. ; } ; 
   /// 3rd order basis function \f$ H^3_{01}(x) =     x^3 - 2 x^2 + x \f$
   static double H3_01( const double x) { return  ((x >= 0.) && (x <= 1.)) ?       pow3(x) - 2. * pow2(x) + x  : 0. ; } ; 
   /// 3rd order basis function \f$ H^3_{11}(x) =     x^3 -   x^2     \f$
   static double H3_11( const double x) { return  ((x >= 0.) && (x <= 1.)) ?       pow3(x) -      pow2(x)      : 0. ; } ; 
    
   /**
   *  @brief Basis functions for 5th-order interpolation
   *  @image html  HermiteInterpolation_5thOrderFunctions.png
   **/
   /// 5th order basis function \f$ H^5_{00} = -6x^5 + 15 x^4 - 10 x^3 + 1 \f$
   static double H5_00( const double x) { return ((x >= 0.) && (x <= 1.)) ? - 6.  * pow5(x) + 15.  * pow4(x) - 10.  * pow3(x) + 1.            : 0. ; };
   /// 5th order basis function \f$ H^5_{10} =  6x^5 - 15 x^4 + 10 x^3     \f$
   static double H5_10( const double x) { return ((x >= 0.) && (x <= 1.)) ?   6.  * pow5(x) - 15.  * pow4(x) + 10.  * pow3(x)                 : 0. ; };
   /// 5th order basis function \f$ H^5_{01} =  -3x^5 + 8 x^4 - 6 x^3 + x  \f$
   static double H5_01( const double x) { return ((x >= 0.) && (x <= 1.)) ? - 3.  * pow5(x) +  8.  * pow4(x) -  6.  * pow3(x) + x             : 0. ; };
   /// 5th order basis function \f$ H^5_{11} =  -3x^5 + 7 x^4 - 4 x^3      \f$
   static double H5_11( const double x) { return ((x >= 0.) && (x <= 1.)) ? - 3.  * pow5(x) +  7.  * pow4(x) -  4.  * pow3(x)                 : 0. ; };
   /// 5th order basis function \f$ H^5_{02} =  -\tfrac{1}{2}x^5 + \tfrac{3}{2} x^4 - \tfrac{3}{2} x^3 + \tfrac{1}{2} x^2      \f$ @bug not correct
   static double H5_02( const double x) { return ((x >= 0.) && (x <= 1.)) ? - 0.5 * pow5(x) +  1.5 * pow4(x) -  1.5 * pow3(x) + 0.5 * pow2(x) : 0. ; };
   /// 5th order basis function \f$ H^5_{12} =  -\tfrac{1}{2}x^5 -              x^4 + \tfrac{1}{2} x^3    \f$
   static double H5_12( const double x) { return ((x >= 0.) && (x <= 1.)) ?   0.5 * pow5(x) -        pow4(x) +  0.5 * pow3(x)                 : 0. ; };

    // add 9th order
};


#endif // HERMITE_INTERPOLATION_H
