/*
 * =====================================================================================
 *
 *       Filename: RootFinding.h
 *
 *    Description: Implements various one dimensional root finding algrotihm.
 *
 *         Author: Paul P. Hilscher (2013), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __GKC_ROOTFINDING_H__
#define __GKC_ROOTFINDING_H__

#include <functional>
#include <cassert>


/**
*
*   @brief Root finding for 1-dimensional functions  
*
*   Implements severl one dimensional root-finding algortihm
*
*
**/
class RootFinding {

 public:

  /**
  *   @brief Bisection method
  *
  *   Reference : http://en.wikipedia.org/wiki/Bisection_method
  *
  *   @param x0 
  *
  *   @note f(x0) and f(x1) have to be different sign
  **/
  static double BiSection(std::function<double (double)> func, double x_a, double x_b, int maxIter=512, double ftol=1.e-9)
  {

    // return True for +, False for - (doesn't handle nan, but -0,+0)
    auto sign = [] (double val) -> bool { return val == std::abs(val); };
     
    // evaluate functions
    double f_a = func(x_a), 
           f_b = func(x_b);

    if(f_a == 0.) return x_a;
    if(f_b == 0.) return x_b;
    
    // Bi-section method needs sign change
    assert(f_a * f_b < 0.);

    double x_m; // mid-point 

    for(int iter = 0; iter < maxIter; iter++) {
        
      // update midpoint
      x_m        = (x_a + x_b)/2.;
      double f_m = func(x_m);
     
      if(std::abs(f_m) < ftol) break;

      f_a = func(x_a);
      
      (sign(f_m) == sign(f_a)) ? x_a = x_m : x_b = x_m;
    }

    return x_m;
  };

  /**
  *   @brief Secant method
  *
  *   Is similar to newtons method but derivative is approximated using
  *   finite differences.
  *
  *   Reference : 
  *
  *   @param x0 and x1 should idealliy be close to the root 
  *
  **/
  static double Secant(std::function<double (double)> func, double x0, double x1, int maxIter=512, double ftol=1.e-9)
  {

    if(func(x0) == 0.) return x0;
    if(func(x1) == 0.) return x1;
   

    double x_nm1 = x0, f_nm1 = func(x_nm1), 
           x_nm2 = x1, f_nm2 = func(x_nm2);
  
    // iterate and update
    for(int iter = 0; iter < maxIter; iter++) {

      double x_n = x_nm1 - func(x_nm1) * (x_nm1 - x_nm2) / (f_nm1 - f_nm2);

      // update
      x_nm2 = x_nm1;
      x_nm1 = x_n  ;

      f_nm2 = f_nm1;
      f_nm1 = func(x_n);

      // break if tolerance is reached
      if(std::abs(f_nm1) <= ftol) break;
    }

    return x_nm1;
  };
  
};
  
#endif // __GKC_ROOTFINDING_H__
