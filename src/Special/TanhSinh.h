/*
 * =====================================================================================
 *
 *       Filename: TanhSinh.cpp
 *
 *    Description: Tanh-Sinh integration scheme for one-dimensional functions
 *
 *         Author: Paul P. Hilscher (2009-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __TANHSINH_H__
#define __TANHSINH_H__

#include <functional>

class TanhSinh {

  public:
	/**
     *      Extracted from wikipedia article 
     *
     *      \int_{-1}^1 f(x) dx \approx \sum_k={-\infty}^\infty \omega_k f(x_k)
     *
     *      with nodes at
     *
     *          x = \tanh\left( \tfrac{1}{2} \pi \sinh t   \right)
     *
     *      and weights
     *         
     *        w_k = \frac{\tfrac{1}{2} h \pi \cosh (k h)}
     *                   {\cosh^2 \left(\tfrac{1}{2} \pi \sinh (k h) \right)
     *
     *
     *
     *
     *
     *
	*/
  static cmplxd integrate(std::function<cmplxd (double)> func, double a, double b, int n=9)
  {

      const double h = 2./n;

      // calculate nodes and weights on the fly
      auto x =  [=] (int k) -> double { return tanh(0.5 * M_PI * sinh(k * h)) ; };
      auto w =  [=] (int k) -> double { return 0.5 * h  * M_PI * cosh(k * h) / pow2(cosh(0.5 * M_PI * sinh( k *h))) ;};

	  const double A = 0.5*(b-a);
	  const double B = 0.5*(b+a);

      cmplxd s = 0. ;
      for (int i=-n ; i<=n; i++) s += w(i) *  func(A*x(i) + B );

	  return A*s;
  }

};

#endif // __TANHSINH_H__
