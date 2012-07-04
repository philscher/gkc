/*
 * =====================================================================================
 *
 *       Filename: Integrate.h
 *
 *    Description: Function to support high order integration of a function 
 *                 or by direct summation.
 *
 *         Author: Paul P. Hilscher (2012-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __GAUSS_LEGENDRE_H__
#define __GAUSS_LEGENDRE_H__

#include "GaussLegendreWeights.h"
#include "GaussLaguerreWeights.h"
#include "GaussRadauWeights.h"

#include <functional>

/**
 *
 *   Integration class for various types
 *
 *
 *
 *   Gauss-Radau : Note that -1 is our fixed node.
 *
 *
 *   Gauss Quadrature formulas only defined between [-1,1], thus for integrating
 *   [ a, b]  we do rescaling
 *    
 *    int_a^b f(x)
 *     
 *    w_i * (b-a)/2 * int_1^{-1} f((b-a)/2 * x_i + (a+b)/2)
 *
 *
 *
 **/
class Integrate {

  public:
 //template<class T, class R> static T integrate(int n, T (*f)(R), R a, R b)
 //template<class T> static T integrate(int n, T (*func)(double), double a, double b)
 //template<class T> static T integrate(int n, std::function func, double a, double b)
private:

  static cmplxd integrate(std::function<cmplxd (double)> func, double a, double b, int n, const double x[], const double w[])
  {
	const double A = 0.5*(b-a);
	const double B = 0.5*(b+a);

    cmplxd s = 0. ;
    for (int i = 0 ; i <= n; i++) s += w[i] *  func(A*x[i] + B );

	return A*s;
  }


  double *weights;
  double *nodes;


  const int order;
public:

  Integrate(std::string type, int _order = 11, double a=0., double b=1.) : order(_order) {
       
         auto setupNodesAndWeights = [=] (double *x, double *w, bool rescale) { 
            
	        const double A = (rescale) ? 0.5*(b-a) : 1.;
	        const double B = (rescale) ? 0.5*(b+a) : 0.;

                nodes   = new double[order];
                weights = new double[order];
                
                // Rescale
                for(int n = 0; n < order; n++)  {
                    nodes[n]   = A * x[n] + B;
                    weights[n] = A * w[n];
                }

         };

        if     (type == "Gauss-Legendre") { setupNodesAndWeights(gauss_weights[order].points   , gauss_weights[order].weights,   true); }
        else if(type == "Gauss-Laguerre") { setupNodesAndWeights(Laguerre_weights[order].points, Laguerre_weights[order].weights, false); }
        else if(type == "Gauss-Radau"   ) { setupNodesAndWeights(Radau_weights[order].points   , Radau_weights[order].weights,   true); }
        else   check(-1, DMESG("No such quadrature rule"));
        

  };

 ~Integrate()
 {
    delete[] nodes;
    delete[] weights;
    
 };

 // return nan else
 double x(const int n) const { return (n < order) ? nodes[n] : 0. ; };
 // return nan else
 double w(const int n) const { return (n < order) ? weights[n] : 0.; };

 /*   
  static cmplxd i/ntegrate(std::function<cmplxd (double)> func, double a, double b, int n=9)
    {

    // Maximum Gauss Weigths/
    if (n > GAUSS_WEIGHTS_MAX) {
      check(-1, DMESG("Maximum order of intergration limited to 47"));
      return 0.;
    }

    // Table has a shift by 1, thus n=0 corresponds to w[1] = { 1.}
    const double *x = gauss_weights[n].points;
    const double *w = gauss_weights[n].weights;

	const double A = 0.5*(b-a);
	const double B = 0.5*(b+a);

    cmplxd s = 0. ;
    for (int i = 0 ; i <= n; i++) s += w[i] *  func(A*x[i] + B );

	return A*s;
  }
  */
  static cmplxd GaussLegendre(std::function<cmplxd (double)> func, double a, double b, int n=9) {

        return integrate(func, a, b, n, gauss_weights[n].points, gauss_weights[n].weights);

  }
  
  static cmplxd GaussRadau(std::function<cmplxd (double)> func, double a, double b, int n=9) {

        return integrate(func, a, b, n, Radau_weights[n].points, Radau_weights[n].weights);
  
  }
	
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
  static cmplxd TanhSinh(std::function<cmplxd (double)> func, double a, double b, int n=9)
    {

      const double h = 2./n;

      // calculate nodes and weights on the fly
      double x[n], w[n];

      for(int k = 0; k < n; k++) {
         x[k] =  tanh(0.5 * M_PI * sinh(k * h)) ; 
         w[k] =  0.5 * h  * M_PI * cosh(k * h) / pow2(cosh(0.5 * M_PI * sinh( k *h))) ;
      }

      return integrate(func, a, b, n, Radau_weights[n].points, Radau_weights[n].weights);

  }



};

#endif // __GAUSS_LEGENDRE_H__
