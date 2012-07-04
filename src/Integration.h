/*
 * =====================================================================================
 *
 *       Filename: Integration.h
 *
 *    Description: Helper file for Gaussian-Integration. 
 *
 *                 However should be 
 *                 done as an interface with proper Gauss-Hermite/Tanh-Sinh or
 *                 Linear intergration.
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 *
 *
 *        Notes  : Improve implementation, looks quite ugly.
 *
 *
 * =====================================================================================
 */


#ifndef __INTEGRATE_H_
#define __INTEGRATE_H_


#include "Special/GaussLegendreWeights.h"

class Integration {
Array1d points;
Array1d weights;
  
  public:

    Integration(Range RLocal, int N, double L, bool useGauss) {
      weights.resize(RLocal);
      points.resize(RLocal);

      if(N==1) {
        weights = 1.;
        points = 0.;
        return;
      }
     if(useGauss == true) {
     
       check((N > GAUSS_WEIGHTS_MAX) ? -1 : 1, DMESG("Gaussian Integration Implemented until 48 Grid Poiints only. Run GaussIntegrate again")); 
  
      // maps points t, this is very stupid as we have ghost cells too 
      double * points_arr = gauss_weights[N].points;
      double *weights_arr = gauss_weights[N].weights;

      for(int i = (RLocal.first()-1); i <= (RLocal.last()-1); i++) {
        int pos = abs(i - (N-1)/2); 
        // fix for even points
        pos += ((!(N%2) && ((i - (N-1)/2) > 0)) ? -1 : 0);
           
        double p = (((i - (N-1)/2) <= 0) ? -1 : 1) * points_arr[pos];
        points(i+1) = L/2. + L/2. * p;
        weights(i+1) = L/2. * weights_arr[pos];
      }
    } else {
      // skip over 0 ? 
      //for(int i = RLocal.first(); i <= RLocal.last(); i++) points (i) = (i-1) * L/((double) N);
      for(int i = RLocal.first(); i <= RLocal.last(); i++) points (i) = (i-1) * L/((double) N);
       weights = L/((double) N);
    }

};
Array1d getPoints() {
    return points;

};
Array1d getWeights() {
    return weights;
};

Array3d sum(Array3d A, int n) {
    const double weight = weights(n);
    A = weight * A;
   return A; 

};

double d(int n) const {
    return weights(n);
};

};

#endif // __INTEGRATE_H_



