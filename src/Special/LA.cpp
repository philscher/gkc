/*
 * =====================================================================================
 *
 *       Filename: LA.cpp
 *
 *    Description: Linear Algebra Interface
 *
 *         Author: Paul P. Hilscher (2009-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#include "Global.h"
#include "Special/LA.h"

//int solveCycCCSymTriDiagonal(const double *b, double *x, const double *cs_alpha, const int *N) {
int solvecycccsymtridiagonal_(const double *b, double *x, const double *cs_alpha, const int *N) {
    return solveCycCCSymTriDiagonal(b, x, *cs_alpha, *N); 
};




// Solve Tridiagonal Cyclic system with constant coefficient sub and sup diagonal and diagonal = 1
int solveCycCCSymTriDiagonal(const double *b, double *x, const double cs_alpha, const int N)
{
      

      double delta[N], gamma[N], alpha[N], c[N], z[N];

      double sum = 0.0;


      if (N < 3) check(-1, DMESG("Use with 3 or more elements")); 

      alpha[0] = 1.;
      gamma[0] = cs_alpha;
      delta[0] = cs_alpha;
      sum += alpha[0] * pow2(delta[0]);

      for (int i = 1; i < N - 2; i++)
        {
          alpha[i] = 1. - cs_alpha * gamma[i - 1];
          gamma[i] = cs_alpha / alpha[i];
          delta[i] = -delta[i - 1] * cs_alpha / alpha[i];
          sum += alpha[i] * pow2(delta[i]);
        }


      alpha[N - 2] = 1. - cs_alpha * gamma[N - 3];
      gamma[N - 2] = (cs_alpha - cs_alpha * delta[N - 3]) / alpha[N - 2];
      alpha[N - 1] = 1. - sum - alpha[(N - 2)] * pow2(gamma[N - 2]);
      

      /* update */
      z[0] = b[0];
      for (int i = 1; i < N - 1; i++) z[i] = b[i] - z[i - 1] * gamma[i - 1];
      sum = 0.0;
      for (int i = 0; i < N - 2; i++)   sum += delta[i] * z[i];
      z[N - 1] = b[N - 1] - sum - gamma[N - 2] * z[N - 2];
      for (int i = 0; i < N; i++) c[i] = z[i] / alpha[i];

      /* backsubstitution */
      x[N - 1] = c[N - 1];
      x[N - 2] = c[N - 2] - gamma[N - 2] * x[N - 1];
          
      for (int i = N - 3, j = 0; j <= N - 3; j++, i--) x[i] = c[i] - gamma[i] * x[i + 1] - delta[i] * x[N - 1];

  return 1;
}




