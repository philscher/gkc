 // Special functions -*- C++ -*-

// small modification by 2011 (c) Paul P. Hilscher

 // Copyright (C) 2006, 2007, 2008 , GPL(v3+)
 // Free Software Foundation, Inc.
 //


#ifndef SPECIAL_H
#define SPECIAL_H


class Special {

  public: 
 //
 // ISO C++ 14882 TR1: 5.2  Special functions
 // Importaned from gcc   @file tr1/poly_hermite.tcc
 // Written by Edward Smith-Rowland based on:
 //   (1) Handbook of Mathematical Functions,
 //       Ed. Milton Abramowitz and Irene A. Stegun,
 //       Dover Publications, Section 22 pp. 773-802
 
     /**
      *   @brief This routine returns the Hermite polynomial
      *          of order n: \f$ H_n(x) \f$ by recursion on n.
      * 
      *   The Hermite polynomial is defined by:
      *   @f[
      *     H_n(x) = (-1)^n e^{x^2} \frac{d^n}{dx^n} e^{-x^2}
      *   @f]
      *
      *   @param n The order of the Hermite polynomial.
      *   @param x The argument of the Hermite polynomial.
      *   @return The value of the Hermite polynomial of order n
      *           and argument x.
      *   @brief This routine returns the Hermite polynomial
      *          of order n: \f$ H_n(x) \f$.
      * 
      *   The Hermite polynomial is defined by:
      *   @f[
      *     H_n(x) = (-1)^n e^{x^2} \frac{d^n}{dx^n} e^{-x^2}
      *   @f]
      *
      *   @param n The order of the Hermite polynomial.
      *   @param x The argument of the Hermite polynomial.
      *   @return The value of the Hermite polynomial of order n
      *           and argument x.
      */
     template<typename T> inline static T HermitePolynomial(const unsigned int n, const T x)
     {
       //  Compute H_0.
       T HP_0 = 1;
       if (n == 0) return HP_0;
 
       //  Compute H_1.
       T HP_1 = 2 * x;
       if (n == 1) return HP_1;
 
       //  Compute H_n.
       T HP_n, HP_nm1, HP_nm2;
       unsigned int i = 2;
       for  (HP_nm2 = HP_0, HP_nm1 = HP_1; i <= n; ++i)
         {
           HP_n = 2 * (x * HP_nm1 + (i - 1) * HP_nm2);
           HP_nm2 = HP_nm1;
           HP_nm1 = HP_n;
         }
 
       return HP_n;
     };

     static int factorial (int n) { return ((n ==1) || (n==0)) ? 1 : n * factorial(n-1); };
     
     template<typename T> static inline T HermiteFunction(const unsigned int n, const T x)
     {
       //normalization term
       const double norm  = 1./sqrt((sqrt(M_PI) * pow(2.,(int)n) * factorial(n)) );
       return norm * HermitePolynomial(n, x) * exp(-(x*x)/2.);
     } 
 
};

#endif // __SPECIAL_H
