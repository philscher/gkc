/*
 * =====================================================================================
 *
 *       Filename:  gkc-iCode.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *
 *         Author:  Paul P. Hilscher (2012) 
 *
 * =====================================================================================
 */


#include "iCode-func.h"

#include "SpecialMath.h"

#include <fenv.h>
#include<iostream>              

template<class T> T pow2(const T z) { return z*z; };
template<class T> T pow4(const T z) { return pow2(pow2(z)); };
template<class T> T pow6(const T z) { return pow4(z) * pow2(z); };


#include<omp.h>


//typedef double float_t ;
//typedef std::complex<double> complex_t; 

inline  std::complex<double> PlasmaDispersion(std::complex<double> z) {
   
        static const double  SQRT_PI = 1.7724538509055160272981674833411451;
        static const std::complex<double> ii = std::complex<double>(0.,1.);
        std::complex<double>  Z;
       

        // sometimes imaginary part is too high, due to test functions, limit
        if(imag(z) > 25.) z = std::complex<double>(real(z), 25.);

        if(real(z*z) > 100) {
            //z = std::min(500., (real(z)*real(z) - imag(z)*imag(z)))  +  ii * 2. * real(z) * imag(z); 
            //Z = - 1./z * ( 1. + 1./(2.*pow2(z)) + 3./(4.*pow4(z)) + 15./(8.*pow6(z)));		
            const double x = real(z); 
            const double y = imag(z);
            /* 
            double sigma;
            if      (y           >  1./std::abs(x))  sigma = 0.;
            else if (std::abs(y) <  1./std::abs(x))  sigma = 1.;
            if    (y           < -1./std::abs(x))  sigma = 2.;
             * */
            
            const double sigma = (y >  0.) ? 0 : (( y == 0) ? 1. : 2.);
            Z = ii *   SQRT_PI * sigma * exp(-z*z) - 1./z * ( 1. + 1./(2.*z*z) + 3./(4.*pow4(z)) + 15./(8.*pow6(z)));		
  
        }
        else Z = ii * SQRT_PI * std::exp(- z*z) * MathFunctions::erfc(- ii * z);
    return Z;
}; 




std::complex<double> getL(std::complex<double> w, double ky, double kx, double kx_, double *X,  double Ls,  double Ln, const int Nx,
            const double q, const double mass, const double T, const double eta) 

{

  static const std::complex<double> ii = std::complex<double>(0.,1.);
  
  std::complex<double> R=std::complex<double>(0.,0.);

  const double t = 0.;

  // We can ignore  terms with kx - kx_ >> 1, as integration over x is 0 due to oscillating term std::exp( ii * (kx - kx_) * X[x]);

  // Integrate over x, using fixed kx, kx_
  for(int x = 0; x < Nx; x++) {

         if(X[x] == 0.) continue;
         // ignore modes to far away as physically no connection but numerical problems
         // if(std::abs(kx-kx_) > 5.) continue;

         const double kp2    = ky*ky + kx * kx ;
         const double kp2_   = ky*ky + kx_* kx_;
     
         const double Omega  =  - q  / mass;

         const double v_th   = sqrt(2. * T / mass);
         const double rho    = v_th / Omega ;

         const double w_star = ky * T / ( Omega * mass * Ln );
         const double w_D    = 0. ;//#- w_star * Ln / LB;
         const double k_p    = (X[x]/Ls) * ky;
     
         const double b  = kp2  * rho*rho / 2.;
         const double b_ = kp2_ * rho*rho / 2.;
     
         const std::complex<double>  zeta   = (w - w_D * t) / (v_th * std::abs(k_p));
         const std::complex<double>  Z      = PlasmaDispersion(zeta) ;
     
         // Note : b_a >= b_g for every b,b_
         const double b_a = (b + b_) / 2.; // arithmetic mean
         const double b_g = sqrt(b * b_);  // geometric  mean
              
         // need asymptoctic approximation otherwise it fails, (points corresponds where
         // errors from double precision are large than due to expansion.
         const double G0 = (b_g < 10.) ? MathFunctions::i0(b_g) * exp(-b_a) 
                                      : exp(b_g - b_a)/sqrt(2.*M_PI*b_g) * ( 1. + 1./(8.*b_g) +  9./(128.*b_g*b_g));
         const double G1 = (b_g < 10.) ? MathFunctions::i1(b_g) * exp(-b_a)
                                      : exp(b_g - b_a)/sqrt(2.*M_PI*b_g) * ( 1. - 3./(8.*b_g) - 15./(128.*b_g*b_g));
         const std::complex<double> z = w_star/w; 

         const std::complex<double>    Lq  =  (1.   - z) * zeta * Z * G0  
                                - eta * z * (pow2(zeta) + (pow2(zeta) -3./2.) * zeta * Z) * G0 
                                - eta * z * zeta * Z * ((1. - b_a) * G0 + b_g * G1);


         // Note : int_0^2\pi exp(1.j * x) = 0.
         R += Lq *  std::exp( ii * (kx_ - kx) * X[x]);
  } 
  
  return R;

};


void setupMatrix(std::complex<double> w, double ky, double *X, double *K,  double Ls,  double Ln, int Nx, std::complex<double> *M, double dh,
                double q, double mass, double T, double eta, double n) 
{
       
      #pragma omp parallel for
      for(int x_k = 0; x_k < Nx; x_k++) {  for(int w_k = 0; w_k < Nx; w_k++) {
         
        const int idx = x_k*Nx + w_k;              
        
          M[idx]  +=  - n * (q*q)  / T * ( (w_k == x_k ? 1. : 0. ) + 1./(2.*M_PI) *  getL(w, ky, K[x_k], K[w_k], X, Ls, Ln, Nx, q, mass, T, eta) * dh );

       } }

      return;
}
