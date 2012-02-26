/*
 * =====================================================================================
 *
 *       Filename: Fourier.h
 *
 *    Description: Helper functions for Fourier space
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __FOURIER_H
#define __FOURIER_H

#include "Global.h"
#include "Setup.h"
#include "Grid.h"
#include "FFTSolver.h"
#include "Geometry.h"
#include "GeometryShear.h"
#include "Geometry2D.h"

#include "Plasma.h"



static inline double BesselI0( const double x )
/* ------------------------------------------------------------
 *  PURPOSE: Evaluate modified Bessel function In(x) and n=0.  
 * ------------------------------------------------------------*/
{
     double ax,ans;


     if ((ax=fabs(x)) < 3.75) {
         const double y=pow2(x/3.75);
         ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
     } else {
         const double  y=3.75/ax;
         ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
          +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
                                                                        +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
                                                                                       +y*0.392377e-2))))))));
                                             }
              return ans;
}

static inline double BesselI1( double x ) {
        double I1;
                                                           
           if (abs(x) <  3.75)  { 
                    const double  p[8] = { 0.0, 0.5, 0.87890594, 0.51498869, 0.15084934, 0.2658733e-1, 0.301532e-2, 0.32411e-3 };
                    const double y=pow2((x/3.75));
                    I1=x*(p[1]+y*(p[2]+y*(p[3]+y*(p[4]+y*(p[5]+y*(p[6]+y*p[7]))))));
           } else {
                 const double  q[10] = { 0.0, 0.39894228e0, -0.3988024e-1, -0.362018e-2, 0.163801e-2, -0.1031555e-1, 
                                         0.2282967e-1, -0.2895312e-1, 0.1787654e-1, -0.420059e-2 } ;
                const double ax=abs(x);
                const double y=3.75/ax;
                I1=(exp(ax)/sqrt(ax))*(q[1]+y*(q[2]+y*(q[3]+y*(q[4]+y*(q[5]+y*(q[6]+y*(q[7]+y*(q[8]+y*q[9]))))))));
                      
           }
            return (x  > 0.) ? I1 : -I1;
}
  // This is valid for half complex order !



//static inline  double _1mGamma0(const double b, bool gyro=true) {
//    if(gyro) return  1. - BesselI0(b)*exp(-b);
//    else     return  b;
// };

static inline  double _1mGamma0(const double b, bool gyro=true) {
    if(gyro) return  b / ( 1.e0 + b);
    else     return  b;
 };

//__declspec(vector) static inline  double _1mGamma0(const double b) {
static inline  double _21mGamma0(const double b) {
    return  b / ( 1.e0 + b);
   // return  1. - BesselI0(b)*exp(-b);
 };

static inline  double Gamma0(const double b, bool gyro=true) {
    if(gyro) return  BesselI0(b)*exp(-b);
    else     return  b;
 };

static inline  double Gamma1(const double b, bool gyro=true) {
    if(gyro) return  BesselI1(b)*exp(-b);
    else     return  b;
 };

// there should be a simpler epxression for this @!!
static inline  double G0mG1(const double b, bool gyro=true) {
    if(gyro) return  Gamma0(b) - Gamma1(b);
    else     return  b;
 };

static inline  double I1sp(const double b) {
	const double i = 1.;
	return  -2. * i *  BesselI1(i * b)/b ;
};



/** Helper function to sum \f[ \sum_\sigma \frac{1}{\lambda_{D\sigma}} \left[ \Gamma_0( k_\perp^2 rho_{t\sigma}^2) - 1 \right] \f] 
 *   note that
 *  \left<  \lambda^2 = -\rho  \nabla_\perp^2 
 *                    = - \rho_t ( k_x^2 + k_y^2 )
 *
 **/
inline  double sum_qqnT_1mG0(const double k2_p) {
    double g0 = 0.;
    for(int s = NsGlD; s <= NsGuD; s++) {
        const double qqnT   = plasma->species(s).n0 * pow2(plasma->species(s).q)/plasma->species(s).T0;
        const double rho_t2 = plasma->species(s).T0  * plasma->species(s).m / pow2(plasma->species(s).q * plasma->B0);
        
        g0 += qqnT * _1mGamma0( rho_t2 * k2_p);
    }
    return g0;
};

/** Helper function to sum \f[ \sum_\sigma \frac{1}{\lambda_{D\sigma}} \left[ \Gamma_0( k_\perp^2 rho_{t\sigma}^2) - 1 \right] \f] 
 *   note that
 *  \left<  \lambda^2 = -\rho  \nabla_\perp^2 
 *                    = - \rho_t ( k_x^2 + k_y^2 )
 *
 **/
inline  double sum_phiG0mG1(const double k2_p) {
    double g0 = 0.;
    for(int s = NsGlD; s <= NsGuD; s++) {
        const double norm   = plasma->species(s).n0 * pow2(plasma->species(s).q)/plasma->species(s).T0;
        const double rho_t2 = plasma->species(s).T0 * plasma->species(s).m / (pow2(plasma->species(s).q) * plasma->B0);

        g0 += norm * _1mGamma0( rho_t2 * k2_p);
    }
    return g0;
};

/** Helper function to sum \f[ \sum_\sigma \frac{1}{\lambda_{D\sigma}} \left[ \Gamma_0( k_\perp^2 rho_{t\sigma}^2) - 1 \right] \f] 
 *   note that
 *  \left<  \lambda^2 = -\rho  \nabla_\perp^2 
 *                    = - \rho_t ( k_x^2 + k_y^2 )
 *
 **/
inline  double sum_BpaG0mG1(const double k2_p) {
    double g0 = 0.;
    for(int s = NsGlD; s <= NsGuD; s++) {
        const double norm   = plasma->species(s).n0 * pow2(plasma->species(s).q)/plasma->species(s).T0;
        const double rho_t2 = plasma->species(s).T0 * plasma->species(s).m / pow2(plasma->species(s).q * plasma->B0);

        g0 += norm * _1mGamma0( rho_t2 * k2_p);
    }
    return g0;
};


/** Helper function to sum \f[ \sum_\sigma \frac{1}{\lambda_{D\sigma}} \left[ \Gamma_0( k_\perp^2 rho_{t\sigma}^2) - 1 \right] \f] 
 *   note that
 *  \left<  \lambda^2 = -\rho  \nabla_\perp^2 
 *                    = - \rho_t ( k_x^2 + k_y^2 )
 *
 **/
inline  double sum_sa2qG0(const double kp_2) {
    double g0 = 0.;
    for(int s = NsGlD; s <= NsGuD; s++) {
        const double sa2q   = plasma->species(s).sigma * pow2(plasma->species(s).alpha) * plasma->species(s).q;
        const double rho_t2 = plasma->species(s).T0 * plasma->species(s).m / (pow2(plasma->species(s).q) * plasma->B0);

        const double b = rho_t2 * kp_2;
        //g0 += sa2q * (1.e0 - _1mGamma0( rho_t2 * kp_2, plasma->species(s).doGyro));
        
        
        //if(plasma->species(s).doGyro == true) g0 += sa2q * I0(b)*exp(-b);
        if(plasma->species(s).doGyro == true) g0  += sa2q * (1.e0 - _1mGamma0(b));
        else g0 += sa2q;
        //g0 += sa2q * (1.e0 - _1mGamma0( rho_t2 * kp_2, plasma->species(s).doGyro));
    }
    return g0;
};

/**
 *   \f[     \Delta = \left( I_0(x) - I_1(x) \right) \exp(-x)  \f]
 *
 * */
inline  double Delta(const double x) {
    return (BesselI0(x)-BesselI1(x))*exp(-x);
}

inline  double sum_qnB_Delta(const double k2_p) {
    double g0 = 0.;
    for(int s = NsGlD; s <= NsGuD; s++) {
        const double qn    = plasma->species(s).q * plasma->species(s).n0;
        const double rho_t2 = plasma->species(s).T0 * plasma->species(s).m / (pow2(plasma->species(s).q) * plasma->B0);

        g0 += qn * Delta(rho_t2 * k2_p);
    }
    return g0/plasma->B0;
};

inline  double sum_2TnBB_Delta(const double k2_p) {
    double g0 = 0.;
    for(int s = NsGlD; s <= NsGuD; s++) {
        const double Tn   = plasma->species(s).T0 * plasma->species(s).n0;
        const double rho_t2 = plasma->species(s).T0 * plasma->species(s).m / (pow2(plasma->species(s).q) * plasma->B0);

        g0 += Tn * Delta(rho_t2 * k2_p);
    }
    return 2. * g0 /pow2(plasma->B0) ;
};

class Fourier3D {



    std::vector<int> suppressModeX, suppressModeY, suppressModeZ;
    std::vector<int> convolveModeX, convolveModeY, convolveModeZ;


double epsilon_0, sigma;
  protected:
 FFTSolver *fft;



  public:
  
  Fourier3D(Setup *setup, Grid *grid, FFTSolver *fftsolver);
  virtual  ~Fourier3D();

  std::string getSolverInfo();
  
   //! If parallel mode is on it boundaries to other CPUs using MPI
   /*!
    * Because we need a 2-dimensional decomposition strategy, we need to send exchange
    * the values to between cpus to for different nodes.
    * See also the Vlasov solver for more details 
    *                
    *
    */ 

  
    //! k get the wavenumber 
    /*! Note that we have a symmetry about the Nyquist fequency which
     *  is the highest frequency one can use (half of the sampling frequency)
     *  and FFT includes positive AND negative frequencies !
     */
    /*!  We can suppress various modes, this is set in the setupuration fields.
     *   and sets the Fourier mode to zero.
     */
    int suppressModes(Array4z k2Out, const int field=1);
    int suppress3DMode(Array4d A);


};


#endif // __Fourier_H
