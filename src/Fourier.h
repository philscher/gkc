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
#include "SpecialMath.h"


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
       
        // doGyro = true
        if(plasma->species(s).doGyro == true) g0 += qqnT * SpecialMath::_1mGamma0( rho_t2 * k2_p);
        else g0 += qqnT * (rho_t2 * k2_p);
        
        //g0 += qqnT * _1mGamma0( rho_t2 * k2_p, plasma->species(s).doGyro);
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

        g0 += norm * SpecialMath::_1mGamma0( rho_t2 * k2_p);
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

        g0 += norm * SpecialMath::_1mGamma0( rho_t2 * k2_p);
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
        if(plasma->species(s).doGyro == true) g0  += sa2q * (1.e0 - SpecialMath::_1mGamma0(b));
        else g0 += sa2q;
    }
    return g0;
};

inline  double sum_qnB_Delta(const double k2_p) {
    double g0 = 0.;
    for(int s = NsGlD; s <= NsGuD; s++) {
        const double qn    = plasma->species(s).q * plasma->species(s).n0;
        const double rho_t2 = plasma->species(s).T0 * plasma->species(s).m / (pow2(plasma->species(s).q) * plasma->B0);

        g0 += qn * SpecialMath::Delta(rho_t2 * k2_p);
    }
    return g0/plasma->B0;
};

inline  double sum_2TnBB_Delta(const double k2_p) {
    double g0 = 0.;
    for(int s = NsGlD; s <= NsGuD; s++) {
        const double Tn   = plasma->species(s).T0 * plasma->species(s).n0;
        const double rho_t2 = plasma->species(s).T0 * plasma->species(s).m / (pow2(plasma->species(s).q) * plasma->B0);

        g0 += Tn * SpecialMath::Delta(rho_t2 * k2_p);
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


};


#endif // __Fourier_H
