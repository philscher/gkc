/*
 * =====================================================================================
 *
 *       Filename: FieldsFFT.h
 *
 *    Description: Fields Solver in Fourier space
 *
 *         Author: Paul P. Hilscher (2009-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __FIELDS_FFT_H
#define __FIELDS_FFT_H

#include "Global.h"

#include "Setup.h"
#include "Parallel.h"
#include "Grid.h"
#include "Fourier.h"
#include "Geometry.h"
#include "Fields.h"
//! Fields - class for solving the Fields's equation
/*!
*  This class us the non-MPI version. But is has been parellized by using
*  OpenMP. Unfortunately as most of the function are memory operation parallelization
*  efficiency is rather low.
*  For the FFT solver, this class rely on FFTW (www.fftw.org). This library is very efficient
*  and if compiled for parallism (pthread, OpenMP or/and CUDA) it can make us of
*  the avaible CPUs/Cores very well.
*
*  We us FFTW to transform into Fourier space using the real2complex function. This function
*  takes real input and produces complex output. Becausefor real transform half of the
*  values are complex conjugant to each other and reduntant.
*
*  After transforming a 3 dimensions n_x * n_y * n_z array of real variables the output
*  is n_x * n_y * (n_z/2 + 1). We should remember that for r2c/c2r transforms the
*  input array is destroyed !
*
*
*  This class uses p3dfft to perform the parallized Fourier transformation.
*  This is done using a 2D decompostion in X and Y. The Vlasov solver
*  is adapted to p3dfft to follow their order
*
*
*  p3dfft proviedes basically 4 functions:
*       p3dfft_setup : Initialiization, we need to specify the decomposition in
*                      (nx, ny, nz) direction. Because decompoisiton is in
*                      ny and nz (!), we need to initialize as (nz, nx, ny) and
*                      switch the variables accordingly
*       get_dims     : We get the decomposition of the grid in real (conf = 1)
*                      and Fourier space (conf1). We need to keep in mind that
*                      the data in Fourier space is in trnasposed form !
*
*  Note that all the array properties for input/output are hidden inside the r3In 
*  k3Out classes, in face they are
*
*  FFTW     :  r3In is normal; k3Out is normal
*  FFTW-MPI :  r3In is normal; k3Out(y,x,z)
*  P3DFFT   :  r3In(y,x,z),    k3Out(y,x,z) (don't set DSTRIDE)
*
*  all this is hidden inside the GeneralArray storage class, so we can access
*  r3In and k3Out as (x,y,z) !!
*
*
*/ 
class FieldsFFT : public Fields, public Fourier3D {
protected:
/** 
 *  Calculation of the ITG potential.
 * \f[
 *      \left< \phi \right>_{yz}(k_x) = \phi_k(k_x, k_y=0, k_z =0)
 *
 * \f]
 * solve the averaged Fields equation over y-z plane only for FFT
 * n this should be OK. For slab geometry this is k(kx,0,0) !
 * Result is given in Fourier space which is only k_x dependent. For solving
 * the fields equation, we need to substract it from the adiabatic repsonse, but
 * only for(ky=0, kz=0) modes !!!.
 **/
virtual Array1z calcFluxSurfAvrg(Array4z rho_kykx);


virtual Array4z solveFieldEquations(Array4z fields, Timing timing);
/*!
 *
 *  Solver the gyrokinetic Fields's equation : 
 * 
 *  \f[
 *      -k^2 \hat{\phi} + \hat{\phi} \sum_\sigma \frac{1}{\lambda_{D\sigma}} 
 *      \left[ \Gamma_0(k_\perp^2 \rho_{t\sigma})  - 1 \right] = \pi \hat{rho}(k_x,k_y,k_z)
 *  \f]
 *
 */
Array3z virtual solvePoissonEquation(Array3z rho, Timing timing);
Array3z virtual solveAmpereEquation (Array3z   j, Timing timing);
Array3z virtual solveBParallelEquation(Array3z phi, Timing timing);
Array1z  phi_yz;
public:


/*! Performs the gyro-average for variable stored in FFT->r3In
 *  
 *  We perform the gyro-average in Fourier-space according to
 *  
 *  \f[ \bar{\phi} = \sum_k J_0(k_\perp^2 \rho_\sigma ) \phi(k_x,k_y,k_z) \f]
 *
 *  or scaled
 *
 *
 *  \f[ \bar{n}(x) = \frac{B_0}{m} \int J_0(\lambda) F dv_\parallel d\mu d\theta \f]
 *
 */
  Array3z gyroAverage(Array3z fields, int m, int s, int nField, bool gyroField=false);
    /**
     *   Lets use Dorland Filtering, which means multiplications with exp(-k_p square)
     *
     *   \f[ \exp \left( - b_\sigma \right) \f]
     *
     *   Gyro-averaging apprximation for \f[ g_{1 sigma} = f(z,\v_\parallel) \cdot f(z, \mu)_M = f(z,v_\parallel) exp(-\frac{\hat{\mu} B}{T}) \f]
     *
     * */
  Array3z gyroFirst(Array3z fields, int m, int s, int nField, bool gyroField=false);
  Array3z gyroFull (Array3z fields, int m, int s, int nField, bool gyroField=false);
 
  
  
  
/**
 * We can make symple approximation 
 *
 *
 *   f = f_M,v_perp f_v_parallel 
 *   and use exp(-k_perp^2) instead as gyro average,
 *   does not incude wave coupling or other kinetic effects
 *   
 *   Ref. : W.Dorland & G.W. Hammet
 *          Gyrofluid turbulence models with kinetic effects
 *          Phys. Fluids B 5 (3.), March 1993
 *
 *   
 *   Averages qs
 *   \f[  
 *         f = f_M(v_\perp) f(v_\parallel)
 *   	   <phi_k> = \phi_k \exp \left[ - \frac{1}{2} k_perp^2 \right]
 *   \f]        
 *
 **/

  //* Constructor
  FieldsFFT(Setup *setup, Grid *grid, Parallel *parallel, FileIO * fileIO, Geometry<HELIOS_GEOMETRY> *geo, FFTSolver *fftsolver);
  //* Destructor
  ~FieldsFFT();
  
protected:

        virtual void printOn(ostream &output) const {
        Fields::printOn(output);
         output   << "Fields     |  FFT "  << std::endl;
        }
};


#endif // __FIELDS_FFT_H
