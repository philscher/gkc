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



/**
*   @brief solved the field equations using the Fourier method
*
*   Source terms Q are transformed to Fourier space in x-direction
*   as a complex-to-complex transformation. 
* 
*   \f[
*        Q(x,y_k,z) \rightarrow \hat{Q}(k_x,k_y,z)
*   \f]
*
*   @note that this  implicitly assumes periodicity in x-direction.
*
*   By gyro-ordering the contributions of the field in z-direction
*   can be neglected.
*
*   A huge advantage of Fourier method is that the gyro-averaging
*   operator \f[ <left \cdot right> \f] takes a particular simple
*   for of
*
*   \f[
*        \left< Q \right> = J_0(\sqrt{2 k_\perp^2 \rho_{t\sigma}\mu}) Q
*   \f]
*   where \f$ J_0 \f$ is the Bessel function of first order. The
*   double averaging which arises in the Poisson equation can
*   be analytically integrated.
*
*   \f[
*       I_0(b) e^{-b} = int_{\mu=0}^\infty J_0^2 e^{-\mu} quad,
*   \f]
*   with \f$ b = k_perp^2 rho_{t\sigma}^2 \f$ and $I_0$ is the
*   modified Bessel function of first kind and zeroth order.
*   Further we use the common definition \f$ \Gamma_0(b) = I_0(b) e^{-b}  \f$.
*   
*   Geometry factors enter the equations by a modification of the
*   perpendicular wavenumber (?) 
*   \f[
*        k_\perp^2 = g_xx k_x^2 + g_{xy} k_x k_y + g_yy k_y^2 
*   \f]
*
*   
*   Fourier transform is implemented by using FFTSolver class.
*   Recommended version is fftw3-mpi/fftw3. 
*   But also wrapper to fftw2 or other exists and can be used 
*   (although will require some fixed)
*
*
*/ 
class FieldsFFT : public Fields, public Fourier3D {

  protected:

  /** 
  *  @brief Calculation of the flux-surface average
  *
  * \f[
  *      \left< \phi \right>_{yz}(k_x) = \phi_k(k_x, k_y=0, k_z =0)
  *
  * \f]
  * solve the averaged Fields equation over y-z plane only for FFT
  * n this should be OK. For slab geometry this is k(kx,0,0) !
  * Result is given in Fourier space which is only k_x dependent. For solving
  * the fields equation, we need to subtract it from the adiabatic response, but
  * only for(ky=0, kz=0) modes !!!.
  **/
  virtual Array1z calcFluxSurfAvrg(Array4z rho_kykx);

  /**
  *   @brief Set if Nyquist frequency is removed from FFT spectra
  *   
  *   The Nyquist mode is unphysical and should be removed
  *   @todo explain why ??
  *
  **/
  bool screenNyquist;

  /**
  *
  *   @brief Solved the field equations
  *
  *   Implements filtering. 
  *   Calling of various solvers (coupled and uncoupled)
  *
  *
  **/
  virtual Array4z solveFieldEquations(Array4z fields, Timing timing);

  /**
  *
  *  @brief Solves the gyrokinetic Poisson equation
  * 
  *  The equation to solve is
  *
  *  \f[
  *      \lambda_D^2 \phi + \sum_\sigma frac{q_\sigma^2 n_\sigma}{T_\sigma}
  *       left( 1 - \Gamma_0(k_\perp^2\rho_{t\sigma}^2) \right) + \frac{q_0^2 n_0}{T_0} 
  *       \left(\phi - \left< right>_{yz} \right)
  *       = \rho 
  *  \f]
  *  where \f$ \lambda_D$ is the Debye length. The part latter part on the LHS
  *  is due to flux-surface averaging effect of an electron-species.
  *
  *  @warning We need to normalize the FFT transform here.
  *
  */
  Array3z virtual solvePoissonEquation(Array3z rho, Timing timing);
  
  /**
  *
  *  @brief Solves the gyrokinetic Ampere equation for the parallel vector potential \f$ A_{1\parallel} \f$
  *       
  *  The equation to solve is
  *
  *  \f[
  *      \left( \nabla_\perp^2 - 
  *      Y \hat{e} \beta \sum_\sigma sigma_sigma alpha_sigma^2 q_\sigma Gamma_0(b_\sigma)
  *      \right)
  *       = j_\parallel 
  *  \f]
  *  where \f$ Y\hat{e} \beta = \tfrac{1}{2}\hat{e} \beta$.
  *
  *  Reference : Eq. (2.40) of @cite Dannert_2006:PhDThesis
  *
  *  @warning We need to normalize the FFT transform here.
  *
  */
  Array3z virtual solveAmpereEquation (Array3z   j, Timing timing);

  /**
  *
  *
  *  @brief Solves the gyrokinetic Ampere equation for the parallel magnetic field \f$ B_{1\parallel} \f$
  *
  *  Reference : @cite Maerz_2008:PhDThesis
  *
  *  In case of including perturbation of the parallel magnetic field the field equations are coupled  
  *  and are modified to
  * 
  *  \f[
  *     C_1 = k_perp^2 lambda_D^2 + \sum_\sigma \frac{q_\sigma^2}{T_{0\sigma}} \left( 1 - \Gamma_0(b_\sigma) \right) \\
  *     C_2 = - \sum_\sigma frac{q_\sigma n_{0\sigma}}{B_0} Delta(b_\sigma)                                          \\
  *     C_3 = \frac{2}{\beta} - \sum_\sigma \frac{2 T_{0\sigma} n_{0\sigma}}{B_0^2} \Delta(b_\sigma)
  *  \f]
  *   
  *  The coupled Poisson-Ampere equation takes the form
  *
  *  \f[
  *      \left( \begin{array}{cc} C_1 & C_2 \\ C_2 & C_3 \right)
  *      \cdot
  *      \left( \begin{array}{c} \phi \\ B_\parallel \right)
  *      =
  *      \left( \begin{array}{c} \rho \\ j_\perp     \right)
  *  \f]
  * 
  *   It can be decoupled to (see reference for details)
  *
  *  \f[
  *      \left( \begin{array}{c} \phi \\ B_\parallel \right)
  *      =
  *      \frac{1}{C_1 C_3 - C_2^2}
  *      \left( \begin{array}{c} C_3 \rho - C_2 j_\perp \\ C_1 j_\perp - C_2 \phi \right)
  *  \f]
  *
  *
  *  Field equations for \f$ \phi \f$ and \f$ B_\parallel \f$ are now decoupled and can be solved
  *  subsequently.
  *    
  *  @warning We need to normalize the FFT transform here.
  *
  */
  Array3z virtual solveBParallelEquation(Array3z phi, Timing timing);

  /**
  *   @brief holds the flux surface average
  *
  *
  *
  *
  */ 
  Array1z  phi_yz;

  public:

  /** 
  *
  *  @brief performs the gyro-averaging of the field. 
  * 
  *  Performs gyro-average of fields. Depending on the settings, this function
  *  calls :
  *  
  *   - gyro-First (first order FLR effects)
  *   - gyroFull (full gyro-kinetic)
  *   - Nothing (drift-kinetic)
  *
  *   @todo link to function.
  *
  */
  Array4z gyroAverage(Array4z A, int m, int s, int nField, bool gyroField=false);
    
  /**
  * 
  *   @brief gyro-averaging for \f$ f_1(v_\parallel,\mu) \propto f_1(v_\parallel) exp(-\mu) \f$
  *          approximation
  *   
  *   Effectively approximate finite Larmor radius effects to first order. FLR damping 
  *   effects are overestimated, however, with the benefit that we can analytically integrate over
  *   the perpendicular velocity direction and thus reduce dimensionality by one.
  *   However, requires modification of the Vlasov equation.
  *
  *   Reference : @cite Dorland_Hammet_1993:GyroFluidwithKineticEffects
  *   
  *  \f{align}{
  *      \left< \phi       \right> &= exp(-k_\perp^2 \rho_\sigma ) \phi       (k_x, k_y, z) \quad, \\
  *      \left< j_parallel \right> &= exp(-k_\perp^2 \rho_\sigma ) j_\parallel(k_x, k_y, z) \quad, \\
  *  \f}
  *  
  *  Note : Not sure if this gyro-averaging is valid for \f$ j_\parallel \f$ and how to do for 
  *         \f$ B_{1\parallel} \f$.
  *
  *  @param A         the field or source terms
  *  @param m         index of magnetic moment
  *  @param s         index to species
  *  @param nField    index of field (phi, jp, jo) not used
  *  @param gyroField forward- or backward transformation
  *
  */
  Array4z gyroFirst(Array4z fields, int m, int s, int nField, bool gyroField=false);
  
  /**
  *   @brief performs gyro-averaging in Fourier space
  *  
  *  \f{align}{
  *      \left< \phi       \right> &= J_0(k_\perp^2 \rho_\sigma ) \phi       (k_x, k_y, z) \quad, \\
  *      \left< j_parallel \right> &= J_0(k_\perp^2 \rho_\sigma ) j_\parallel(k_x, k_y, z) \quad, \\
  *      \left< \phi       \right> &= I_1(k_\perp^2 \rho_\sigma ) B_\parallel(k_x, k_y, z) \quad. \\
  *  \f}
  *
  *  \f$ J_0 \f$ is the Bessel function of first kind and \f$ I_1 \f$ is the second order modified 
  *  Bessel function zeroth kind
  *
  *  @todo rename gyro-field to something better, e.g isForwardAverage
  *
  *  @param A         the field or source terms
  *  @param m         index of magnetic moment
  *  @param s         index to species
  *  @param nField    index of field (phi, jp, jo) not used
  *  @param gyroField forward- or backward transformation
  *
  */
  Array4z gyroFull (Array4z fields, int m, int s, int nField, bool gyroField=false);
  
  /**
  *   Constructor
  */ 
  FieldsFFT(Setup *setup, Grid *grid, Parallel *parallel, FileIO * fileIO, Geometry<HELIOS_GEOMETRY> *geo, FFTSolver *fftsolver);
  /**
  *   Constructor
  */ 
  ~FieldsFFT();
  
  /**
  *    @brief calculates the field energy
  *
  **/ 
  void calculateFieldEnergy(Array4z Q, double& phi, double& Ap, double& Bp);
  
protected:

  /**
  *   @brief Print out some runtime information
  */ 
  virtual void printOn(ostream &output) const;


  /**
  *  @brief helper function for 
  *  \f[
  *     r = \sum_\sigma \frac{q_\sigma^2 n_{0\sigma}}{T_{0\sigma}} (1 - \Gamma_0(b_\sigma)) 
  *  \f]
  *  
  *
  *  @note is it worth to precalculate this values ?
  *
  *  @param   \f$ k_\perp^2 \f$ (not normalized) perpendicular wavenumber
  *  @return  \f$ r \f$
  *
  **/
  inline  double sum_qqnT_1mG0(const double k2_p) ;

  /**
  *  @brief helper function for 
  *  
  *  \f[
  *     r = \sum_\sigma \sigma_\sigma alpha_\sigma^2 q_\sigma \Gamma_0(b_\sigma) 
  *  \f]
  *
  *  @note is it worth to pre-calculate this values ?
  *
  *  @param   \f$ k_\perp^2 \f$ (not normalized) perpendicular wavenumber
  *  @return  \f$ r \f$
  *
  **/
  inline  double sum_sa2qG0(const double kp_2) ;

  /**
  *
  *  @brief helper function for 
  *  
  *  \f[
  *     r = frac{1}{B_0} \sum_\sigma q_\sigma n_{0\sigma} \Delta(b_\sigma) 
  *  \f]
  *  where we used \f$ \Delta = I_0 - I_1 \f$
  *  
  *  @note is it worth to pre-calculate this values ?
  *
  *  @param   \f$ k_\perp^2 \f$ (not normalized) perpendicular wavenumber
  *  @return  \f$ r \f$
  **/
  inline  double sum_qnB_Delta(const double k2_p) ;

  /**
  *
  *  @brief helper function for 
  *  
  *  \f[
  *     r = frac{2}{B_0^2} \sum_\sigma T_{0\sigma} n_{0\sigma} \Delta(b_\sigma) 
  *  \f]
  *  where we used \f$ \Delta = I_0 - I_1 \f$
  *  
  *  @note is it worth to pre-calculate this values ?
  *
  *  @param   \f$ k_\perp^2 \f$ (not normalized) perpendicular wavenumber
  *  @return  \f$ r \f$
  **/
  inline  double sum_2TnBB_Delta(const double k2_p) ;


};


#endif // __FIELDS_FFT_H
