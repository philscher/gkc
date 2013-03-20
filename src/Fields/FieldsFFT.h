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
#include "Grid.h"
#include "Parallel/Parallel.h"
#include "FFTSolver/FFTSolver.h"
#include "Geometry/Geometry.h"
#include "Fields.h"



/**
*   @brief solved the field equations using the Fourier method
*
*   Source terms Q are transformed to Fourier space in x-direction
*   as a complex-to-complex transformation. 
* 
*   \f[
*        Q(x,k_y,z) \rightarrow \hat{Q}(k_x,k_y,z)
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
*       I_0(b) e^{-b} = \int_{\mu=0}^\infty J_0^2 e^{-\mu} \quad,
*   \f]
*   with \f$ b = k_\perp^2 \rho_{t\sigma}^2 \f$ and \f$ I_0 \f$ is the
*   modified Bessel function of first kind and zeroth order.
*   Further we use the common definition \f$ \Gamma_0(b) = I_0(b) e^{-b}  \f$.
*   
*   Geometry factors enter the equations by a modification of the
*   perpendicular wavenumber (?) 
*   \f[
*        k_\perp^2 = g_{xx} k_x^2 + g_{xy} k_x k_y + g_{yy} k_y^2 
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
class FieldsFFT : public Fields {

 protected:
   
  FFTSolver *fft;

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
  *
  *
  * Nx/2+1 ,  we do not know the exact size in advance thus have to
  * included whole size and set to zero.
  **/
  void calcFluxSurfAvrg(CComplex kXOut[Nq][NzLD][Nky][FFTSolver::X_NkxL],
                        CComplex phi_yz[Nx]);

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
  *   solve Eq. in Fourier space, using periodic boundary conditions in X and Y
  *
  **/
  void solveFieldEquations(const CComplex Q     [Nq][NzLD][Nky][NxLD],
                                 CComplex Field0[Nq][NzLD][Nky][NxLD]);

  /**
  *
  *  @brief Solves the gyrokinetic Poisson equation
  * 
  *  The equation to solve is
  *
  *  \f[
  *      \lambda_D^2 \phi + \sum_\sigma frac{q_\sigma^2 n_\sigma}{T_\sigma}
  *       \left( 1 - \Gamma_0(k_\perp^2\rho_{t\sigma}^2) \right) + \frac{q_0^2 n_0}{T_0} 
  *       \left(\phi - \left< \phi \right>_{yz} \right)
  *       = \rho 
  *  \f]
  *  where \f$ \lambda_D \f$ is the Debye length. The part latter part on the LHS
  *  is due to flux-surface averaging effect of an electron-species.
  *
  *  @warning We need to normalize the FFT transform here.
  *
  */
  virtual void solvePoissonEquation(CComplex kXOut[FFTSolver::X_NkxL][Nky][NzLD][Nq],
                                    CComplex kXIn [FFTSolver::X_NkxL][Nky][NzLD][Nq]);
  
  /**
  *
  *  @brief Solves the gyrokinetic Ampere equation for the parallel vector potential \f$ A_{1\parallel} \f$
  *       
  *  The equation to solve is
  *
  *  \f[
  *       \left( \nabla_\perp^2 - 
  *      Y \hat{e} \beta \sum_\sigma \sigma_\sigma \alpha_\sigma^2 q_\sigma \Gamma_0(b_\sigma)  \right)
  *       = j_\parallel 
  *  \f]
  *  where \f$ Y\hat{e} \beta = \tfrac{1}{2}\hat{e} \beta \f$.
  *
  *  Reference : Eq. (2.40) of @cite Dannert_2006:PhDThesis
  *
  *  @warning We need to normalize the FFT transform here.
  *
  **/
  void solveAmpereEquation(CComplex kXOut[FFTSolver::X_NkxL][Nky][NzLD][Nq],
                           CComplex kXIn [FFTSolver::X_NkxL][Nky][NzLD][Nq]);

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
  *      \left( \begin{array}{cc} C_1 & C_2 \\ C_2 & C_3 \end{array} \right)
  *      \cdot
  *      \left( \begin{array}{c} \phi \\ B_\parallel     \end{array} \right)
  *      =
  *      \left( \begin{array}{c} \rho \\ j_\perp         \end{array} \right)
  *  \f]
  * 
  *   It can be decoupled to (see reference for details)
  *
  *  \f[
  *      \left( \begin{array}{c} \phi \\ B_\parallel     \end{array} \right)
  *      =
  *      \frac{1}{C_1 C_3 - C_2^2}
  *      \left( \begin{array}{c} C_3 \rho - C_2 j_\perp \\ C_1 j_\perp - C_2 \phi \end{array} \right)
  *  \f]
  *
  *
  *  Field equations for \f$ \phi \f$ and \f$ B_\parallel \f$ are now decoupled and can be solved
  *  subsequently.
  *    
  *  @warning We need to normalize the FFT transform here.
  *
  */
  void solveBParallelEquation(CComplex kXOut[FFTSolver::X_NkxL][Nky][NzLD][Nq],
                              CComplex kXIn [FFTSolver::X_NkxL][Nky][NzLD][Nq]);

 public:
  
  /**
  *   Constructor
  */ 
  FieldsFFT(Setup *setup, Grid *grid, Parallel *parallel, FileIO * fileIO, Geometry *geo, FFTSolver *fftsolver);
  /**
  *   Destructor
  */ 
 ~FieldsFFT();
  

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
   virtual void gyroAverage(const CComplex In [Nq][NzLD][Nky][NxLD], 
                                  CComplex Out[Nq][NzLD][Nky][NxLD],
                            const int m, const int s, const bool forward, const bool stack=false);
    
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
  *      \left< \phi       \right> &= exp(-k_\perp^2 \rho_{t\sigma} ) \phi       (k_x, k_y, z) \quad, \\
  *      \left< j_parallel \right> &= exp(-k_\perp^2 \rho_{t\sigma} ) j_\parallel(k_x, k_y, z) \quad, \\
  *  \f}
  *  
  *  Note : Not sure if this gyro-averaging is valid for \f$ j_\parallel \f$ and how to do for 
  *         \f$ B_{1\parallel} \f$.
  *
  *  @param fields    the field or source terms
  *  @param s         index to species
  *  @param nField    index of field (phi, jp, jo) not used
  *  @param gyroField forward- or backward transformation
  *
  */
  void gyroFirst(const CComplex   In [Nq][NzLD][Nky][NxLD], 
                       CComplex   Out[Nq][NzLD][Nky][NxLD],
                       CComplex kXOut[Nq][NzLD][Nky][FFTSolver::X_NkxL],
                       CComplex kXIn [Nq][NzLD][Nky][FFTSolver::X_NkxL],
                 const int s, const bool gyroField=false) ; 
  
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
  *
  *  @todo rename gyro-field to something better, e.g isForwardAverage
  *
  *  @param fields    the field or source terms
  *  @param m         index of magnetic moment
  *  @param s         index to species
  *  @param nField    index of field (phi, jp, jo) not used
  *  @param gyroField forward- or backward transformation
  *
  */
  void gyroFull(const CComplex In [Nq][NzLD][Nky][NxLD], 
                      CComplex Out[Nq][NzLD][Nky][NxLD],
                      CComplex kXOut[Nq][NzLD][Nky][FFTSolver::X_NkxL],
                      CComplex kXIn [Nq][NzLD][Nky][FFTSolver::X_NkxL],
                const int m, const int s, bool stack=false);  

  /**
  *   @brief performs the double gyro-averaging over Maxwellian in Fourier space
  * 
  *   Note, that in the local, periodic case, following indentities hold
  *
  *  \f{align}{
  *         
  *    \int_0^\infty J_0^2(\sqrt{(lambda mu)}) e^{-\mu}  d\mu  
  *                &= I_0(b) e^{-b} \\
  *    \int_0^\infty J_0^2(\sqrt{(lambda mu)}) e^{-\mu}\mu d\mu
  *                &= \Gamma_0-b\left(\Gamma_0-\Gamma_1\right) \\
  *    \int_0^\infty J_0^2(\sqrt{(lambda mu)}) e^{-\mu}\mu^{1/2} 
  *                &= ? \\
  *    \int_0^\infty J_0^2(\sqrt{(lambda mu)}) e^{-\mu}\mu^2 d\mu
  *                &= ? 
  *
  *  \f}
  *
  *  for, first and second equation, see Callen.
  *
  **/  
  void doubleGyroExp(const CComplex In [Nq][NzLD][Nky][NxLD], 
                           CComplex Out[Nq][NzLD][Nky][NxLD], const int m, const int s);
  
  /**
  *    @brief calculates the field energy
  * 
  *    @todo according to what ? Cite idomura
  **/ 
  void getFieldEnergy(double& phiEnergy, double& ApEnergy, double& BpEnergy);
  
protected:

  /**
  *   @brief Print out some runtime information
  */ 
  virtual void printOn(std::ostream &output) const;


  /**
  *  @brief helper function for 
  *  \f[
  *     r = \sum_\sigma \frac{q_\sigma^2 n_{0\sigma}}{T_{0\sigma}} (1 - \Gamma_0(b_\sigma)) 
  *  \f]
  *  
  *
  *  @note is it worth to precalculate this values ?
  *
  *  @param  k2_p \f$ k_\perp^2 \f$ (not normalized) perpendicular wavenumber
  *  @return      \f$ r \f$
  *
  **/
  inline double sum_qqnT_1mG0(const double k2_p) ;

  /**
  *  @brief helper function for 
  *  
  *  \f[
  *     r = \sum_\sigma \sigma_\sigma alpha_\sigma^2 q_\sigma \Gamma_0(b_\sigma) 
  *  \f]
  *
  *  @note is it worth to pre-calculate this values ?
  *
  *  Note :
  *
  *    In case of drift-kinetic 1-(1-Gamma_0)) becomes 1 thus
  *  
  *  \f[
  *     r = \sum_\sigma \sigma_\sigma alpha_\sigma^2 q_\sigma 
  *  \f]
  *
  *  @param   k2_p \f$ k_\perp^2 \f$ (not normalized) perpendicular wavenumber
  *  @return       \f$ r \f$
  *
  **/
  inline double sum_sa2qG0(const double k2_p) ;

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
  *  @param  k2_p  \f$ k_\perp^2 \f$ (not normalized) perpendicular wavenumber
  *  @return      \f$ r \f$
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
  *  @param  k2_p \f$ k_\perp^2 \f$ (not normalized) perpendicular wavenumber
  *  @return     \f$ r \f$
  **/
  inline  double sum_2TnBB_Delta(const double k2_p) ;


};


#endif // __FIELDS_FFT_H
