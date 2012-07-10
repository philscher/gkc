/*
 * =====================================================================================
 *
 *       Filename: Fields.h
 *
 *    Description: Main Interface for solving fields equations. 
 *                 Implements source terms calculations.
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __FIELDS_H
#define __FIELDS_H

#include "Global.h"

#include "Setup.h"
#include "Parallel.h"
#include "Grid.h"
#include "Geometry.h"
#include "GeometrySlab.h"
#include "GeometryShear.h"
#include "Geometry2D.h"

#include "FileIO.h"

#include "Timing.h"
#include "Plasma.h"

/**
*  @brief defined offset to access field variables
*
*         \f$ \phi          \f$ electric potential
*         \f$ A_{1\parallel}\f$ parallel magnetic vector potential 
*         \f$ B_{1\parallel}\f$ parallel magnetic field
*/
namespace Field { const int phi=1, Ap=2, Bp=3, Bpp=4; }

/**
*  @brief defined offset to access source term variables S
*         \f$ \rho, j_\parallel, j_\perp \f$
*         with 
*         \f$ \rho        \f$ charge density
*         \f$ j_\parallel \f$ parallel current density
*         \f$ j_\perp     \f$ perpendicular current density
*/
namespace Q     { const int rho=1, jp=2, jo=3; }

/**
*
*  @file Fields.h
*
*  @brief 
*
*  Governs the calculation of the gyro-averaged potentials 
*  \f$ <\phi> \f$, \f$<A_\parallel>\f$, \f$ \left<B_\perp\right> \f$
*  from the phase space function \f[ f_{1\sigma} \f].
*
*  The calculation is split up into 4 steps.
*
*  1. caluclates the gyro-averaged source density (integration over v)
*  2. Back-transformation
*  3. Solves the field equation
*  4. Forward-transformation of the field equation.
*
*  Calculating the field terms are provided as pure virtual function.
*  Thus the solution of the corresponding Laplace type equation needs
*  to be implmented using the method of choice.
*
*/
class Fields : public IfaceHelios {
  
   /**
   *  @brief buffers for field equations.
   *
   *  @todo : current version does not need to communicate Y
   *          direction. Thus remove. 
   */
   Array6z  SendXu, SendYu, SendZu, SendXl, SendYl, SendZl; 
   /**
   *  @brief please document me
   */
   Array6z  RecvXu, RecvYu, RecvZu, RecvXl, RecvYl, RecvZl;

protected:

   /**
   *  @brief please document me
   */
   Grid *grid;
   
   /**
   *  @brief please document me
   */
   Parallel *parallel;
   
   /**
   *  @brief Geometry interface
   */
   Geometry<HELIOS_GEOMETRY> *geo;
  
   /**
   *  @brief please document me
   */
   int solveEq;

   /** 
   *  @brief Performs the gyro-average for variable stored in FFT->r3In
   */
   Array1z  phi_yz;

   /** 
   *  Calculation of the ITG potential.
   * solve the averaged Poisson equation over y-z plane only for FFT
   * n this should be OK. For slab geometry this is k(kx,0,0) !
   * Result is given in Fourier space which is only k_x dependent. For solving
   * the fields equation, we need to substract it from the adiabatic repsonse, but
   * only for(ky=0, kz=0) modes !!!.
   *
   * @todo : This is implementation dependent. Does need to be as an extra function.
   */
   virtual Array1z calcFluxSurfAvrg(Array4z n) = 0;

public:
   /** 
   *  @brief Calculates the charge density for \f$ \rho(x,y_k,z; \mu, \sigma) \f$.
   *
   *  The charge density is caluclated according to
   *
   *  \f[   
   *   \rho(x,y_k,z;\mu,\sigma) = q_\sigma \int_{v_\parallel}  g_{1\sigma}(x,y_k,z,v_\parallel,\mu,\sigma) \textrm{d}alpha 
   *  \f]
   *  with \f$\textrm{d}\alpha=n_{0;\sigma} \pi \hat{B}_0 \textrm{d}v_\parallel \textrm{d}\mu\f$.
   *
   *  @param f0 The Maxwellian phase-space background distribution
   *  @param f  Currect phase-space distribution
   *  @param m  index for perpendcular velocity \f$ \mu = M(m)  \f$
   *  @param s  index for species
   *
   */
   virtual Array3z calculateChargeDensity         (Array6z f0,  Array6z f, const int m, const int s);

   /**
   * @brief Calculates the parallel current density \f$ j_\parallel(x,y_k,z;\mu,\sigma) \f$.
   *
   * The parallel current density is calculated according to 
   *
   *  \f[   
   *       j_\parallel(x,y_k,z;\mu,\sigma) = q_\sigma \int_{v_\parallel} g_{1\sigma}(x,y_k,z;\mu,\sigma) \textrm{d}\beta
   *  \f]
   *  with \f$\textrm{d}\beta=n_{0;\sigma} \alpha_\sigma \pi \hat{B}_0 \textrm{d}v_\parallel \textrm{d}\mu\f$.
   * 
   *  Reference : 
   *  
   *  @param f0 The Maxwellian phase-space background distribution
   *  @param f  Currect phase-space distribution
   *  @param m  index for perpendcular velocity \f$ \mu = M(m)\f$
   *  @param s  index for species
   *
   */
   virtual Array3z calculateParallelCurrentDensity(Array6z f0, Array6z f, const int m, const int s );

   /**
   *  Calculates the perpendicular current density \f$ j_\perp(x,y_k,z;\mu,\sigma) \f$
   *
   *  Calculates
   *  \f[
   *      j_\perp(x,y_k,z;\mu,\sigma) = q_\sigma \int_{v_\parallel=-\infty}^\infty \mu g_{1\sigma} \textrm{d}\gamma
   *  \f]
   *  with \f$\textrm{d}\gamma=n_{0;\sigma} \alpha_\sigma \pi B_0 \textrm{d}v_\parallel \textrm{d}\mu \f$.
   *
   *  See Equation ?? @cite DannertPhD
   *
   *  @note Defined virtual, as there may be more effective implementations
   *
   *  @param f0 The Maxwellian phase-space background distribution
   *  @param f  Currect phase-space distribution
   *  @param m  index for perpendcular velocity \f$ \mu = M(m)\f$
   *  @param s  index for species
   */
   virtual Array3z calculatePerpendicularCurrentDensity(Array6z f0, Array6z f, const int m, const int s );


   /**
   *  @brief Solves the field equation from the source terms.
   *  
   *  @param Q current source terms
   * 
   */
   virtual Array4z solveFieldEquations(Array4z Q, Timing timing) = 0;

   /**
   *    @brief virtual function to solve the Poisson equation for \f$ \phi \f$
   * 
   */
   virtual Array3z solvePoissonEquation(Array3z rho, Timing timing)= 0;

   /**
   *
   *    @brief virtual function to solve the Ampere equation for \f$ A_{1\parallel} \f$
   *  
   */
   virtual Array3z solveAmpereEquation(Array3z j, Timing timing) = 0;

   /**
   *    @brief virtual function to solve the Ampere equation for \f$ B_{1\parallel} \f$
   */
   virtual Array3z solveBParallelEquation(Array3z j, Timing timing) = 0;

public:
  
   /** Sets the boundary
   *  @brief  update boundaries for the gyro-averaged field which arises due to domain
   *          decomposition
   *  @param  Current gyro-averaged field  \f[ A(x,y,z,\mu, \sigma) \f] 
   */ 
   int setBoundary(Array6z  A);
   
   /**
   *  @brief please document me
   */
   std::string  ApPerturbationStr;
  
   /**
   *  @brief performs gyro-averaging (back and forward) of the fields.
   *
   *  \f[
   *        \begin{array}{l}
   *        \left< \phi        \right>_{\mu\sigma} \\
   *        \left< A_\parallel \right>_{\mu\sigma} \\ 
   *        \left< B_\perp     \right>_{\mu\sigma} \\
   *        \end{array}
   *        =
   *        \mathcal{G}_{\mu\sigma}{\left( 
   *        \begin{array}{l}
   *         \phi \\
   *         A_\parallel\\
   *         B_\perp  
   *        \end{array} \right) }
   *  \f]
   *  where \f$ \mathcal{G}_{\mu\sigma} \f$ the gyro-averaging operator, for magnetic
   *  moment \f$ \mu \f$ and species \f$ \sigma \f$.
   *
   *
   *  @param fields    the corresponding fields as 4 dimensional array.
   *  @param m         the corresponding index to the magnetic moment
   *  @param s         for species \f$ \sigma \f$
   *  @param nField    specify which field (not used)
   *  @param gyroField true if forward transformation, false if backward transformation
   *
   */
   virtual Array4z gyroAverage(Array4z fields, int m, int s, int nField, bool gyroField=false) = 0;
  
   /**
   *  @brief four-dimensional array hold the source terms in drift-coordinates.
   *
   *  @Note  Q stand for the German Quellterme (translated : source terms).
   */
   Array4z  Q;
   /**
   *  @brief four-dimensional array hold the field terms in drift-coordinates.
   */
   Array4z  Field0;
   
   /**
   *  @brief gyro-averaged field quantities as \f$ \left<\phi \right, <A_\parallel>, <B_\perp>)\f$
   */
   Array6z  Field;
   
   /**
   *  @brief please document me
   */
   Array5z  phi, Ap, Bp;
   
   /**
   * 
   *  @brief \f$ \tfrac{1}{2} \hat{e} \hat{B} \f$ normalization constants 
   *
   *  Brackets should be 1/2 but due to numerical error we should include calucalte ourselves, see Dannert[2] 
   *  Yeb = (1./sqrt(M_PI) * sum(pow2(V) * exp(-pow2(V))) * dv) * geo->eps_hat * plasma->beta; 
   * 
   *  \f[ Y = \frac{1}{\sqrt{\pi}} \int_{-\infty}^\infty v_\parallel^2 e^{-v_\parallel^2} dv \equiv \frac{1}{2} \f] 
   *
   *   but due to discretization errors, this is not exact, thus needs to be calculated numerically.
   *
   *   \f Yeb = Y \hat{\eps}_{geo} \beta  \f]
   *
   *   @ref see Dannert (2008), PhD Thesis], where which chapter ... ?
   */
   double Yeb;
  
   /**
   *
   * Constructor
   *
   */
   Fields(Setup *setup, Grid *grid, Parallel *parallel, FileIO *fileIO, Geometry<HELIOS_GEOMETRY> *geo);
  
   /**
   *
   *   Destructor
   *
   */
   virtual ~Fields();
  
   /**
   *  @brief solve the field equations
   *
   *
   */
   int solve(Array6z f0, Array6z  f, Timing timing=0, int rk_step=0);

   /**
   *  @brief set which equation to solve if one field is assumed to be fixed
   *
   */
   int getSolveEq() const { return solveEq; };

   /**
   *  @brief write field values put to data file
   * 
   *  Thus function saves the state of all fields at the current timestep.
   *  you are using, which is either \f$ (\phi) \f$, \f$ (\phi, A_\parallel) \f$,
   *  \f$ (\phi, A_\parallel, B_\parallel) \f$. Values are stored at
   *  
   *  Can be accessed by HDF-5 using (assuming fileh5 is your the file) 
   *   
   *  fileh5.root.Fields.phi[z, y_k, x, frame]
   *  fileh5.root.Fields.Ap [z, y_k, x, frame]
   *  fileh5.root.Fields.Bp [z, y_k, x, frame]
   *  fileh5.root.Fields.Timing[frame]
   *
   */ 
   virtual void writeData(Timing timing, double dt);

protected:

   /**
   *  @brief please document me
   */
  	FileAttr *FA_phi, *FA_Ap, *FA_Bp, *FA_phiTime;

   /**
   *  @brief please document me
   */
  	Timing dataOutputFields;

   /**
   *  @brief please document me
   */
  	void initDataOutput(Setup *setup, FileIO *fileIO);
   
   /**
   *  @brief please document me
   */
   void closeData();

   /**
   *  @brief please document me
   */
   virtual void printOn(ostream &output) const ;
};


#endif // __FIELDS_H
