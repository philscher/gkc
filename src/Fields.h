/*
 * =====================================================================================
 *
 *       Filename: Fields.h
 *
 *    Description: Main Interface for solving fields equations. Implements
 *                 source terms calculations.
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

namespace Field {  const int phi=1, Ap=2, Bp=3, Bpp=4; }
namespace Q     {  const int rho=1, jp=2, jo=3; }
/**
 *
 *  Poisson base
 *
 *  Governs the calculation of the gyro-averaged potential, which
 *  should be same to all solution methods
 *
 *  basically it does
 *
 *  1. caluclates the gyro-averaged charge density
 *  2. solves the Poisson equation
 *  3. calculates the gyro-averaged electric potential
 *  4. calculates the gyro-averaged electric fields.
 * 
 *
 *
 *  Field equations for $\phi\$ and $B_\parallel$ are coupled. We can decouple
 *  them by the method of 
 *  Original Reference : Merz, PhD Thesis
 *  Backgroun Informations : wikipedia
 *
 * Define : 
 * \f[
 *     C_1 = \\
 *     C_2 = \\
 *     C_3 = '' 
 *    
 *     \hat{rho}     = \epsilon_\sigma n_{0\sigma} \pi q_sigma B_{0i} \int Gyro<g_j> d\mu dv_\parallel \\ 
 *     \hat{j}_\perp = \epsilon_\sigma n_{0\sigma} \pi q_sigma B_{0i} \int \sqrt{\mu} J_1(\lambda_sigma) g_j d\mu dv_\parallel \\ 
 *
 *  \f]
 *  NOTE: \hat{j} can be more effectively caluclated in Fourier space, see reference there..
 *  and 
 *  \f[
 *      \phi          = ... \\
 *      B_\parallel  = ...
 *  \f]
 *
 *   write as, where  
 *
 *   \f[
 * %        \left[  array{ll} { C_1 & C_2 \\ C_2 & C_3  } \right] 
 *         \left[ \array{l}  {\phi \\ B_\parallel  }   \right]
 *      =  \left[ \array{l} { <\rho> \\ <j_\perp>    } \right]
 *   \f]
 *   multiplying the inverse of C we get the following equations, which are now decoupled
 *
 *   \f[
 *      \left[ \array{l} { \phi \\ B_\parallel } \right]  = 
 *      \left[  array{ll} { C_1 & C_2 \\ C_2 & C_3  } \right] \cdot  \left[ \array{ll} { <\rho> \\ <j_\perp>    } \right]
 *   \f]
 *    we can now solve for $\phi$ and $B_\parallel$ subsequently.
 **/
class Fields : public IfaceHelios {
  
  Array6z  SendYu, SendXu, SendYl, SendXl, SendZu, SendZl; 
  Array6z  RecvYu, RecvXu, RecvYl, RecvXl, RecvZu, RecvZl;

protected:
  Grid *grid;
  Parallel *parallel;
  Geometry<HELIOS_GEOMETRY> *geo;
  
  /** Solve Equation system */
  int solveEq;
/*! Performs the gyro-average for variable stored in FFT->r3In
 *  
 *
 */
  Array1z  phi_yz;

/** 
 *  Calculation of the ITG potential.
 * solve the averaged Poisson equation over y-z plane only for FFT
 * n this should be OK. For slab geometry this is k(kx,0,0) !
 * Result is given in Fourier space which is only k_x dependent. For solving
 * the fields equation, we need to substract it from the adiabatic repsonse, but
 * only for(ky=0, kz=0) modes !!!.
 **/
virtual Array1z calcFluxSurfAvrg(Array4z n) = 0;

public:
/** 
 *  Calculates the charge density, by performing gyro-average in Fourier space.
 *  This function works in 3 steps.
 *
 *
 *  \f[   \left< rho \right>_G = \int_v \sum_\sigma \alpha_sigma \hat{q} \pi \hat{B}  \f]
 *
 *
 *
 * */
virtual Array3z calculateChargeDensity         (Array6z f0,  Array6z f);


/**
 * Calcultaes the current density
 *
 * Calculates
 *  \f[   \left< j_\parallel \right>_G = \int_v \sum_\sigma \alpha_sigma \hat{q} \pi \hat{B}
 *       int_m \widetilde{ int v_\parallel g_j1 dv_\parallel} d\mu
 *  \f]
 *
 */
virtual Array3z calculateParallelCurrentDensity(Array6z f0, Array6z f );

/**
 *  Calculates the perpendicular current density
 *
 *  Calculates
 *  \f[
 *          \left< j_\perp  \right> = \int  \mu I_1(\lambda_\sigma) f_\sigma d^3 v 
 *  \f]
 *
 *  Reference :
 *
 *  Note : Defined virtual, as there may be more effective implementations
 */
virtual Array3z calculatePerpendicularCurrentDensity(Array6z f0, Array6z f );




virtual Array4z solveFieldEquations(Array4z fields, Timing timing) = 0;

/*!
 *
 *  Solver the gyrokinetic Poisson's equation :
 *  
 *  This is the decoupled version and includes \f[ \beta_\parallel \f] effects if enabled
 *  for high-beta plasmas.
 * 
 */
virtual Array3z solvePoissonEquation(Array3z rho, Timing timing)= 0;
/*!
 *
 *  Solver the gyrokinetic Poisson's equation :
 *  
 *  This is the decoupled version and includes \f[ \beta_\parallel \f] effects if enabled
 *  for high-beta plasmas.
 * 
 */
virtual Array3z solveAmpereEquation(Array3z j, Timing timing) = 0;
/*!
 *
 *  Solver the gyrokinetic Poisson's equation :
 *  
 *  This is the decoupled version and includes \f[ \beta_\parallel \f] effects if enabled
 *  for high-beta plasmas.
 * 
 */
virtual Array3z solveBParallelEquation(Array3z j, Timing timing) = 0;




public:
  
  /** Sets the boundary
   *
   *
   *
   *
   */ 
  int setBoundary(Array6z  A);
   
  
  std::string  ApPerturbationStr;
  virtual Array3z gyroAverage(Array3z fields, int m, int s, int nField, bool gyroField=false) = 0;
  //* Electric potential and Electric fields
  Array4z  Q, Field0;
  Array6z  Field;
  Array5z  phi, Ap, Bp;
  /** BlabBla
   *
   *
       //  brackets should be 1/2 but due to numerical error we should include calucalte ourselves, see Dannert[2] 
       Yeb = (1./sqrt(M_PI) * sum(pow2(V) * exp(-pow2(V))) * dv) * geo->eps_hat * plasma->beta; 
   * 
   *  \f[ Y_{analyitic} = \frac{1}{\sqrt{pi}} \int_-infty^infty_ v^2 e^{-v^2} dv \equiv \frac{1}{2} \f] 
   *
   *   but due to discretization errors, this is not exactly
   *
   *   \f Yeb = Y \hat{\eps}_{geo} \beta  \f]
   *
   *   see Dannert
   */

  double Yeb;
  //* Constructor
  Fields(Setup *setup, Grid *grid, Parallel *parallel, FileIO *fileIO, Geometry<HELIOS_GEOMETRY> *geo);
  //* Destructor
  virtual ~Fields();
   
  std::string gyroAverageModel;
  
  int solve(Array6z f0, Array6z  f, Timing timing=0, int rk_step=0);

  FileAttr *FA_phi, *FA_Ap, *FA_Bp, *FA_phiTime;

     
  void initDataOutput(Setup *setup, FileIO *fileIO);
  virtual void writeData(Timing timing, double dt);
  
  void closeData();
  Timing dataOutputPhi;

  int getSolveEq() const { return solveEq; };
protected:

        virtual void printOn(ostream &output) const {
         output   << "Poisson    |  " << "Base class" << std::endl;
         output   << "Ampere     |  " << ((plasma->nfields >= 2) ? "beta :  " + Num2String(plasma->beta) : " --- no electromagnetic effects ---") << std::endl;
         output   << "B_parallel |  " << ((plasma->nfields >= 3) ? "beta :  " + Num2String(plasma->beta) : " --- no electromagnetic effects ---") << std::endl;
         output   << "           |  Pert. : " << ((!(solveEq & Field::Ap) && (plasma->nfields >= 2))     ? ApPerturbationStr  : " ") << std::endl;
         output   << "           | Gyro-Average model : " << gyroAverageModel << std::endl;
        }
};


#endif // __FIELDS_H
