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
 *  @file Fields.h
 *
 *  @brief 
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
 *     C_3 = 
 *    
 *     \hat{rho}     = \epsilon_\sigma n_{0\sigma} \pi q_\sigma B_{0i} \int Gyro<g_j> d\mu dv_\parallel 
 *     \hat{j}_\perp = \epsilon_\sigma n_{0\sigma} \pi q_\sigma B_{0i} \int \sqrt{\mu} J_1(\lambda_\sigma) g_j d\mu dv_\parallel 
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
 *         \left[ \begin{array}{ll}  C_1 & C_2 \\ C_2 & C_3    \end{array} \right] \cdot 
 *         \left[ \begin{array}{l}  \phi       \\ B_\parallel  \end{array} \right]
 *      =  \left[ \begin{array}{l}   <\rho>    \\ <j_\perp>    \end{array} \right]
 *   \f]
 *   multiplying the inverse of C we get the following equations, which are now decoupled
 *
 *   \f[
 *      \left[ \begin{array}{l}   \phi \\ B_\parallel      \end{array} \right]  = 
 *      \left[ \begin{array}{ll}  C_1 & C_2 \\ C_2 & C_3   \end{array} \right] \cdot
        \left[ \begin{array}{ll}  <\rho> \\ <j_\perp>      \end{array} \right]
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
 *  \f[   \left< rho \right>_G = \int_v \sum_\sigma \alpha_\sigma \hat{q} \pi \hat{B}  \f]
 *
 *  @param f0 The Maxwellian phase-space background distribution
 *  @param f  Currect phase-space distribution
 *  @param m  index for perpendcular velocity \f[ \mu = M(m)  \f]
 *  @param s  index for species
 *
 * */
virtual Array3z calculateChargeDensity         (Array6z f0,  Array6z f, const int m, const int s);


/**
 * Calcultaes the current density
 *
 * Calculates
 *  \f[   \left< j_\parallel \right> = \int_v \sum_\sigma \alpha_\sigma \hat{q} \pi \hat{B}
 *        int_m \widetilde{ int v_\parallel g_{j1} dv_\parallel} d\mu
 *  \f]
 *
 *  @param f0 The Maxwellian phase-space background distribution
 *  @param f  Currect phase-space distribution
 *  @param m  index for perpendcular velocity \f[ \mu = M(m)  \f]
 *  @param s  index for species
 *
 */
virtual Array3z calculateParallelCurrentDensity(Array6z f0, Array6z f, const int m, const int s );

/**
 *  Calculates the perpendicular current density
 *
 *  Calculates
 *  \f[
 *      \left< j_\perp \right> = \int  \mu I_1(\lambda_\sigma) f_\sigma d^3 v 
 *  \f]
 *
 *  Reference :
 *
 *  Note : Defined virtual, as there may be more effective implementations
 *
 *  @param f0 The Maxwellian phase-space background distribution
 *  @param f  Currect phase-space distribution
 *  @param m  index for perpendcular velocity \f[ \mu = M(m)  \f]
 *  @param s  index for species
 */
virtual Array3z calculatePerpendicularCurrentDensity(Array6z f0, Array6z f, const int m, const int s );




/*!
 *  @brief Solve the field equation for the current phase-space distribution  
 *  This is the decoupled version and includes \f[ \beta_\parallel \f] effects if enabled
 *  for high-beta plasmas.
 *
 *  @param Q current quell terms
 * 
 */
virtual Array4z solveFieldEquations(Array4z Q, Timing timing) = 0;

/*!
 *
 *  This is the decoupled version and includes \f[ \beta_\parallel \f] effects if enabled
 *  for high-beta plasmas.
 * 
 */
virtual Array3z solvePoissonEquation(Array3z rho, Timing timing)= 0;
/*!
 *
 *  This is the decoupled version and includes \f[ \beta_\parallel \f] effects if enabled
 *  for high-beta plasmas.
 * 
 */
virtual Array3z solveAmpereEquation(Array3z j, Timing timing) = 0;
/*!
 *
 *  This is the decoupled version and includes \f[ \beta_\parallel \f] effects if enabled
 *  for high-beta plasmas.
 * 
 */
virtual Array3z solveBParallelEquation(Array3z j, Timing timing) = 0;




public:
  
  /** Sets the boundary
   *  @brief  update boundaries for the gyro-averaged field which arises due to domain
   *          decomposition
   *  @param  Current gyro-averaged field  \f[ A(x,y,z,\mu, \sigma) \f] 
   */ 
  int setBoundary(Array6z  A);
   
  
  std::string  ApPerturbationStr;
  virtual Array4z gyroAverage(Array4z fields, int m, int s, int nField, bool gyroField=false) = 0;
  //* Electric potential and Electric fields
  Array4z  Q, Field0;
  Array6z  Field;
  Array5z  phi, Ap, Bp;
  /**
   *
   *
   *    //  brackets should be 1/2 but due to numerical error we should include calucalte ourselves, see Dannert[2] 
   *    Yeb = (1./sqrt(M_PI) * sum(pow2(V) * exp(-pow2(V))) * dv) * geo->eps_hat * plasma->beta; 
   * 
   *  \f[ Y = \frac{1}{\sqrt{\pi}} \int_{-\infty}^\infty v_\parallel^2 e^{-v_\parallel^2} dv \equiv \frac{1}{2} \f] 
   *
   *   but due to discretization errors, this is not exact, thus needs to be calculated numerically.
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
  
  /**
   *  @brief solve the field equations
   *
   *
   */
  int solve(Array6z f0, Array6z  f, Timing timing=0, int rk_step=0);


  int getSolveEq() const { return solveEq; };
  /*
   *  @brief write field values put to data file
   * 
   *  
   */ 
  virtual void writeData(Timing timing, double dt);
protected:
  	FileAttr *FA_phi, *FA_Ap, *FA_Bp, *FA_phiTime;
  	Timing dataOutputFields;
  	void initDataOutput(Setup *setup, FileIO *fileIO);
        void closeData();

        virtual void printOn(ostream &output) const {
         output   << "Poisson    |  " << "Base class" << std::endl;
         output   << "Ampere     |  " << ((plasma->nfields >= 2) ? "beta :  " + Num2String(plasma->beta) : " --- no electromagnetic effects ---") << std::endl;
         output   << "B_parallel |  " << ((plasma->nfields >= 3) ? "beta :  " + Num2String(plasma->beta) : " --- no electromagnetic effects ---") << std::endl;
         output   << "           |  Pert. : " << ((!(solveEq & Field::Ap) && (plasma->nfields >= 2))     ? ApPerturbationStr  : " ") << std::endl;
        }
};


#endif // __FIELDS_H
