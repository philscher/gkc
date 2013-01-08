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
#include "Parallel/Parallel.h"
#include "Grid.h"
#include "Geometry.h"
#include "FileIO.h"

#include "Timing.h"
#include "Plasma.h"


extern int GC2, ///< Number of boundary points (needed ?)
           GC4, ///< Number of extended boundary points (needed?)
            Nq; ///< Number of fields

/**
*  @brief defined offset to access field variables
*
*         \f$ \phi          \f$ electric potential
*         \f$ A_{1\parallel}\f$ parallel magnetic vector potential 
*         \f$ B_{1\parallel}\f$ parallel magnetic field
**/
namespace Field { const int  phi =0 , ///< electric field 
                             Ap  =1 , ///< parallel magnetic vector potential component
                             Bp  =2 , ///< parallel magnetic field component
                             Iphi=4 , ///< electric field (logic) index
                             IAp =8 , ///< parallel magnetic potential component (logic) index
                             IBp =16; ///< parallel magnetic field component     (logic) index
                };
//namespace Field { constexpr const int phi=1, Ap=2, Bp=3, Bpp=4; }

/**
*
*  @brief Calculate sources and interface to field solvers
*
*  Governs the calculation of the gyro-averaged potentials 
*  \f$ <\phi> \f$, \f$<A_\parallel>\f$, \f$ \left<B_\perp\right> \f$
*  from the phase space function \f[ f_{1\sigma} \f].
*
*  The calculation is split up into 4 steps.
*
*  1. calculates the gyro-averaged source density (integration over v)
*  2. Back-transformation
*  3. Solves the field equation
*  4. Forward-transformation of the field equation.
*
*  Calculating the field terms are provided as pure virtual function.
*  Thus the solution of the corresponding Laplace type equation needs
*  to be implemented using the method of choice.
*
*  Integration
*          
*           The integration in \f$ v_\parallel \f$ is performed using
*           the Trapezoidal rule, which due to 
*           \f$f_1(-L_v) \approx f_1(L_v) \approx 0 \f$ should give
*           approximately 2nd order integration. 
*
*           @note Alternative modified Simpson rule can also be used, 
*                 however, no benefit is found. Re-confirm.
*
*
*           The integration in \f$ \mu \f$ is performed using either
*           Gaussian-Legendre integration or Trapezoidal rule.
*           @note : is Gauss-Laguerre not the better choice ? Check. 
*           The perpendicular interation can be setup using 
*           Integration.cpp (move to grid).
*   
*
**/
class Fields : public IfaceGKC {
  
  /**
  *  @brief buffers for field equations.
  *
  *  @todo : current version does not need to communicate Y
  *          direction. Thus remove.
  *
  *  Send to down in our hypercube means from 
  *  X = Coord[DIR_X] -> Coord[DIR_X] + 1
  *
  **/
  CComplex *SendXl,  ///< Send Buffer in X-direction (send to down)
           *SendXu,  ///< Send buffer in X-direction (send to up)
           *SendZl,  ///< Send buffer in Z-direction (send to down)
           *SendZu,  ///< Send buffer in Z-direction (send to up) 
           *RecvXl,  ///< Recv buffer in X-direction (recv from down)
           *RecvXu,  ///< Recv buffer in X-direction (recv from up)
           *RecvZl,  ///< Recv buffer in Z-direction (recv from down) 
           *RecvZu;  ///< Recv buffer in Z-direction (recv from up)
   
  nct::allocate ArrayBoundX, ///< Array class for Boundary in X
                ArrayBoundZ; ///< Array class for Boundary in Z

 protected:

  /**
  *  @brief please document me
  **/
  Grid *grid;
   
  /**
  *  @brief please document me
  **/
  Parallel *parallel;
   
  /**
  *  @brief Geometry interface
  **/
  Geometry *geo;
  
  /**
  *  @brief please document me
  **/
  int solveEq;


 public:
  
  /**
  *
  *   @brief Constructor
  *
  *   Accepts following setup parameters ....
  *
  **/
  Fields(Setup *setup, Grid *grid, Parallel *parallel, FileIO *fileIO, Geometry *geo);
  

  /**
  *
  *   Destructor
  *
  **/
  virtual ~Fields();

  /** 
  *  @brief Calculates the charge density for \f$ \rho(x,y_k,z; \mu, \sigma) \f$.
  *
  *  The charge density is calculated according to
  *
  *  \f[   
  *   \rho(x,y_k,z;\mu,\sigma) = q_\sigma \int_{v_\parallel}  g_{1\sigma}(x,y_k,z,v_\parallel,\mu,\sigma) \textrm{d}alpha 
  *  \f]
  *  with \f$\textrm{d}\alpha=n_{0;\sigma} \pi \hat{B}_0 \textrm{d}v_\parallel \textrm{d}\mu\f$.
  *
  *  @param f0 The Maxwellian phase-space background distribution
  *  @param f  Current phase-space distribution
  *  @param m  index for perpendicular velocity \f$ \mu = M[m]  \f$
  *  @param s  index for species
  *
  **/
  void calculateChargeDensity(const CComplex f0 [NsLD][NmLD][NzLB][Nky][NxLB][NvLB],
                              const CComplex f  [NsLD][NmLD][NzLB][Nky][NxLB][NvLB],
                              CComplex Field0           [Nq][NzLD][Nky][NxLD]      ,
                              const int m, const int s) ;

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
  *  @param m  index for perpendicular velocity \f$ \mu = M[m]\f$
  *  @param s  index for species
  *
  **/
  void calculateParallelCurrentDensity(const CComplex f0 [NsLD][NmLD][NzLB][Nky][NxLB][NvLB],
                                       const CComplex f  [NsLD][NmLD][NzLB][Nky][NxLB][NvLB],
                                       CComplex Field0           [Nq][NzLD][Nky][NxLD]      ,
                                       const int m, const int s) ;

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
  *  @param f  Current phase-space distribution
  *  @param Q  Current phase-space distribution
  *  @param m  index for perpendicular velocity \f$ \mu = M[m]\f$
  *  @param s  index for species
  *
  **/
  virtual void  calculatePerpendicularCurrentDensity(const CComplex f0 [NsLD][NmLD][NzLB][Nky][NxLB][NvLB],
                                                     const CComplex f  [NsLD][NmLD][NzLB][Nky][NxLB][NvLB],
                                                     CComplex Field0           [Nq][NzLD][Nky][NxLD]      ,
                                                     const int m, const int s); 

  /**
  *  @brief Solves the field equation from the source terms.
  *  
  *  @param Q current source terms
  *
  *  @note pure virtual function
  * 
  **/
  virtual void solveFieldEquations(const CComplex Q     [Nq][NzLD][Nky][NxLD],
                                         CComplex Field0[Nq][NzLD][Nky][NxLD]) = 0;

  /**
  *    @brief calculate Field energy
  *
  *    calculates The field energy
  *
  *    @out phiEnergy  Total Energy of electric field
  *    @out ApEnergy   total energy of perpendicular magnetic fieled (perturbation)
  *    @out BpEnergy   total energy of parallell magnetic field (perturbation)
  **/
  virtual void getFieldEnergy(double& phiEnergy, double& ApEnergy, double& BpEnergy) = 0;

  /** Sets the boundary
  *  @brief  update boundaries for the gyro-averaged field which arises due to domain
  *          decomposition
  **/ 
  void updateBoundary(); 
   
  /**
  *  @brief please document me
  **/
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
  **/
  virtual void gyroAverage(const CComplex In [Nq][NzLD][Nky][NxLD], 
                                 CComplex Out[Nq][NzLD][Nky][NxLD],
                           const int m, const int s, const bool forward, const bool stack=false) = 0;
  
   
  /**
  *    @brief performed double-gyroaverage over Maxwellian background
  *
  *    \f[
  *        \int_0^\infty <<\phi>> e^{-\mu} d\mu
  *    \f]
  *
  *    @params In  Input array
  *    @params Out Output array
  *    @params m   Integration m
  *    @params s   Species index
  *
  **/
  virtual void doubleGyroExp(const CComplex In [Nq][NzLD][Nky][NxLD], 
                                   CComplex Out[Nq][NzLD][Nky][NxLD], 
                              const int m, const int s) = 0;
 

  /**
  *  @brief four-dimensional array hold the source terms in drift-coordinates.
  *
  *  @note  Q stand for the German Quellterme (translated : source terms).
  *
  *
  **/
  CComplex *Q,  ///< Hold source terms (Quellterme)
           *Qm; ///< Array to held intermediate Q results

   
  /**
  *  @brief four-dimensional array hold the field terms in drift-coordinates.
  *
  *  Holds the (guiding center) field values. 
  *
  *  with boundaries [NzLlD:NzLD][0:Nky][NxLlD:NxLD]
  *
  *
  **/
  CComplex *Field0;

   
  /**
  *  @brief gyro-averaged field quantities as \f$ \left( <\phi>, <A_{1\parallel}>, <B_{1\parallel}> \right)\f$
  *
  *  Hold the gyro-averaged (gyro-center) field quantities 
  *
  **/
  CComplex *Field;

   
  nct::allocate ArrayField0, ///< Array allocator for Field0 variable
                ArrayField ; ///< Array allocator for Field variable
   
  /**
  * 
  *  @brief \f$ \tfrac{1}{2} \hat{e} \hat{B} \f$ normalization constants 
  *
  *  Brackets should be 1/2 but due to numerical error we should include calculate ourselves, see Dannert[2]
  *
  *  \f[ Y = \frac{1}{\sqrt{\pi}} \int_{-\infty}^\infty v_\parallel^2 e^{-v_\parallel^2} dv \equiv \frac{1}{2} \f] 
  *
  *   but due to discretization errors, this is not exact, thus needs to be calculated numerically.
  *
  *   \f[ Yeb = Y \hat{\epsilon}_{geo} \beta  \f]
  *
  *   @cite Dannert_2006:PhDThesis, where which chapter ... ?
  **/
  double Yeb;
  

  
  /**
  *  @brief solve the field equations
  *
  *
  **/
  void solve(const CComplex *f0, CComplex  *f, Timing timing=0);


  /**
  *  @brief set which equation to solve if one field is assumed to be fixed
  *
  *
  *
  **/
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
  **/ 
  virtual void writeData(const Timing &timing, const double dt);


 protected:

  /**
   *
  *  @brief please document me
  *
  *
  **/
  FileAttr *FA_fields    , ///< FileAttribute for fields (phi,Ap,Bp) 
           *FA_fieldsTime; ///< FileAttribute for time 


  /**
  *  @brief please document me
  *
  **/
  Timing dataOutputFields;


  /**
  *  @brief please document me
  *
  **/
  void initData(Setup *setup, FileIO *fileIO);
  

  /**
  *
  * 
  *  @brief please document me
  *
  *
  *
  **/
  void closeData();


  /**
  * 
  *  @brief please document me
  *
  *
  **/
  virtual void printOn(std::ostream &output) const ;


};


#endif // __FIELDS_H
