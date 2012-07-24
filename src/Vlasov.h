/*
 * =====================================================================================
 *
 *       Filename: Vlasov.h
 *
 *    Description: Vlasov Solver Interface
 *
 *         Author: Paul P. Hilscher (2009-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __VLASOV_H
#define __VLASOV_H

#include "Global.h"
#include "Parallel.h"
#include "Setup.h"
#include "Fields.h"
#include "Grid.h"
#include "FFTSolver.h"

#include "Collisions.h"


/**
*
*
*
**/
enum Boundary  { BOUNDARY_DIRTY=0,  BOUNDARY_CLEAN=1};


class Event;
class PETScMatrixVector;
/**
*
* @brief  class for solving the Vlasov equation
*
* We solve the 6 dimensional Vlasov equation (X, Y, Z, Vp), m, s 
*
**/
class Vlasov : public IfaceGKC {

   /**
   *    Please Document Me !
   *
   **/
   friend class Event;

   /**
   *    Please Document Me !
   *
   **/
   friend class PETScMatrixVector;

 protected:
   /**
   *    Please Document Me !
   *
   **/
   double collisionBeta;

 private:
   /**
   *    Please Document Me !
   *
   *    Set to periodic  or zero boundary
   *
   **/
   std::string useBoundary;

        
   /**
   *  @brief Buffer for exchange ghost cells in each direction
   */ 
   Array6z  SendYu, SendXu, SendYl, SendXl, SendVl, SendVu, SendZl, SendZu; 
  
   /**
   *  @brief Buffer for exchange ghost cells in each direction
   */ 
   Array6z  RecvYu, RecvXu, RecvYl, RecvXl, RecvVl, RecvVu, RecvZl, RecvZu;
      
   /**
   *    Please Document Me !
   *
   **/
   bool boundary_isclean;
  
   /**
   *   @brief Update boundary conditions in case non-blocking MPI is used
   *
   **/
   int cleanBoundary(Array6z A);
        
   /**
   *   @brief reference to last used array used for non-blocking
   *         boundaries
   *
   **/
   Array6z f_boundary;


 protected:
   /**
   *    Please Document Me !
   **/
   bool useAntiAliasing;

   FFTSolver *fft;
   Parallel *parallel;
   Grid *grid;
   Setup *setup;
   Geometry *geo;
   
   /**
   *    Please Document Me !
   **/
   Collisions     *collisions;

   /**
   *   Stabilize simulation by adding a small amount of hyper-viscosity.
   *   Especially needed in x-direction to avoid odd-even decoupling.
   *     
   *   Is this similar effect as mentioned in 
   *    JJ. Quirk , A contribution to the great Riemann solver debate, 1994 ?
   *
   *   Values between 1.e-5 - 1.e-3 are OK, otherwise it will have impact on
   *   physical results.
   *
   **/
   double hyper_visc[DIR_ALL];        
  
   /**
   *    @brief name of Vlasov equation to solve 
   *
   **/
   std::string equation_type;

   /**
   *    @brief set if non-linear simulations are performed
   *
   **/
   bool calculate_nonLinear;
    
   /**
   *
   *  Set the Krook operator 
   *  \f[
   *     \frac{\partial g_{1sigma}}{\partial t} = \dots - \nu(x) g_{1\sigma}
   *  \f]
   * Is used to damp oscillations close to the simulation boundary.
   *
   *  Note : 
   *          * Is this the no-slip boundary condition ?
   *          * Violates conservation of particles, energy and momentum and
   *            needs to be fixed by modifing the fields. See Lapillone.
   *
   *
   **/
   double krook_nu;

        
   /**
   *
   *  Interface to solve Vlasov equation of name equation_type
   *
   *
   **/
   virtual int solve(std::string equation_type, Fields *fields, Array6z fs, Array6z fss, double dt, int rk_step, int user_boundary_type=BOUNDARY_CLEAN) = 0;
 public:
  
   /**
   *    Please Document Me !
   *
   **/
   Array1d  Xi_max;
  
   /**
   *   @brief 6-dimensional phase-space function
   *
   *   @note should this be allocated in TimeIntegation or set to 7dim ?
   *
   **/
   Array6z f0, f, fs, fss, ft, f1;
   
   /**
   *    Please Document Me !
   *
   **/
   Array4z G, Xi;

   /**
   *    Please Document Me !
   *
   **/
   Vlasov(Grid *grid, Parallel *parallel, Setup *setup, FileIO *fileIO, Geometry *geo, FFTSolver *fft);

   /**
   *
   **/
   virtual ~Vlasov();

   /**
   *
   *  Calls Vlasov::solve(equation type ...)
   *  Handles boundary conditions
   *
   **/
   int solve(Fields *fields, Array6z fs, Array6z fss, double dt, int rk_step, int user_boundary_type=BOUNDARY_CLEAN);

   /**
   *    Please Document Me !
   *
   **/
   int setBoundary(Array6z  A, int boundary_type=BOUNDARY_CLEAN);


   /**
   *    Please Document Me !
   *
   **/
   std::string getEquationType() const { return equation_type;};

   
   /**
   *    for global f simulation we update f
   *
   *
   *   ToDo : Skeleton, update new Naxwellian background
   *   this is a dummy functions only
   */ 
   int updateMaxwellian() { return GKC_SUCCESS;};

   /**
   *   Calculates the timestep according to a defined CFL number. For the Vlasov equation
   *   several terms needs to be calculated for the highest possible timestep. These are
   *
   **/
   double getMaxTimeStep(int dir, const double maxCFL);


   /**
   *
   *
   *   Updated the CFL (Courant-Friedrich-Levy number). For explicit time-stepening
   *   the CFL value has to be always < 1 to ensure stability of the system
   *   (practically < 0.4).
   *   
   *   Note : Stability is still not guranteed. As the system is unstable. Thus the
   *          time-steppening scheme needs to allows imaginary values e.g.
   *          (RK-3, RK-4, Heun method).
   *
   *   Calculated using ....
   *  @note get linking error if defined inline. Check Performance !
   **/
   void updateCFL(const cmplxd dphi_dx, const cmplxd dphi_dy, const cmplxd dphi_dz);
 
 protected :

   /**
   *    Please Document Me !
   *
   **/
   void printOn(ostream &output) const;

   /**
   *    Please Document Me !
   *
   **/
   virtual void initDataOutput(FileIO *fileIO);

   /**
   *    Please Document Me !
   *
   **/
   virtual void writeData(Timing *timing) {};

   /**
   *    Please Document Me !
   *
   **/
   virtual void closeData() {};


};


#endif //  __VLASOV_H

