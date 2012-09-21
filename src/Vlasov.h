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
#include "Analysis/Benchmark.h"

#include "Collisions.h"




class Event;
class PETScMatrixVector;
/**
*
* @brief  class for solving the Vlasov equation
*
* We solve the 6 dimensional Vlasov equation 
*
* \f[
*     f_{1\sigma}(x, k_y, z, v_\parallel, \mu, \sigma)
* \f]
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

   /**
   *    Please Document Me !
   *
   **/
   friend class Benchmark;

 protected:
  /**
  *
  **/
  enum Boundary  { BOUNDARY_SEND=1,  BOUNDARY_RECV=2, BOUNDARY_SENDRECV=3};


   /**
   *    Please Document Me !
   *
   **/
   double collisionBeta;

 private:
   /**
   *  @brief Buffer for exchange ghost cells in each direction
   */ 
   Array6C  SendYu, SendXu, SendYl, SendXl, SendVl, SendVu, SendZl, SendZu; 
  
   /**
   *  @brief Buffer for exchange ghost cells in each direction
   */ 
   Array6C  RecvYu, RecvXu, RecvYl, RecvXl, RecvVl, RecvVu, RecvZl, RecvZu;
      
   /**
   *    Please Document Me !
   *
   **/
   bool boundary_isclean;
  
   /**
   *   @brief Update boundary conditions in case non-blocking MPI is used
   *
   **/
   void cleanBoundary(Array6C A);
        
   /**
   *   @brief reference to last used array used for non-blocking
   *         boundaries
   *
   **/
   Array6C f_boundary;


 protected:
   FFTSolver *fft;
   Parallel *parallel;
   Grid *grid;
   Setup *setup;
   Geometry *geo;
   Benchmark *bench;  

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
   bool nonLinear;
    

        
   /**
   *
   *  Interface to solve Vlasov equation of name equation_type
   *
   *
   **/
   virtual int solve(std::string equation_type, Fields *fields, Array6C fs, Array6C fss, double dt, int rk_step, const double rk[3]) = 0;
 public:
  
   /**
   *    Please Document Me !
   *
   **/
   double  Xi_max[3];
  
   /**
   *   @brief 6-dimensional phase-space function
   *
   *   @note should this be allocated in TimeIntegation or set to 7dim ?
   *
   **/
   Array6C f0, f, fs, fss, ft, f1;
   
   /**
   *    Please Document Me !
   *
   **/
   Array4C G, Xi;

   /**
   *    Please Document Me !
   *
   **/
   Vlasov(Grid *grid, Parallel *parallel, Setup *setup, FileIO *fileIO, Geometry *geo, FFTSolver *fft, Benchmark *bench);

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
   int solve(Fields *fields, Array6C fs, Array6C fss, double dt, int rk_step, const double rk[3],  bool useNonBlockingBoundary=true);

   /**
   *    Please Document Me !
   *
   **/
   void setBoundary(Array6C);


   void setBoundary(Array6C, int boundary_type);

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
   *    @brief Holds the non-linear terms from the ExB
   *
   **/
   CComplex *nonLinearTerms; 

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
   *
   *   This needs only to be caluclated in the non-linear terms
   *
   *  @note get linking error if defined inline. Check Performance !
   *
   *  @depracated This is not done direclty in the non-linear term
   **/
   /*
   void inline updateCFL(const Complex dphi_dx, const Complex dphi_dy, const Complex dphi_dz)
{
  Xi_max[DIR_X] = max(Xi_max[DIR_X], abs(dphi_dx));
  Xi_max[DIR_Y] = max(Xi_max[DIR_Y], abs(dphi_dy));
  Xi_max[DIR_Z] = max(Xi_max[DIR_Z], abs(dphi_dz));
};
 */
 
 protected :

   /**
   *    Please Document Me !
   *
   **/
   void printOn(std::ostream &output) const;

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

