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
#include "Parallel/Parallel.h"
#include "Setup.h"
#include "Fields/Fields.h"
#include "Grid.h"
#include "FFTSolver/FFTSolver.h"
#include "Analysis/Benchmark.h"

#include "Collisions/Collisions.h"




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
  enum class Boundary  { SEND=1,  RECV=2, SENDRECV=3};


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
   *   @brief Update boundary conditions in case non-blocking MPI is used
   *
   **/
   void cleanBoundary(Array6C A);
        

 protected:
   FFTSolver *fft;
   Parallel *parallel;
   Grid *grid;
   Setup *setup;
   Geometry *geo;
   Benchmark *bench;  
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
   virtual int solve(std::string equation_type, Fields *fields, CComplex *fs, CComplex *fss, double dt, int rk_step, const double rk[3]) = 0;
 public:
  
  /**
   *   @brief Pre-calculate often used factors
   *
   *   Pre-calculate often used factors which e.g. include divisions,
   *   as division requires much more cycles than multiplications.
   *         
   *   KW is Kehrwert (German for Multiplicative Inverse)
   *
   **/
   const double _kw_12_dx_dx,    ///< \f$ \frac{1}{12 dx^2} \f$
                _kw_12_dv   ,    ///< \f$ \frac{1}{16 dv  } \f$
                _kw_12_dx   ,    ///< \f$ \frac{1}{12 dx  } \f$
                _kw_12_dv_dv,    ///< \f$ \frac{1}{12 dv^2} \f$
                _kw_16_dx4  ;    ///< \f$ \frac{1}{16 dx^4} \f$



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
   nct::allocate ArrayPhase;
   CComplex *f0, *f, *fs, *fss, *ft, *f1;
   
   /**
   *    Please Document Me !
   *
   **/
   CComplex *G, *Xi;

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
   int solve(Fields *fields, CComplex *fs, CComplex *fss, double dt, int rk_step, const double rk[3],  bool useNonBlockingBoundary=true);

   /**
   *    Please Document Me !
   *
   **/
   void setBoundary(CComplex *f);


   void setBoundary(CComplex *f, const Boundary boundary_type);

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

