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
  enum class Boundary : int { SEND=1,  RECV=2, SENDRECV=3};


 private:
   /**
   *  @brief Buffer for exchange ghost cells in each direction
   */ 
   CComplex *SendXl, *SendXu, *SendZl, *SendZu, *SendVl, *SendVu;
   CComplex *RecvXl, *RecvXu, *RecvZl, *RecvZu, *RecvVl, *RecvVu;
   nct::allocate ArrayBoundX, ArrayBoundZ, ArrayBoundV;
    

 protected:
   FFTSolver *fft;
   Parallel *parallel;
   Grid *grid;
   Setup *setup;
   Geometry *geo;
   Benchmark *bench;  
   Collisions     *coll;


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
   bool doNonLinear;
    
        
   /**
   *
   *  Interface to solve Vlasov equation of name equation_type
   *
   *  set fs to const f_in 
   **/
   virtual void solve(std::string equation_type, Fields *fields, CComplex *fs, CComplex *fss, double dt, int rk_step, const double rk[3]) = 0;
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
                _kw_12_dz   ,    ///< \f$ \frac{1}{12 dz  } \f$
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
   CComplex *f0,         ///< Maxwellian 
            *f ,         ///< Perturbed Phase Space Function
            *fs,         ///<
            *fss,
            *ft, 
            *f1,         ///< f1
            *Coll;       ///< Collisional corrections
   
   /**
   *    Please Document Me !
   *
   **/
   CComplex *G, *Xi;

   nct::allocate ArrayG,  ///< Allocation class for G
                 ArrayXi, ///< Allocation class for Xi
                 ArrayNL; ///< Allocation class for non-linear term
   /**
   *    Please Document Me !
   *
   **/
   Vlasov(Grid *grid, Parallel *parallel, Setup *setup, FileIO *fileIO, Geometry *geo, FFTSolver *fft, Benchmark *bench, Collisions *coll);

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
   void solve(Fields *fields, CComplex *fs, CComplex *fss, double dt, int rk_step, const double rk[3],  bool useNonBlockingBoundary=true);

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
   //int updateMaxwellian() { return GKC_SUCCESS;};

   /**
   *   @brief get maximum non-linear time step 
   *
   *   Calculates the timestep according to a defined CFL number. For the Vlasov equation
   *   several terms needs to be calculated for the highest possible timestep. These are
   *
   *   @param maxCFL maximum CFL value for non-linear terms
   *
   *   @return maximum time step value
   *
   **/
   double getMaxNLTimeStep(const double maxCFL);

   
   
   /**
   *    @brief Holds the non-linear terms from the ExB
   *
   **/
   CComplex *nonLinearTerm; 
   
   /**
   *    Please Document Me !
   *
   **/
   virtual void writeData(const Timing &timing, const double dt);


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
   virtual void initData(Setup *setup, FileIO *fileIO);
//   virtual void initData(hid_t groupID, FileIO *fileIO);

   void loadData(FileIO *fileIO);

   Timing dataOutputF1;

   /**
   *    Please Document Me !
   *
   **/
   virtual void closeData();

 private:
   FileAttr *FA_f1, *FA_f0, *FA_psfTime; 

};


#endif //  __VLASOV_H

