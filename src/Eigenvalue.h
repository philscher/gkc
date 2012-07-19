/*
 * =====================================================================================
 *
 *       Filename: Eigenvalue.h
 *
 *    Description: Interface for eigenvalue calulcations
 *
 *         Author: Paul P. Hilscher (2011-)
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#ifndef __Eigenvalue_H
#define __Eigenvalue_H

#include "Global.h"
#include "Parallel.h"
#include "Setup.h"
#include "Control.h"
#include "Vlasov.h"
#include "Fields.h"

#include "Visualization.h"

/**
*    @brief  Interface for eigenvalue calculations
*
*    Can be used to calculate the eigenvalues of the
*    current setup, use is rather limited due to the
*    problem size - or to calculate the maximum absolute
*    eigenvalue in order to estimate the time step.
*
*    @cite Merz2012:MultiDimensionalGKParameterStudy
*    @cite Roman2010:FastEigenvalueMassivelyParallel
*
**/
class Eigenvalue : public IfaceHelios {

  protected:
   /**
   *  @brief structu to store the eigenvalue
   *
   **/
   struct EigenValue {
    
     cmplxd EigenValue   ; ///< the eigenvalue
     double AbsoluteError; ///< the absolute error

   };


   Parallel *parallel;
   Grid     *grid;

  public:
   /**
   *    Please Document Me !
   *
   **/
   Eigenvalue(FileIO *fileIO, Setup *setup, Grid *_grid, Parallel *_parallel) : parallel(_parallel), grid(_grid) {};

   /**
   *    Please Document Me !
   *
   **/
   virtual void solve(Vlasov *vlasov, Fields *fields, Visualization *visual, Control *control) = 0;
  
   /**
   *    Please Document Me !
   *
   **/
   virtual ~Eigenvalue() {};

  
   /**
   *  @brief Get Maximum Absolute Eigenvalue for determining
   *         the maximum stable linear timestep
   *
   **/
   virtual cmplxd getMaxAbsEigenvalue(Vlasov *vlasov, Fields *fields) = 0;

  protected :
  
   /**
   *    Please Document Me !
   *
   **/
   void printOn(ostream &output) const { };

   /**
   *    Please Document Me !
   *
   **/
   void InitDataOutput(Setup *setup, FileIO *fileIO);

};


#endif // __Eigenvalue_H
