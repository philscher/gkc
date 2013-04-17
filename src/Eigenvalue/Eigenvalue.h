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


#ifndef __GKC_EIGENVALUE_H
#define __GKC_EIGENVALUE_H

#include<array>

#include "Global.h"
#include "Parallel/Parallel.h"
#include "Setup.h"
#include "Control.h"
#include "Vlasov/Vlasov.h"
#include "Fields/Fields.h"

#include "Visualization/Visualization.h"

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
class Eigenvalue : public IfaceGKC {

 protected:

  Parallel *parallel;
  Grid     *grid;
    
  double growth_min;  ///< Minimum growth rate to save eigenvector

  /**
  *  @brief structure to store the eigenvalue
  *
  **/
  struct EigenValue {
    
    Complex EigenValue  ; ///< the eigenvalue
    double AbsoluteError; ///< the absolute error

  };

  bool includeZF; ///< true if Zonal Flow is included

 public:

  /**
  *    Please Document Me !
  *
  **/
  Eigenvalue(FileIO *fileIO, Setup *setup, Grid *_grid, Parallel *_parallel) : parallel(_parallel), grid(_grid) 
  {
    includeZF  = setup->get("Eigenvalue.IncludeZF"    , 0);
    growth_min = setup->get("Eigenvalue.MinGrowthrate", 0.);
  };

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
  *         the maximum stable linear time step
  *
  **/
  virtual Complex getMaxAbsEigenvalue(Vlasov *vlasov, Fields *fields) = 0;

 protected :
  
  /**
  *    Please Document Me !
  *
  **/
  void InitDataOutput(Setup *setup, FileIO *fileIO);

};


#endif // __GKC_EIGENVALUE_H
