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
   *    Please Document Me !
   *
   **/
class Eigenvalue : public IfaceHelios {

  protected:
  /**
   *    Please Document Me !
   *
   **/
struct EigenValue {
  /**
   *    Please Document Me !
   *
   **/
    cmplxd EigenValue;
  /**
   *    Please Document Me !
   *
   **/
    double AbsoluteError;
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
   *    Please Document Me !
   * Get Maximum Absolute Eigenvalue 
   *
   **/
  virtual cmplxd getMaxAbsEigenvalue(Vlasov *vlasov, Fields *fields) = 0;
protected :
  /**
   *    Please Document Me !
   *
   **/
   void printOn(ostream &output) const { };

private:
  /**
   *    Please Document Me !
   *
   **/
   void InitDataOutput(Setup *setup, FileIO *fileIO);
};


#endif // __Eigenvalue_H
