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

class Eigenvalue : public IfaceHelios {

struct EigenValue {
    cmplxd EigenValue;
    cmplxd AbsoluteError;
};


protected:
  Parallel *parallel;
  Grid     *grid;

public:
  Eigenvalue(FileIO *fileIO, Setup *setup, Grid *_grid, Parallel *_parallel) : parallel(_parallel), grid(_grid) {};
  virtual void solve(Vlasov *vlasov, Fields *fields, Visualization *visual, Control *control) = 0;
  virtual ~Eigenvalue() {};

  /* Get Maximum Absolute Eigenvalue */
  virtual cmplxd getMaxAbsEigenvalue(Vlasov *vlasov, Fields *fields) = 0;
protected :
   void printOn(ostream &output) const { };

private:
   void InitDataOutput(Setup *setup, FileIO *fileIO);
};


#endif // __Eigenvalue_H
