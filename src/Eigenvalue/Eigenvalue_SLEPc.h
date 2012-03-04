/*
 * =====================================================================================
 *
 *       Filename: Eigenvalue_SLEPc.h
 *
 *    Description: Eigenvalue calculation using SLEPc
 *
 *         Author: Paul P. Hilscher (2011), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include "Global.h"
#include "config.h"


#ifndef __Eigenvalue_SLEPc_H
#define __Eigenvalue_SLEPc_H


#include "Eigenvalue.h"

#include <slepceps.h>

/**
 *    Calculates the eigenvlaues of the system using SLEPc
 *
 *
 *
 *
 * */
class Eigenvalue_SLEPc : public Eigenvalue {


    EPS EigvSolver;

    cmplxd tolerance;
    Mat A_F1;
    
    TableAttr *EigVal_TableAttr;

protected:

public:
  Eigenvalue_SLEPc(FileIO *fileIO, Setup *setup, Grid *grid, Parallel *_parallel);
  void solve(Vlasov *vlasov, Fields *fields, Visualization *visual, Control *control);
  virtual ~Eigenvalue_SLEPc();

  /** Get Maximum Absolute Eigenvalue */
  cmplxd getMaxAbsEigenvalue(Vlasov *vlasov, Fields *fields);
protected :
   void print2On(ostream &output);

private:
   void initDataOutput(Setup *setup, FileIO *fileIO);
   //void write();
};


#endif // __Eigenvalue_SLEPc_H
