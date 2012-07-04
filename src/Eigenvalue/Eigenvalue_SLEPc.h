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


  /**
   *    Please Document Me !
   *
   **/
    EPS EigvSolver;

  /**
   *    Please Document Me !
   *
   **/
    cmplxd tolerance;
  /**
   *    Please Document Me !
   *
   **/
    Mat A_F1;
    
  /**
   *    Please Document Me !
   *
   **/
    TableAttr *EVTable;
  /**
   *    Please Document Me !
   *
   **/
	hid_t eigvGroupID;

protected:

public:
  /**
   *    Please Document Me !
   *
   **/
  Eigenvalue_SLEPc(FileIO *fileIO, Setup *setup, Grid *grid, Parallel *_parallel);
  /**
   *    Please Document Me !
   *
   **/
  void solve(Vlasov *vlasov, Fields *fields, Visualization *visual, Control *control);
  /**
   *    Please Document Me !
   *
   **/
  virtual ~Eigenvalue_SLEPc();

  /**
   *    Please Document Me !
   * Get Maximum Absolute Eigenvalue 
   *
   **/
  cmplxd getMaxAbsEigenvalue(Vlasov *vlasov, Fields *fields);
protected :
  /**
   *    Please Document Me !
   *
   **/
   void print2On(ostream &output);

private:
  /**
   *    Please Document Me !
   *
   **/
   void initDataOutput(Setup *setup, FileIO *fileIO);
};


#endif // __Eigenvalue_SLEPc_H
