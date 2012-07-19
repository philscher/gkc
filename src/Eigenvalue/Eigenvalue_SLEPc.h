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
*   @brief Calculates the eigenvalues using SLEPc
*          
*   SLEPc  uses Krylov Matrix-Free methods in order
*   to extract the eigenvalues of the system. 
*  
*   See <a href http://www.grycap.upv.es/slepc/>SLEPc Homepage</a>
*
**/
class Eigenvalue_SLEPc : public Eigenvalue {


    EPS EigvSolver;  ///< The eigensolver context

    cmplxd tolerance; ///< Set Tolerence of the solution
    
    Mat A_F1;         ///< The Matrix Context
    
    TableAttr *EVTable; ///< Eigenvalue Table
	
    hid_t eigvGroupID; ///< HDF-5 reference to Eigenvalue group

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
   void printOn(ostream &output);

  /**
   *    Please Document Me !
   *
   **/
   void initDataOutput(Setup *setup, FileIO *fileIO);
};


#endif // __Eigenvalue_SLEPc_H
