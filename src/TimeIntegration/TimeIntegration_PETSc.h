/*
 * =====================================================================================
 *
 *       Filename: TimeIntegration_PETSc.h
 *
 *    Description: Implicit Time Integration Interface for gkc using PETSc. 
 *
 *         Author: Paul P. Hilscher (2012), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#include "config.h"
#include "Global.h"

#ifndef TIMEINTEGRATION__PETSC_H__
#define TIMEINTEGRATION__PETSC_H__


#include "TimeIntegration.h"
#include "petscts.h"  

/**
 *
 *    Implicit Time Integration Interface for PETSc.
 *    Uses SLEPc Interface.
 *
 *
 *
 *
 **/
class TimeIntegration_PETSc : public TimeIntegration {
    
    TS ts;
    Mat A_F1;

    public:
        TimeIntegration_PETSc(Setup *setup, Grid *grid, Parallel *parallel, Vlasov *vlasov, Fields *fields, Eigenvalue *eigenvalue);
        double solveTimeStep(Vlasov *vlasov, Fields *fields, TestParticles *particles, Timing &timing);

        ~TimeIntegration_PETSc();
  private:

};
       



#endif // TIMEINTEGRATION_PETSC
