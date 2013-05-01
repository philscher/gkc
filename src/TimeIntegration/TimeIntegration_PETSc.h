/*
 * =====================================================================================
 *
 *       Filename: TimeIntegration_PETSc.h
 *
 *    Description: Implicit Time Integration Interface using PETSc. 
 *
 *         Author: Paul P. Hilscher (2012), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#ifndef TIMEINTEGRATION__PETSC_H__
#define TIMEINTEGRATION__PETSC_H__

#include "Global.h"
#include "TimeIntegration.h"

#include "petscts.h"  

/**
*   @brief Implicit Time Integration Interface for PETSc.
*
*   Using PETSc TS (time steeper) to integration gyrokinetic
*   equation system. The integration is performed implicitly
*   by finding the inverse of the gyrokinetic operator using
*   an iterative Krylov method. Note that PETSc performs
*   the time steeping using Euler backward, which should be
*   unstable for a growing system ! However, maybe it is
*   worth checking out the integrators which come along with
*   Sundials and which PETSc can interface.
*
*   @todo cleanup and interface with sundials
*   @todo show that it actually works for a small problem ..
*
**/
class TimeIntegration_PETSc : public TimeIntegration 
{
    
  TS ts; ///< PETSc time steeper object
  
  Mat A_F1; ///< Matrix (in matrix-free object)

 public:
  /**
  *    Please Document Me !
  *
  **/
  TimeIntegration_PETSc(Setup *setup, Grid *grid, Parallel *parallel, Vlasov *vlasov, Fields *fields,
                       TestParticles *particles, Eigenvalue *eigenvalue, Benchmark *bench);
  
  /**
  *    Please Document Me !
  **/
 ~TimeIntegration_PETSc();

  /**
  *    Please Document Me !
  *
  **/
  double solveTimeStep(Vlasov *vlasov, Fields *fields, TestParticles *particles, Timing &timing);

};
       
#endif // TIMEINTEGRATION_PETSC

