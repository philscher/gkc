/*
 * =====================================================================================
 *
 *       Filename: TimeIntegration_PETSc.cpp
 *
 *    Description: Implicit Time Integration Interface for gkc using PETSc. 
 *
 *         Author: Paul P. Hilscher (2012), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include "TimeIntegration_PETSc.h"
#include "Special/PETScMatrixVector.h"

static char help[] = "Help for PETSc Interface not available, please look up gkc & PETSc manual.";

extern PetscErrorCode MatrixVectorProduct(Mat A, Vec x, Vec y);

TimeIntegration_PETSc::TimeIntegration_PETSc(Setup *setup, Grid *grid, Parallel *parallel, Vlasov *vlasov, Fields *fields, Eigenvalue *eigenvalue, Benchmark *bench)
     : TimeIntegration(setup, grid, parallel, vlasov, fields,eigenvalue, bench)
{
      PetscInitialize(&setup->argc, &setup->argv, (char *) 0,  help);
      
      // create Matrix Operations
      
      int subDiv=0;

      // create Matrix, a shell for Matrix-free methods
      MatCreateShell(parallel->Comm[DIR_ALL], grid->getLocalSize(), grid->getLocalSize(), grid->getGlobalSize(), grid->getGlobalSize(), &subDiv, &A_F1);
      MatSetFromOptions(A_F1);
      MatShellSetOperation(A_F1, MATOP_MULT, (void(*)()) PETScMatrixVector::MatrixVectorProduct);
      
      // Initialize implicit solver
      TSCreate(parallel->Comm[DIR_ALL], &ts);
      TSSetProblemType(ts, TS_LINEAR);
      TSSetRHSFunction(ts, PETSC_NULL,TSComputeRHSFunctionLinear, PETSC_NULL);
      TSSetRHSJacobian(ts, A_F1, A_F1, TSComputeRHSJacobianConstant, PETSC_NULL);
      TSSetType(ts, TSBEULER);

      TSSetFromOptions(ts); 

      // Setup inital vector
      Vec Vec_init;
      Complex *init_x = PETScMatrixVector::getCreateVector(grid, Vec_init);
    
      for(int x = NxLlD, n = 0; x <= NxLuD; x++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { for(int z = NzLlD; z <= NzLuD; z++) {
      for(int v = NvLlD       ; v <= NvLuD; v++) { for(int m   = NmLlD ; m   <= NmLuD ; m++  ) { for(int s = NsLlD; s <= NsLuD; s++) {
                init_x[n++] = vlasov->f(x,y_k,z,v,m,s);

      }}} }}}

      VecRestoreArray(Vec_init, (PetscScalar **) &init_x);
      TSSetSolution(ts, Vec_init);
      
      VecDestroy(&Vec_init);


}

TimeIntegration_PETSc::~TimeIntegration_PETSc()
{
        TSDestroy(&ts);

}



double TimeIntegration_PETSc::solveTimeStep(Vlasov *vlasov, Fields *fields, TestParticles *particles, Timing &timing)
{
            if(timeIntegrationScheme == "Implicit_RK4") {
                    PETScMatrixVector pMV(vlasov, fields);
                    TSSetTimeStep(ts, dt);
                    TSStep(ts); 
                    // Get Solution
    Vec    Vec_F1;   
    Complex  *x_F1; 
    

                    TSGetSolution(ts, &Vec_F1);
                    VecGetArray(Vec_F1, (PetscScalar **) &x_F1);

            // copy whole phase space function (waste but starting point) (important due to bounday conditions
           // we can built wrapper around this and directly pass it
   for(int x = NxLlD, n = 0; x <= NxLuD; x++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { for(int z = NzLlD; z <= NzLuD; z++) {
   for(int v = NvLlD       ; v <= NvLuD; v++) { for(int m   = NmLlD ; m   <= NmLuD ; m++  ) { for(int s = NsLlD; s <= NsLuD; s++) {

                vlasov->fs(x,y_k,z,v,m,s) = x_F1[n];
                vlasov->f(x,y_k,z,v,m,s) = x_F1[n++];

           }}} }}}
           
           fields->solve(vlasov->f0,  vlasov->f); 
           VecRestoreArray    (Vec_F1, (PetscScalar **) &x_F1);
            
           timing.time += dt;
            timing.step++;
            writeTimeStep(timing, maxTiming, dt);
        
            }
            else TimeIntegration::solveTimeStep(vlasov, fields, particles, timing);

      return dt;
}



