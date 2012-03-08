/*
 * =====================================================================================
 *
 *       Filename: FieldsHypre.cpp
 *
 *    Description: Interface to Hypre Matrix Solver routine
 *
 *         Author: Paul P. Hilscher (2011), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#include "FieldsPETSc.h"


FieldsPETSc::FieldsPETSc(Setup *setup, Grid *grid, Parallel *parallel, FileIO *fileIO, Geometry<HELIOS_GEOMETRY> *geo, FFTSolver *fft) : FieldsFFT(setup, grid, parallel,fileIO, geo, fft)
{

  std::string linearSolver = "PETSc";
  std::string solverType   = "LU";
  
  PetscErrorCode ierr;

  PetscInitialize(&setup->argc, &setup->argv, (char *) 0,  help);
  
  

  int local_size  = NxLD * NkyLD ;
  int global_size = grid->NxGD * NkyLD ;

  // Create Matrix
  
  ierr = MatCreate(PETSC_COMM_WORLD,&A);
  ierr = MatSetSizes(A,local_size, local_size, global_size, global_size);
  ierr = MatSetFromOptions(A);
  
  ierr = MatMPIAIJSetPreallocation(A,5,PETSC_NULL,5,PETSC_NULL);
  ierr = MatSeqAIJSetPreallocation(A,5,PETSC_NULL);

  setupMatrix();

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  ierr = MatAssemblyEnd  (A,MAT_FINAL_ASSEMBLY);

  /* A is symmetric. Set symmetric flag to enable ICC/Cholesky preconditioner */
  ierr = MatSetOption(A,MAT_SYMMETRIC,PETSC_TRUE);

    VecCreateMPI(MPI_COMM_WORLD, local_size, global_size, &Vec_x);
    VecCreateMPI(MPI_COMM_WORLD, local_size, global_size, &Vec_b);
    
    
    VecAssemblyBegin(Vec_b); VecAssemblyEnd(Vec_b);
    VecAssemblyBegin(Vec_x); VecAssemblyEnd(Vec_x);

    //Create linear solver context
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);
    ierr = KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);

    //  '-ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps -mat_mumps_icntl_7 2' 
    //  '-ksp_type preonly -pc_type ilu -pc_factor_mat_solver_package superlu -mat_superlu_ilu_droptol 1.e-8' 
    //  '-ksp_type preonly -pc_type lu/ilu -pc_factor_mat_solver_package petsc' 
  
    ierr = KSPSetType(ksp,KSPPREONLY);
    ierr = KSPGetPC(ksp,&pc);
    ierr = KSPSetFromOptions(ksp);
  
  
    if      (linearSolver == "SuperLU") ierr = PCFactorSetMatSolverPackage(pc,MATSOLVERSUPERLU);
    else if (linearSolver == "PETSc"  ) ierr = PCFactorSetMatSolverPackage(pc,MATSOLVERPETSC);
    else if (linearSolver == "MUMPS"  ) ierr = PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS);
    else check(-1, DMESG("No such solver"));
  
    if      (solverType == "LU"      )  ierr = PCSetType(pc,PCLU);
    else if (solverType == "ILU"     )  ierr = PCSetType(pc,PCILU); 
    else if (solverType == "Cholesky")  ierr = PCSetType(pc,PCCHOLESKY); 
    else check(-1, DMESG("No such solver"));

    ierr = PCFactorSetMatSolverPackage(pc,MATSOLVERSUPERLU);
    ierr = PCFactorSetUpMatSolverPackage(pc); /* call MatGetFactor() to create F */
    ierr = PCFactorGetMatrix(pc,&F);
    //ierr = MatSuperluSetILUDropTol(F,1.e-8);
  

    
}

FieldsPETSc::~FieldsPETSc()
{
  PetscErrorCode ierr;
  
  ierr = KSPDestroy(&ksp);
  ierr = VecDestroy(&Vec_x);
  ierr = VecDestroy(&Vec_b); 
  ierr = MatDestroy(&A);

  ierr = PetscFinalize();
}




Array3z FieldsPETSc::solvePoissonEquation(Array3z rho, Timing timing)
{

   cmplxd *x_phi, *b_rho; 
   
   
   for(int z = NzLlD; z <= NzLuD; z++) {
   

        calcRhoStar(Q);
         
        
        VecGetArray    (Vec_b, (PetscScalar **) &b_rho);
        
        for(int x = NxLlD, k = 0; x <= NxLuD; x++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { 
            b_rho[k++] = rho_star(x, y_k, z);
        }}

        VecRestoreArray    (Vec_b, &b_rho);


        // Solve Equation System
        PetscErrorCode ierr = KSPSolve(ksp,Vec_b,Vec_x);
   
        VecGetArrayRead(Vec_x, (const PetscScalar **) &x_phi);

        // back to real imaginary
        for(int x = NxLlD, k = 0; x <= NxLuD; x++)  { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) {
            Field0(x, y_k, z, Field::phi) = x_phi[k++];
        } }
    
        VecRestoreArrayRead(Vec_x, &x_phi);
  
   }



}
   
   
Array4z FieldsPETSc::solveFieldEquations(Array4z Q, Timing timing) 
{
   
  
   // for 3 fields phi and B_perp are coupled and need to be solved differently
   if( solveEq & Field::phi) solvePoissonEquation  (Q(RxLD, RkyLD, RzLD, Field::phi), timing);
   if( solveEq & Field::Ap ) check(-1, DMESG("Not Implemented"));
   if( solveEq & Field::Bpp) check(-1, DMESG("Not Impleemnted"));
    
   return Field0(RxLD, RkyLD, RzLD, RFields);

}



void FieldsPETSc::calcRhoStar(Array4z Q) 
{
    
    //if(plasma->species(0).doGyro)  phi_yz = calcFluxSurfAvrg(fft->kXOut);
    
    fft->rXIn(RxLD, RkyLD, RzLD, RFields) = Q(RxLD, RkyLD, RzLD, RFields);
    fft->solve(FFT_X, FFT_FORWARD, NkyLD * NzLD * plasma->nfields);
    const double rho_t2 = plasma->species(1).T0  * plasma->species(1).m / pow2(plasma->species(1).q * plasma->B0);
    
    #pragma omp parallel for
    for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int z=NzLlD; z<=NzLuD;z++) { for(int x_k=fft->K1xLlD; x_k<= fft->K1xLuD;x_k++) {

          if((x_k == 0) && (y_k == 0)) { fft->kXIn(x_k,y_k,z, Field::phi) = (cmplxd) 0.e0 ; continue; }
          
          const double b = rho_t2 * fft->k2_p(x_k,y_k,z);
          
          //const cmplxd rhs  = (fft->kXOut(x_k,y_k,z, Q::rho) + adiab * phi_yz(x_k))/fft->Norm_X; 
          const cmplxd rhs  = fft->kXOut(x_k,y_k,z, Q::rho)/fft->Norm_X; 
          
          fft->kXIn(x_k,y_k,z, Field::phi) = (1. + b) * rhs; 
         
    } } }
   
    fft->solve(FFT_X, FFT_BACKWARD, NkyLD * NzLD * plasma->nfields);
   
    rho_star(RxLD, RkyLD, RzLD) = fft->rXOut(RxLD, RkyLD, RzLD, Field::phi);

    return;
}

/* 
void FieldsPETSc::calcRhoStar(Array4z Q) 
{


      for(int z = NzLlD; z <= NzLuD; z++) {
      for(int x = NxLlD; x <= NxLuD; x++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) {
                
                const cmplxd d2phi_dx = (16. * (Q(x+1,y_k,z,Field::phi) + Q(x-1,y_k,z,Field::phi)) - (Q(x+2,y_k,z,Field::phi) + Q(x-2,y_k,z,Field::phi)) - 30. * Q(x,y_k,z,Field::phi)) / (12. * pow2(dx));
                const cmplxd d2phi_dy = pow2(fft->ky(y_k)) * Q(x,y_k,z,Field::phi);

                rho_star(x,y_k,z) = Q(x,y_k,z,Field::phi) + (d2phi_dx + d2phi_dy);
      }}

    }

}
 * */
  
void FieldsPETSc::setupMatrix()
{
  
  int idx_l, idx_u;

  PetscErrorCode ierr = MatGetOwnershipRange(A,&idx_l,&idx_u);
 
  
  for (int i=idx_l; i<idx_u; i++) {
           
                int y_k = 1;
    
                const double h = dx*dx;
               
                const double ky2 = pow2(fft->ky(y_k));
                
                const double val_diag  = 1. + (2.+plasma->debye2) * (2./h + ky2) ;
                const double val_right =  - (2.+plasma->debye2)/h ;

                // Set D [ S_ 0 D S]
                ierr = MatSetValues(A,1,&i,1,&i,(const PetscScalar *) &val_diag,INSERT_VALUES);

                int J = i + 1;
                ierr = MatSetValues(A,1,&i,1,&J,(const PetscScalar *) &val_right,INSERT_VALUES);
  
  
  }

    
}



