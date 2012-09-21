/*
 * =====================================================================================
 *
 *       Filename: MatrixSolver.h
 *
 *    Description: Interface for PETSc matrix solver routines
 *
 *         Author: Paul P. Hilscher (2012), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef PETSC_MATRIX_SOLVER
#define PETSC_MATRIX_SOLVER


#include <petscksp.h>

#include "Parallel/Parallel.h"

/**
*
*  @brief Solves a linear equation system using PETSc (http://www.mcs.anl.gov/petsc/).
*
**/
class MatrixSolver
{
 
   KSP      ksp; 
   PC pc    ;
   Mat     F;

  public:
  
   /**
   *    Please Document Me !
   *
   **/
   MatrixSolver( Parallel *parallel, Mat Matrix, int dir, bool solverTypeDirect=false, std::string linearSolver="PETSc")
   {

        //Create linear solver context
        //
        PetscErrorCode ierr;

        ierr = KSPCreate(parallel->Comm[dir], &ksp);
        ierr = KSPSetOperators(ksp, Matrix, Matrix, SAME_PRECONDITIONER);


        if(solverTypeDirect == true) ierr = KSPSetType(ksp,KSPPREONLY);

        ierr = KSPGetPC(ksp, &pc);
  
        std::string pcType = "LU";

        if(solverTypeDirect == true) {

          if      (pcType == "LU"      )  ierr = PCSetType(pc,PCLU);
          else if (pcType == "ILU"     )  ierr = PCSetType(pc,PCILU); 
          else if (pcType == "Cholesky")  ierr = PCSetType(pc,PCCHOLESKY); 
          else check(-1, DMESG("No such solver"));
 

          //if      (linearSolver == "SuperLU") ierr = PCFactorSetMatSolverPackage(pc,MATSOLVERSUPERLU_DIST);
          // ToDo : Distinguish between SuperLU_disc, Super_LU, super
          if      (linearSolver == "SuperLU") ierr = PCFactorSetMatSolverPackage(pc,MATSOLVERSUPERLU);
          else if (linearSolver == "PETSc"  ) ierr = PCFactorSetMatSolverPackage(pc,MATSOLVERPETSC);
          else if (linearSolver == "MUMPS"  ) ierr = PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS);
          else if (linearSolver == "Spooles") ierr = PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS);
          else check(-1, DMESG("No such solver"));
  
          PCFactorSetUpMatSolverPackage(pc);
          PCFactorSetColumnPivot(pc, 1.);
          //PCFactorReorderForNonzeroDiagonal(pc, 1.e-6);
          ierr = PCFactorGetMatrix(pc, &F);

        } else {
       
           if (pcType == "Jacobi")  ierr = PCSetType(pc,PCJACOBI); 
           else check(-1, DMESG("No such solver"));
     
           KSPSetTolerances(ksp,1.e-9,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
        
        }    
  
        KSPSetFromOptions(ksp);
        // No Need to call this ?!
    //    KSPSetUp(ksp);
  
   };
  
   /**
   *    Please Document Me !
   *
   **/
   bool solve(Vec &Vec_b, Vec &Vec_x) 
   {
        
        PetscErrorCode ierr = KSPSolve(ksp, Vec_b, Vec_x); 
        return (ierr == PETSC_TRUE) ? true : false;
   }

   /**
   *    Please Document Me !
   *
   **/
   ~MatrixSolver() 
   {

    PetscErrorCode ierr = KSPDestroy(&ksp);
 
   }

};



#endif // PETSC_MATRIX_SOLVER
