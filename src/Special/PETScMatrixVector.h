/*
 * =====================================================================================
 *
 *       Filename: PETScMatrixVector.cpp
 *
 *    Description: Calculates y = F(x), with x some input vector. 
 *                 Required for PETSc/SLEPc Matrix-free methods. 
 *
 *         Author: Paul P. Hilscher (2012), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#include "config.h"
#include "Global.h"


#ifndef __GKC_PETSC_MATRIX_VECTOR__H_
#define __GKC_PETSC_MATRIX_VECTOR__H_

#include "Vlasov.h"
#include "Fields.h"
#include "Grid.h"

#include "petsc.h"

PetscErrorCode MatrixVectorProduct(Mat A, Vec x, Vec y);


int petc_signal_handler(int sig, void *ctx);




class PETScMatrixVector
{
  public:
        PETScMatrixVector(Vlasov *vlasov, Fields *fields);
        static PetscErrorCode MatrixVectorProduct(Mat A, Vec Vec_x, Vec Vec_y) ;
        static cmplxd* getCreateVector(Grid *grid, Vec &Vec_x);
};


#endif // __GKC_PETSC_MATRIX_VECTOR__H_

