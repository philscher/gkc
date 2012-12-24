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

#include "Vlasov/Vlasov.h"
#include "Fields/Fields.h"
#include "Grid.h"

#include "petsc.h"

PetscErrorCode MatrixVectorProduct(Mat A, Vec x, Vec y);


int petc_signal_handler(int sig, void *ctx);

/**
*    @brief PETSc interface for Matrix-Vector multiplication
*
*
**/
class PETScMatrixVector
{

 public:
    
  /** 
  *    Note : Needs to be called to initlized global variables
  **/
  PETScMatrixVector(Vlasov *vlasov, Fields *fields, bool includeZF=false);
   
  /** 
  *  
  *  @brief Matrix-Vector multiplication with the Vlasov equation
  *
  *  Interface for matrix-free matrix-vector multiplication.
  *  With 
  *  
  *  \f[
  *        f_{1\sigma} = L_{gy} f_{1\sigma}'
  *  \f]
  *       
  *
  **/
  static PetscErrorCode MatrixVectorProduct(Mat A, Vec Vec_x, Vec Vec_y);

  /**
  *
  * 
  **/
  static CComplex* getCreateVector(Grid *grid, Vec &Vec_x);

};


#endif // __GKC_PETSC_MATRIX_VECTOR__H_

