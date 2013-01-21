/*
 * =====================================================================================
 *
 *       Filename: PETScMatrixVector.cpp
 *
 *    Description: Matrix-Scalar multiplication helper for PETSc
 *
 *         Author: Paul P. Hilscher (2012-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include "Fields.h"

#include "PETScMatrixVector.h"

Vlasov *GL_vlasov; /// global variable to use for MatrixVectorProduct
Fields *GL_fields;

bool GL_includeZF;
int  GL_iter;

extern int process_rank;

int petc_signal_handler(int sig, void *ctx) 
{
  // hard exit ( try to improve using control)
  check(-1, DMESG("PETSc signal received"));

  return 0;
}



PETScMatrixVector::PETScMatrixVector(Vlasov *vlasov, Fields *fields, bool includeZF)
{
  // Set global variables to use for ::MatrixVectorProduct
  GL_includeZF = includeZF;
  GL_vlasov    = vlasov;
  GL_fields    = fields;
  GL_iter      = 0;
}


PetscErrorCode PETScMatrixVector::MatrixVectorProduct(Mat A, Vec Vec_x, Vec Vec_y) 

{

  [=] (CComplex  fs [NsLD][NmLD][NzLB][Nky][NxLB][NvLB],  
       CComplex  fss[NsLD][NmLD][NzLB][Nky][NxLB][NvLB])
  {
      
    if(process_rank == 0 ) std::cout << "\r"   << "Iteration  : " << GL_iter++ << std::flush;

    CComplex *x_F1, *y_F1; 
 
    VecGetArrayRead(Vec_x, (const PetscScalar **) &x_F1);
    VecGetArray    (Vec_y, (      PetscScalar **) &y_F1);

    // copy f1 to vector (important due to bounday conditions, thus cannot pass directly)
    int n = 0;
    for(int s = NsLlD; s <= NsLuD; s++) { for(int m = NmLlD; m <= NmLuD; m++) { for(int z = NzLlD; z <= NzLuD; z++) {
    for(int y_k = (GL_includeZF ? 0 : 1); y_k < Nky-1; y_k++) {
    for(int x = NxLlD; x <= NxLuD; x++) { for(int v = NvLlD; v <= NvLuD; v++) { 
      
      fs[s][m][z][y_k][x][v] = x_F1[n++];

    }}} }}}

    GL_fields->solve(GL_vlasov->f0, GL_vlasov->fs); 
   
    // Set zero time integration coefficient so that fss = F_gy(fs) 
    const double rk_0[] = { 0., 0., 0.}; 
    // Calculate the collision operator
    GL_vlasov->solve(GL_fields, GL_vlasov->fs, GL_vlasov->fss, 1., 0, rk_0, false);
   
    // copy whole phase space function to PETSc vector
    //#pragma omp parallel for, collapse private(n) 
    n = 0;
    for(int s = NsLlD; s <= NsLuD; s++) { for(int m = NmLlD; m <= NmLuD; m++) { for(int z = NzLlD; z <= NzLuD; z++) {
    for(int y_k = (GL_includeZF ? 0 : 1); y_k < Nky-1; y_k++) {   // iterate from y_k=0 only if Zonal Flow is included
    for(int x = NxLlD; x <= NxLuD; x++) { for(int v = NvLlD; v <= NvLuD; v++) { 

      y_F1[n++] = fss[s][m][z][y_k][x][v];

    }}} }}}
 
    VecRestoreArrayRead(Vec_x, (const PetscScalar **) &x_F1);
    VecRestoreArray    (Vec_y, (      PetscScalar **) &y_F1);

    } ((A6zz) GL_vlasov->fs, (A6zz) GL_vlasov->fss);
   
  return 0; // return 0 (success) required for PETSc

}

CComplex* PETScMatrixVector::getCreateVector(Grid *grid, Vec &Vec_x, int NkyRed) 
{

  CComplex *xp;

  // Usually, Zonal flow/Nyquiest frequency is negelcted
  int getGlobalSize =    Nx * (Nky - NkyRed) * Nz   * Nv   * Nm   * Ns; 
  int getLocalSize  =  NxLD * (Nky - NkyRed) * NzLD * NvLD * NmLD * NsLD;
  

  VecCreateMPI(MPI_COMM_WORLD, grid->getLocalSize(), grid->getGlobalSize(), &Vec_x);
  VecAssemblyBegin(Vec_x);
  VecAssemblyEnd(Vec_x);

  VecGetArray    (Vec_x, (PetscScalar **) &xp);

  return xp;

}
 


