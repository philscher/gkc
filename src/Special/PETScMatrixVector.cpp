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

Vlasov *GL_vlasov;
Fields *GL_fields;

bool GL_includeZF;
int  GL_iter;

extern int process_rank;

int petc_signal_handler(int sig, void *ctx) {
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
};


PetscErrorCode PETScMatrixVector::MatrixVectorProduct(Mat A, Vec Vec_x, Vec Vec_y) 

{

  [=] (CComplex  fs  [NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB],  
       CComplex  fss [NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB])
  {
      
    if(process_rank == 0 ) std::cout << "\r"   << "Iteration  : " << GL_iter++ << std::flush;

    CComplex *x_F1, *y_F1; 
 
    VecGetArrayRead(Vec_x, (const PetscScalar **) &x_F1);
    VecGetArray    (Vec_y, (      PetscScalar **) &y_F1);

    // copy whole phase space function (waste but starting point) (important due to bounday conditions
    // we can built wrapper around this and directly pass it
    int n = 0;
    for(int s = NsLlD; s <= NsLuD; s++) { for(int m   = NmLlD   ; m   <= NmLuD   ; m++  ) { for(int z = NzLlD; z <= NzLuD; z++) {
    for(int y_k = (GL_includeZF ? 0 : 1); y_k <= NkyLuD-1; y_k++) {   // iterate from NkyLlD if Zonal Flow is included
    for(int x = NxLlD; x <= NxLuD; x++) { for(int v = NvLlD       ; v <= NvLuD; v++) { 
      
           fs[s][m][z][y_k][x][v] = x_F1[n++];

   }}} }}}

   GL_vlasov->setBoundary(GL_vlasov->fs); 
   GL_fields->solve(GL_vlasov->f0,  GL_vlasov->fs); 
   
   const double rk_0[] = { 0., 0., 0.};
   GL_vlasov->solve(GL_vlasov->getEquationType(), GL_fields, GL_vlasov->fs, GL_vlasov->fss, 1., 0, rk_0);
   
   // copy whole phase space function (waste but starting point) (important due to bounday conditions
   //#pragma omp parallel for, collapse private(n) 
    n = 0;
    for(int s = NsLlD; s <= NsLuD; s++) { for(int m   = NmLlD   ; m   <= NmLuD   ; m++  ) { for(int z = NzLlD; z <= NzLuD; z++) {
    for(int y_k = (GL_includeZF ? 0 : 1); y_k <= NkyLuD-1; y_k++) {   // iteratte from NkyLlD if Zonal Flow is included
    for(int x = NxLlD; x <= NxLuD; x++) { for(int v = NvLlD       ; v <= NvLuD; v++) { 

       y_F1[n++] = fss[s][m][z][y_k][x][v];

   }}} }}}
 
   VecRestoreArrayRead(Vec_x, (const PetscScalar **) &x_F1);
   VecRestoreArray    (Vec_y, (      PetscScalar **) &y_F1);

   } ((A6zz) GL_vlasov->fs, (A6zz) GL_vlasov->fss);
   
   return 0; // return 0 (success) required for PETSc
}

CComplex* PETScMatrixVector::getCreateVector(Grid *grid, Vec &Vec_x) {

    CComplex *xp;

    VecCreateMPI(MPI_COMM_WORLD, grid->getLocalSize(), grid->getGlobalSize(), &Vec_x);
    VecAssemblyBegin(Vec_x);
    VecAssemblyEnd(Vec_x);

    VecGetArray    (Vec_x, (PetscScalar **) &xp);

    return xp;
}
 


