/*
 * =====================================================================================
 *
 *       Filename:  PETScMatrixVector.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  02/02/2012 11:50:28 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include "PETScMatrixVector.h"

Vlasov *GL_vlasov;
Fields *GL_fields;

int GL_iter;

PETScMatrixVector::PETScMatrixVector(Vlasov *vlasov, Fields *fields)
{

            GL_vlasov = vlasov;
            GL_fields = fields;
            GL_iter   = 0;
};


PetscErrorCode PETScMatrixVector::MatrixVectorProduct(Mat A, Vec Vec_x, Vec Vec_y) 

{

  // get Constext
    std::cout << "\r"   << "Iteration  : " << GL_iter++ << std::flush;

    cmplxd *x_F1, *y_F1; 
   
 
    VecGetArrayRead(Vec_x, (const PetscScalar **) &x_F1);
    VecGetArray    (Vec_y, (PetscScalar **) &y_F1);

    // copy whole phase space function (waste but starting point) (important due to bounday conditions
   // we can built wrapper around this and directly pass it
   for(int x = NxLlD, n = 0; x <= NxLuD; x++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { for(int z = NzLlD; z <= NzLuD; z++) {
   for(int v = NvLlD       ; v <= NvLuD; v++) { for(int m   = NmLlD ; m   <= NmLuD ; m++  ) { for(int s = NsLlD; s <= NsLuD; s++) {
   	
	        GL_vlasov->fs(x,y_k,z,v,m,s) = x_F1[n++];

   }}} }}}

   GL_vlasov->setBoundary(GL_vlasov->fs, BOUNDARY_CLEAN); 
   GL_fields->solve(GL_vlasov->f0,  GL_vlasov->fs); 
   GL_vlasov->solve(GL_vlasov->getEquationType() + "_Eigenvalue", GL_fields, GL_vlasov->fs, GL_vlasov->fss, 0., 0);
   
    // copy whole phase space function (waste but starting point) (important due to bounday conditions
   for(int x = NxLlD, n = 0; x <= NxLuD; x++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { for(int z = NzLlD; z <= NzLuD; z++) {
   for(int v = NvLlD       ; v <= NvLuD; v++) { for(int m   = NmLlD ; m   <= NmLuD ; m++  ) { for(int s = NsLlD; s <= NsLuD; s++) {

	    y_F1[n++] = GL_vlasov->fss(x,y_k,z,v,m,s); 

   }}} }}}
 
   VecRestoreArrayRead(Vec_x, &x_F1);
   VecRestoreArray    (Vec_y, &y_F1);


   return 0; // return 0 (success) required for PETSc
}

cmplxd* PETScMatrixVector::getCreateVector(Grid *grid, Vec &Vec_x) {

    cmplxd *xp;

    VecCreateMPI(MPI_COMM_WORLD, grid->getLocalSize(), grid->getGlobalSize(), &Vec_x);
    VecAssemblyBegin(Vec_x);
    VecAssemblyEnd(Vec_x);

    VecGetArray    (Vec_x, &xp);

    return xp;
}
 


