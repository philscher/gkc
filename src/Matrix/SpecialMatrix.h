/*
 * =====================================================================================
 *
 *       Filename:  SpecialMatrix.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/11/2012 09:20:39 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */


#ifndef MATRIX_FD_STENCIL
#define MATRIX_FD_STENCIL

#include "Geometry.h"

#include "MatrixPETSc.h"

#include <functional>

class Matrix_FD_Stencil
{
 public:
  
  static void LaplaceOp(Mat &A, Geometry *geo, int y_k, int z, bool periodic=false) {
  
    int idx_l, idx_u;
    
    PetscErrorCode ierr = MatGetOwnershipRange(A,&idx_l,&idx_u);
        
    for (int idx=idx_l; idx<idx_u; idx++) {
     
      int x = idx + NxGlD  ;
      
      // use 2nd order Central Difference stencil h^{-2} [ 1 -2 1]
      const double h = 2 * pow2(dx);
      
      // partial_xx 
      
      // use 4th order Central Difference stencil (1./12 h^2) [ -1 16 -30 16 -1 ]
      const double h_2   = 12. * pow2(dx);
                    
      
      const CComplex val_dx2_diag = geo->g_xx(x  ,z) * ( 30./h_2);
      const CComplex val_dx2_p1   = geo->g_xx(x+1,z) * (-16./h_2);
      const CComplex val_dx2_m1   = geo->g_xx(x-1,z) * (-16./h_2);
      const CComplex val_dx2_p2   = geo->g_xx(x+2,z) * (  1./h_2);
      const CComplex val_dx2_m2  =  geo->g_xx(x-2,z) * (  1./h_2);
                
      ierr = MatSetValue(A, idx, idx, val_dx2_diag, INSERT_VALUES);
                
      if ( (x < NxGuD-1)) ierr = MatSetValue(A, idx, idx+1, val_dx2_p1, INSERT_VALUES);
      if ( (x > NxGlD+1)) ierr = MatSetValue(A, idx, idx-1, val_dx2_m1, INSERT_VALUES);
      if ( (x < NxGuD-2)) ierr = MatSetValue(A, idx, idx+2, val_dx2_p2, INSERT_VALUES);
      if ( (x > NxGlD+2)) ierr = MatSetValue(A, idx, idx-2, val_dx2_m2, INSERT_VALUES);
                 
      // partial_x ky 
      // use 4th order Central Difference stencil (1./12 h^2) [ -1 16 -30 16 -1 ]
      
      const double h_1   = 12. * pow2(dx);
                    
      const CComplex val_dx1_p1  = geo->g_xy(x+1, z) * (  8./h_1);
      const CComplex val_dx1_m1  = geo->g_xy(x-1, z) * (- 8./h_1);
      const CComplex val_dx1_p2  = geo->g_xy(x+2, z) * (- 1./h_1);
      const CComplex val_dx1_m2  = geo->g_xy(x-2, z) * (  1./h_1);
                
      if ( (x < NxGuD-1)) ierr = MatSetValue(A, idx, idx+1, val_dx1_p1, ADD_VALUES);
      if ( (x > NxGlD+1)) ierr = MatSetValue(A, idx, idx-1, val_dx1_m1, ADD_VALUES);
      
      if ( (x < NxGuD-2)) ierr = MatSetValue(A, idx, idx+2, val_dx1_p2, ADD_VALUES);
      if ( (x > NxGlD+2)) ierr = MatSetValue(A, idx, idx-2, val_dx1_m2, ADD_VALUES);
                    
      // Finally (BUG y_k)
      const double val_dy2_diag = - geo->g_yy(x,z) * y_k;
      
      MatSetValue(A, idx, idx, val_dy2_diag, ADD_VALUES);
      
    }

    
  };
    
  static Complex D_2O2(Mat &A, std::function<double (double)> g, int order = 4 , bool periodic=false) {
    
    int idx_l, idx_u;
    
    PetscErrorCode ierr = MatGetOwnershipRange(A,&idx_l,&idx_u);
    
    for (int idx=idx_l; idx<idx_u; idx++) {
     
      int x   = idx + NxGlD  ;
           
      // use 2nd order Central Difference stencil h^{-2} [ 1 -2 1]
      const double h   = 2 * pow2(dx);
               
      if(order == 2) {

        const CComplex val_p1  =  g(x+1) * ( 1./h);
        const CComplex val_m1  =  g(x-1) * (-1./h);
        
        // why we have to set sub / sup when Matrix is symmetric ? Recheck
        if ( (x < NxGuD-1)) ierr = MatSetValue(A, idx, idx+1, val_p1, INSERT_VALUES);
        if ( (x > NxGlD+1)) ierr = MatSetValue(A, idx, idx-1, val_m1, INSERT_VALUES);
                
        // Matrix entries needs not to be local ... or ? (PETSc does take care of ?)
        
        if(periodic == true) {
        
          if ( x == NxGlD) ierr = MatSetValue(A, 0      , NxGuD-3, val_p1, INSERT_VALUES);
          if ( x == NxGuD) ierr = MatSetValue(A, NxGuD-3,       0, val_m1, INSERT_VALUES);
        }

      } 
      
      else if(order == 4) {

        // use 4th order Central Difference stencil (1./12 h^2) [ -1 16 -30 16 -1 ]
        
        const double h   = 12. * pow2(dx);
        
        const CComplex val_p1   = g(x+1) * (  8./h);
        const CComplex val_m1   = g(x-1) * (- 8./h);
        const CComplex val_p2   = g(x+2) * (- 1./h);
        const CComplex val_m2  =  g(x-2) * (  1./h);
                
        if ( (x < NxGuD-1)) ierr = MatSetValue(A, idx, idx+1, val_p1, INSERT_VALUES);
        if ( (x > NxGlD+1)) ierr = MatSetValue(A, idx, idx-1, val_m1, INSERT_VALUES);
                    
        if ( (x < NxGuD-2)) ierr = MatSetValue(A, idx, idx+2, val_p2, INSERT_VALUES);
        if ( (x > NxGlD+2)) ierr = MatSetValue(A, idx, idx-2, val_m2, INSERT_VALUES);

      } else {
      
        check(-1, DMESG("No such order"));
      }
    }
    
    return 0.;
    
  };
    
  
  static Complex D_1O2(Mat &A, std::function<double (double)> g, int order =4 , bool periodic=false) {
    
  
    int idx_l, idx_u;
    
    PetscErrorCode ierr = MatGetOwnershipRange(A,&idx_l,&idx_u);
    
    for (int idx=idx_l; idx<idx_u; idx++) {
      
      int x   = idx + NxGlD  ;
      
      // use 2nd order Central Difference stencil h^{-2} [ 1 -2 1]
      const double h   = pow2(dx);
      
      if(order == 2) {

        const CComplex val_diag = g(x  ) * ( 2./h);

        const CComplex val_p1   = g(x+1) * (-1./h);
        const CComplex val_m1   = g(x-1) * (-1./h);

        // Set D [ S_ D S^]
        ierr = MatSetValue(A, idx, idx, val_diag, INSERT_VALUES);

        // why we have to set sub / sup when Matrix is symmetric ? Recheck
        if ( (x < NxGuD-1)) ierr = MatSetValue(A, idx, idx+1, val_p1, INSERT_VALUES);
        if ( (x > NxGlD+1)) ierr = MatSetValue(A, idx, idx-1, val_m1, INSERT_VALUES);
        
        // Matrix entries needs not to be local ... or ? (PETSc does take care of ?)
                  
        if(periodic == true) {
        
          if ( x == NxGlD) ierr = MatSetValue(A, 0      , NxGuD-3, val_p1, INSERT_VALUES);
          if ( x == NxGuD) ierr = MatSetValue(A, NxGuD-3,       0, val_m1, INSERT_VALUES);
        }
      } 
      
      else if(order == 4) {
        
        // use 4th order Central Difference stencil (1./12 h^2) [ -1 16 -30 16 -1 ]
        
        const double h   = 12. * pow2(dx);
        
        const CComplex val_diag = g(x  ) * ( 30./h);
        const CComplex val_p1   = g(x+1) * (-16./h);
        const CComplex val_m1   = g(x-1) * (-16./h);
        const CComplex val_p2   = g(x+2) * (  1./h);
        const CComplex val_m2  =  g(x-2) * (  1./h);
                
        ierr = MatSetValue(A, idx, idx, val_diag, INSERT_VALUES);
        
        if ( (x < NxGuD-1)) ierr = MatSetValue(A, idx, idx+1, val_p1, INSERT_VALUES);
        if ( (x > NxGlD+1)) ierr = MatSetValue(A, idx, idx-1, val_m1, INSERT_VALUES);
        
        if ( (x < NxGuD-2)) ierr = MatSetValue(A, idx, idx+2, val_p2, INSERT_VALUES);
        if ( (x > NxGlD+2)) ierr = MatSetValue(A, idx, idx-2, val_m2, INSERT_VALUES);
      } 
      else  {
      
        check(-1, DMESG("No such order"));
      }
      
    }
    
    return 0.;

  };

  static void setLaplace_x_ky(Mat &A, Geometry *geo, double ky2, std::string order = "CD-4", bool periodic_in_x=false) {
  
    // b = k_x^2 + k_y^2 = - \partial_x^2 + k_y^2
    
    int idx_l, idx_u;
    
    PetscErrorCode ierr = MatGetOwnershipRange(A,&idx_l,&idx_u);
    
    for (int idx = idx_l; idx < idx_u; idx++) {
     
      int x   = idx + NxGlD  ;
      
      if (order == "CD-2") {
      
        // use 2nd order Central Difference stencil h^{-2} [ 1 -2 1]
        
        const double h   = pow2(dx);
        
        const CComplex val_diag = (2./h - ky2);
        const CComplex val_pm1  = -1./h;

        // Set D [ S_ D S^]
        
        ierr = MatSetValue(A, idx, idx, val_diag, INSERT_VALUES);

        // why we have to set sub / sup when Matrix is symmetric ? Recheck
        if ( (x < NxGuD-1)) ierr = MatSetValue(A, idx, idx+1, val_pm1, INSERT_VALUES);
        if ( (x > NxGlD+1)) ierr = MatSetValue(A, idx, idx-1, val_pm1, INSERT_VALUES);
                
        //if ( x == NxGlD) ierr = MatSetValue(A, idx, ip1, 0., INSERT_VALUES);
        //if ( x == NxGuD) ierr = MatSetValue(A, idx, im1, 0., INSERT_VALUES);
  
        // Matrix entries needs not to be local ... or ? (PETSc does take care of ?)
        if(periodic_in_x) {
        
          if ( x == NxGlD) ierr = MatSetValue(A, 0      , NxGuD-3, val_pm1, INSERT_VALUES);
          if ( x == NxGuD) ierr = MatSetValue(A, NxGuD-3,       0, val_pm1, INSERT_VALUES);
          
        };
                
        
      } 
      
      else if(order == "CD-4") {

      
        // use 4th order Central Difference stencil (1./12 h^2) [ -1 16 -30 16 -1 ]
        
        const double h   = 12. * pow2(dx);
                    
        
        const CComplex val_diag = (30./h - ky2);
        const CComplex val_pm1  = - 16./h;
        const CComplex val_pm2  = 1./h;
                
        
        ierr = MatSetValue(A, idx, idx, val_diag, INSERT_VALUES);
                
        
        if ( (x < NxGuD-1)) ierr = MatSetValue(A, idx, idx+1, val_pm1, INSERT_VALUES);
        if ( (x > NxGlD+1)) ierr = MatSetValue(A, idx, idx-1, val_pm1, INSERT_VALUES);
        
        if ( (x < NxGuD-2)) ierr = MatSetValue(A, idx, idx+2, val_pm2, INSERT_VALUES);
        if ( (x > NxGlD+2)) ierr = MatSetValue(A, idx, idx-2, val_pm2, INSERT_VALUES);

        
      } else if(order == "CD-6") {

        // use 6th order Central Difference stencil (1./180 h^2) [ 2 -27 270 -490 270 -27 2 ]
        
        const double h   = 180. * pow2(dx);
                    
        
        const CComplex val_diag =  (490./h - ky2);
        const CComplex val_pm1  = - 270./h;
        const CComplex val_pm2  = +  27./h;
        const CComplex val_pm3  = -   2./h;
                
        ierr = MatSetValue(A, idx, idx, val_diag, INSERT_VALUES);
                
        if ( (x < NxGuD-1)) ierr = MatSetValue(A, idx, idx+1, val_pm1, INSERT_VALUES);
        if ( (x > NxGlD+1)) ierr = MatSetValue(A, idx, idx-1, val_pm1, INSERT_VALUES);
                    
        if ( (x < NxGuD-2)) ierr = MatSetValue(A, idx, idx+2, val_pm2, INSERT_VALUES);
        if ( (x > NxGlD+2)) ierr = MatSetValue(A, idx, idx-2, val_pm2, INSERT_VALUES);
        
        if ( (x < NxGuD-3)) ierr = MatSetValue(A, idx, idx+3, val_pm2, INSERT_VALUES);
        if ( (x > NxGlD+3)) ierr = MatSetValue(A, idx, idx-3, val_pm2, INSERT_VALUES);

      } else check(-1, DMESG("No such Order"));
      
    }
  
  };


};



#endif // MATRIX_FD_STENCIL
