/*
 * =====================================================================================
 *
 *       Filename: FieldsHermite.h
 *
 *    Description: Interface to Hermite Matrix Solver routine. Solves Field
 *                 equation using Hermite interpolation.
 *
 *         Author: Paul P. Hilscher (2012), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#include "config.h"
#include "Global.h"


#ifndef FIELDS_HERMITE_H
#define FIELDS_HERMITE_H


#include <petscksp.h>

#include "FieldsFFT.h"
#include "external/Array.h"
#include "Matrix/MatrixSolver.h"
#include "Matrix/MatrixPETSc.h"





/**
*
*    @brief Solving Fields equations using Hermite interpolation method
*    
*    Based according to
*      
*          (Goerler, PhD, 2008) and (Lapillione, PhD, 2009)
*          which is also used in the GENE code (http://www.ipp.mpg.de/~fsj/)
*
*    Gyro-averaging is done by using an interpolation matrix.
*
*
*    A(x_i,ky) = A_{x_i} e^{i k_y y}. Values are only available at the grid points at x_i.
*    So we apply Hermite interpolation, by interpolating between x_i and x_{i+1}.
*
*    A'(x,ky)  = Lambda(x) e^{i ky} 
*
*    where Lambda is the piecewise polynomial Hermitian basis.
*
*
*    To solve for the Poisson equation the double gyro-averaging operator needs to be
*    approximated by
*
*      \f[    \Gamma_0 = \int_0^\infty  J^2 \exp{(-\mu)} d\mu  \f]
*
*    This integral can be calculated accurately by using high-order Gauss-Laguerre
*    quadrature (a temporary gyro-averaging matrix is setup for this particular point).
*
*    Geometrical effects are included by linearized the metric. As discussed by Goerler
*    & Lapillone.
*
*
*
**/
class FieldsHermite : public Fields {

  Vec   GyroVector_X    ,  ///< PETSc Vector for field 
        GyroVector_XAvrg;  ///< PETSc Vector for gyro-averaged field
  
  Array4<Matrix*>  GyroMatrix;   ///< Contains Matrix to solve array 
  
  Array2<MatrixSolver*>  MatrixPoissonSolverLHS;

  /**
  *    @brief Setup double-gyroaverage matrix
  *
  **/
  void setDoubleGyroAverageMatrix(Setup *setup);

protected:   
  
  /**
  *    Please Document Me !
  *
  **/
  void solveFieldEquations(const CComplex Q     [Nq][NxLD][NkyLD][Nz],
                                 CComplex Field0[Nq][NxLD][NkyLD][Nz]);
  /**
  *    Please Document Me !
  *
  **/
  void solvePoissonEquation(const CComplex Q     [Nq][NxLD][NkyLD][Nz],
                                  CComplex Field0[Nq][NxLD][NkyLD][Nz]);
  /**
  *    Please Document Me !
  *
  **/
  void gyroAverage(const CComplex In [Nq][NzLD][NkyLD][NxLD], 
                         CComplex Out[Nq][NzLD][NkyLD][NxLD],
                   const int m, const int s, const bool forward, const bool stack=false);
   
  /**
  *
  *
  **/
  void doubleGyroExp(const CComplex In [Nq][NzLD][NkyLD][NxLD], 
                           CComplex Out[Nq][NzLD][NkyLD][NxLD], const int m, const int s) {};
 
    
private:
  /**
  *    Please Document Me !
  *
  **/
  int interpolationOrder;
    
  /**
  *    Please Document Me !
  *
  **/
  Matrix* getGyroAveragingMatrix(const double mu, const int y_k, const int z, int s);

  /**
  *    Please Document Me !
  *
  **/
  CComplex getElements(const int i, const int n, const double r, const int y_k, const int z);

  /**
  *    Please Document Me !
  *
  **/
  double Lambda(const double x, const int n);

 public:

  /**
  *    Please Document Me !
  *
  **/
  FieldsHermite(Setup *setup, Grid *grid, Parallel *parallel, FileIO *fileIO, Geometry *geo);

  /**
  *    Please Document Me !
  *
  **/
 ~FieldsHermite();
  
  /**
  *    @brief calculates the field energy
  *
  **/ 
  void getFieldEnergy(double& phiEnergy, double& ApEnergy, double& BpEnergy);
 
protected:
  /**
  *    Please Document Me !
  *
  **/
  virtual void printOn(std::ostream &output) const;

};


#endif // __FIELDS_HERMITE_H





