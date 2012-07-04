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

#include "MatrixSolver.h"
#include "MatrixPETSc.h"


static char FieldsHermite_help[] = "Help for PETSc Interface not available, please look up gkc & PETSc manual.";



/**
 *
 *    Used implementation according to 
 *    
 *    Solving Fields equations using Hermite interpolation method (in x-direction).
 *    Based according to
 *      
 *          (Goerler, PhD, 2008) and (Lapillione, PhD, 2009)
 *          which is also used in the GENE code (http://www.ipp.mpg.de/~fsj/)
 *
 *    Gyro-averagin is done by using an interpolation matrix.
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
 *          \Gamma_0 = \int_0^\infty  J^2 \exp{(-\mu)} d\mu  .
 *
 *    This integral can be calculated accuratley by using high-order Gauss-Laguerre
 *    quadrature (a temporary gyro-averaing matrix is setup for this particular point).
 *
 *    Geometrical effects are included by linearizing the metric. As discussed by Goerler
 *    & Lapillone.
 *
 *    More advanced model may be required.
 *
 *
 *
 *
 *
 * */
class FieldsHermite : public FieldsFFT {

  /**
   *    Please Document Me !
   *
   **/
  Array4z              B4;
  /**
   *    Please Document Me !
   *
   **/
  Vec                  GyroVector_X, GyroVector_XAvrg;
  /**
   *    Please Document Me !
   *
   **/
  blitz::Array<Matrix*      ,4>  GyroMatrix;    
  /**
   *    Please Document Me !
   *
   **/
  blitz::Array<MatrixSolver*,2>  MatrixPoissonSolverLHS;


protected:   
  
  /**
   *    Please Document Me !
   *
   **/
  virtual Array4z solveFieldEquations(Array4z fields, Timing timing) ;// override;
  /**
   *    Please Document Me !
   *
   **/
  virtual Array3z solvePoissonEquation(Array3z rho, Timing timing)   ;//override;
  /**
   *    Please Document Me !
   *
   **/
  virtual Array4z g2yroAverage(Array4z A4, int m, int s, const int nField, bool gyroField);// override;
    
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
    cmplxd getElements(const int i, const int n, const double r, const int y_k, const int z);
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
 FieldsHermite(Setup *setup, Grid *grid, Parallel *parallel, FileIO *fileIO, Geometry<HELIOS_GEOMETRY> *geo, FFTSolver *fftsolver);
  /**
   *    Please Document Me !
   *
   **/
 ~FieldsHermite();
 
protected:
  /**
   *    Please Document Me !
   *
   **/
   virtual void printOn(ostream &output) const;

};


#endif // __FIELDS_HERMITE_H





