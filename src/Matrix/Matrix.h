/*
 * =====================================================================================
 *
 *       Filename: Matrix.h
 *
 *    Description: Defined Matrix Interface
 *
 *         Author: Paul P. Hilscher (2012), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __GKC_Matrix_H
#define __GKC_Matrix_H

#include <petscksp.h>

#include "Parallel/Parallel.h"

/**
 *
 *   Class Defined Matrix implementation using PETSc (http://www.mcs.anl.gov/petsc/)
 *
  *  From the PETSc documenation (do we take care of this ? ) : 
  *     
  *     Currently, all PETSc parallel matrix formats are partitioned by
  *     contiguous chunks of rows across the processors.  Determine which
  *     rows of the matrix are locally owned.
 *   
 *
 *
 */
class Matrix
{

  public:
  virtual Matrix(int _n_local, int _n_global, int dir, Parallel *_parallel);
  virtual void checkAndAssemble() ;
  virtual void setValue(int col, int row, Complex value);
  virtual void assemble() ;
  virtual void setZero() ;
  
  
  //// Operator overloading /////////////
  virtual friend Matrix& operator+=(Matrix &A, Complex a) ;
  virtual friend Matrix& operator+=(Matrix &A, Matrix &B) ;
  virtual friend Matrix& operator*(Matrix &A, Matrix &B) ;
  virtual friend Matrix& operator+(Matrix &A, Complex a) ;
  virtual friend Matrix& operator*(Complex a, Matrix &A) ;
  virtual friend Matrix& operator*(double &a, Matrix &A) ;
  virtual friend Matrix& operator*(Matrix &A, double &a) ;
  virtual friend Matrix& operator=(double a) {};
  inline Complex& Matrix::operator() (int row, int col)
  virtual void addDiagonal(Complex a) {
  virtual void reduce(int dir);

};


#endif // __GKC_MatrixP_H
