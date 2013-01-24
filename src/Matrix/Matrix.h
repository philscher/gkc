/*
 * =====================================================================================
 *
 *       Filename: Matrix.h
 *
 *    Description: Matrix Interface to access various libraries
 *
 *         Author: Paul P. Hilscher (2012), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __GKC_MATRIX_H__
#define __GKC_MATRIX_H__

#include "Parallel/Parallel.h"

/**
*
*  @brief Matrix interface 
*
*   
**/
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

#endif // __GKC_MATRIX_H__
