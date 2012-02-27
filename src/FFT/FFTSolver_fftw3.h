/*
 * =====================================================================================
 *
 *       Filename : FFTSolver_fftw3.h
 *
 *    Description : Implementation of FFTSolver using fftw-3(+mpi/threads) 
 *                  from www.fftw.org
 *        Author  : Paul P. Hilscher (2012), 
 *
 * =====================================================================================
 */


#include "config.h"
#include "Global.h"

#ifdef FFTW3

#ifndef __FFTW3_H
#define __FFTW3_H

#include "FFTSolver.h"
#include "Geometry/GeometryShear.h"
#include "Geometry/Geometry2D.h"



class FFTSolver_fftw3 : public FFTSolver
{
  /**
   *   Arrays for real-complex transform
   */ 
  double *data_Y_rIn, *data_Y_rOut;
  cmplxd *data_Y_kIn, *data_Y_kOut;

  cmplxd *data_X_rIn , *data_X_rOut;  
  cmplxd *data_X_kIn , *data_X_kOut; 

  /**
   *    Anti-Aliased Fourier transform used by
   *    multiply.
   *
   * */
  Range AA_RkyLD, AA_RyLD;
  int AA_NyLD, AA_NyLlD, AA_NyLuD, AA_NkyLD, AA_NkyLlD, AA_NkyLuD;
  Array4d AA_rYIn, AA_rYOut;
  Array4z AA_kYOut, AA_kYIn;
  
  public:
    FFTSolver_fftw3(Setup *setup, Parallel *parallel, Geometry<HELIOS_GEOMETRY> *geo);
    ~FFTSolver_fftw3();

    int solve(int FFTtype, int direction, const int N);
    std::string getLibraryName();
    static int    getDecomposition();
    Array3z  multiply(Array3z &A, Array3z &B, Array3z  &R);
    
    virtual void printOn(ostream &output) const;

};


#endif // __FFTW3_H


#endif // FFTW3


