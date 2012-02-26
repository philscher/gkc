/*
 * =====================================================================================
 *
 *       Filename:  m_FFTSolver_fftw3-mpi.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/29/2010 12:52:03 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
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
                
  cmplxd *data_Y_kIn, *data_Y_kOut;
  double *data_Y_rIn, *data_Y_rOut;

  //cmplxd *data_X_rIn_kOut, *data_X_kIn_rOut;  
  
  cmplxd *data_X_kOut, *data_X_kIn;  
  cmplxd *data_X_rIn, *data_X_rOut;  

  int perf_flag;
 
  Range RkyLD_AA;

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


