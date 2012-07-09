/*
 * =====================================================================================
 *
 *       Filename:  m_FFTSolver_fftw3_fftw.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/29/2010 12:51:15 PM
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

#ifndef FFTSOLVER_FFTW3_DST_H
#define FFTSOLVER_FFTW3_DST_H


#include "FFTSolver.h"
#include "Geometry/GeometryShear.h"
#include "Geometry/Geometry2D.h"

namespace Helios_fftw3 {
#include <fftw3.h>
};


class FFTSolver_fftw3_dst : public FFTSolver
{
  Helios_fftw3::fftw_plan plan_Phi1DAvrgForward;
  Helios_fftw3::fftw_plan plan_Phi1DAvrgBackward;

  Helios_fftw3::fftw_plan plan_Phi3DForward;
  Helios_fftw3::fftw_plan plan_Phi3DBackward;

  Helios_fftw3::fftw_plan plan_Phi2DManyForward;
  Helios_fftw3::fftw_plan plan_Phi2DManyBackward;

  /** Discrete Transformation type */
  int DTtype;

  enum { HELIOS_FFT_NONE=0, HELIOS_FFT_DFT=1, HELIOS_FFT_DST=2, HELIOS_FFT_DCT=4};
    
  Range Rk2xL2; 
  
  public:
    FFTSolver_fftw3_dst(Setup *setup, Parallel *parallel, Geometry<HELIOS_GEOMETRY> *geo);
    ~FFTSolver_fftw3_dst();

    int solve(const int FFTtype, const int direction, const int nstacked=1);
    std::string getLibraryName();
    static int    getDecomposition();
       /* 
    virtual void printOn(ostream &output) const {
         FFTSolver::printOn(output);
         output   << "FFTSolver  |  using FFTW3 with ";
         

         // we can also use anti-aliasing 
         if     (DTtype == HELIOS_FFT_DFT)  output << " DFT in X (Nx zero padding)" << std::endl;
         else if(DTtype == HELIOS_FFT_DST)  output << " DST in X                  " << std::endl;
         else if(DTtype == HELIOS_FFT_DCT)  output << " DCT in X                  " << std::endl;
         else check(-1, DMESG("No such type"));
        }
        * */ 
    Array3z  multiply(Array3z &A, Array3z &B, Array3z  &R) { 
      check(-1, DMESG("Please help : Not implemented"));
    }
};

#endif // FFTSOLVER_FFTW3__DST_H


#endif // FFTW3
