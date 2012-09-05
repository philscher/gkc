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
 *         Author:  Paul P. Hilscher 
 *        Company:  
 *
 * =====================================================================================
 */
#include "config.h"
#include "Global.h"

#ifdef FFTW3

#ifndef FFTSOLVER_FFTW3_H
#define FFTSOLVER_FFTW3_H


#include "FFTSolver.h"
#include "Geometry/GeometryShear.h"
#include "Geometry/Geometry2D.h"

namespace Helios_fftw3 {
#include <fftw3.h>
}


class FFTSolver_fftw3 : public FFTSolver
{
  Helios_fftw3::fftw_plan plan_Phi1DAvrgForward;
  Helios_fftw3::fftw_plan plan_Phi1DAvrgBackward;

  Helios_fftw3::fftw_plan plan_Phi3DForward;
  Helios_fftw3::fftw_plan plan_Phi3DBackward;

  Helios_fftw3::fftw_plan plan_Phi2DManyForward;
  Helios_fftw3::fftw_plan plan_Phi2DManyBackward;

  bool isDST;

  public:
    FFTSolver_fftw3(Setup *setup, Parallel *parallel, Geometry<HELIOS_GEOMETRY> *geo);
    ~FFTSolver_fftw3();

    int solve(const int FFTtype, const int direction, const int nstacked=1);
    std::string getLibraryName();
    static int    getDecomposition();
    Array3z  multiply(Array3z &A, Array3z &B, Array3z  &R) { 
      check(-1, DMESG("Please help : Not implemented"));
    }
};

#endif // FFTSOLVER_FFTW3_H


#endif // FFTW3
