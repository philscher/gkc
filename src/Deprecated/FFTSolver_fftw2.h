/*
 * =====================================================================================
 *
 *       Filename:  m_FFTSolver_fftw2-mpi.h
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


///////////////////////////////////////// FFTW-MPI /////////////////////////////////////////////////////


#include "config.h"
#include "Global.h"

#ifdef FFTW2

#ifndef FFTSOLVER_FFTW2_H
#define FFTSOLVER_FFTW2_H


#include "FFTSolver.h"
#include "Geometry/GeometryShear.h"

namespace Helios_fftw2 {
#include <rfftw.h>
};

class FFTSolver_fftw2 : public FFTSolver
{
  double *data_1, *data_2, *data_3, *data_4;
  /*
   * 
  fftw_plan plan_Phi1DAvrgForward;
  fftw_plan plan_Phi1DAvrgBackward;
  fftw_plan plan_Phi3DBackward;
   */

  rfftwnd_plan plan_Phi2DManyForward;
  rfftwnd_plan plan_Phi2DManyBackward;
  rfftwnd_plan plan_Phi3DForward;
  rfftwnd_plan plan_Phi3DBackward;

  public:
    FFTSolver_fftw2(Setup *setup, Parallel *parallel, Geometry *geo);
    ~FFTSolver_fftw2();

    int solve(int FFTtype, int direction);
    std::string getLibraryName();
    static int    getDecomposition();
};







#endif // FFTW2_MPI_H

#endif // FFTW2


