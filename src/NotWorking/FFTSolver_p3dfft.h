
#include "config.h"
#include "Global.h"



#ifndef P3DFFT_H
#define P3DFFT_H


#include "FFTSolver.h"

#include <fftw3-mpi.h>




class FFTSolver_p3dfft : public FFTSolver
{
  fftw_plan plan_Phi1DAvrgForward;
  fftw_plan plan_Phi1DAvrgBackward;
  
  public:
       FFTSolver_p3dfft(Setup *setup, Parallel *parallel);
       virtual ~FFTSolver_p3dfft();
    
       int solve(const int FFTtype, const int direction, const int nstacked=1);
       std::string getLibraryName();
       static int    getDecomposition();
       
};


#endif // P3DFFT_H


