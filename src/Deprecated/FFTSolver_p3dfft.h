
#include "config.h"
#include "Global.h"



#ifndef P3DFFT_H
#define P3DFFT_H


#include "FFTSolver.h"

#include <fftw3-mpi.h>




/*!
*  For the FFT solver, this class rely on FFTW (www.fftw.org). This library is very efficient
*  and if compiled for parallism (pthread, OpenMP or/and CUDA) it can make us of
*  the avaible CPUs/Cores very well.
*
*  We us FFTW to transform into Fourier space using the real2complex function. This function
*  takes real input and produces complex output. Becausefor real transform half of the
*  values are complex conjugant to each other and reduntant.
*
*  After transforming a 3 dimensions n_x * n_y * n_z array of real variables the output
*  is n_x * n_y * (n_z/2 + 1). We should remember that for r2c/c2r transforms the
*  input array is destroyed !
*
*
*  This class uses p3dfft to perform the parallized Fourier transformation.
*  This is done using a 2D decompostion in X and Y. The Vlasov solver
*  is adapted to p3dfft to follow their order
*
*
*  p3dfft proviedes basically 4 functions:
*       p3dfft_setup : Initialiization, we need to specify the decomposition in
*                      (nx, ny, nz) direction. Because decompoisiton is in
*                      ny and nz (!), we need to initialize as (nz, nx, ny) and
*                      switch the variables accordingly
*       get_dims     : We get the decomposition of the grid in real (conf = 1)
*                      and Fourier space (conf1). We need to keep in mind that
*                      the data in Fourier space is in trnasposed form !
*
*  Note that all the array properties for input/output are hidden inside the r3In 
*  k3Out classes, in face they are
*
*  FFTW     :  r3In is normal; k3Out is normal
*  FFTW-MPI :  r3In is normal; k3Out(y,x,z)
*  P3DFFT   :  r3In(y,x,z),    k3Out(y,x,z) (don't set DSTRIDE)
*
*  all this is hidden inside the GeneralArray storage class, so we can access
*  r3In and k3Out as (x,y,z) !!
*
*
*/ 
class FieldsFFT : public Fields, public Fourier3D {
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


