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

#ifndef FFTW2_MPI_H
#define FFTW2_MPI_H

#include "FFTSolver.h"
#include "Geometry/GeometryShear.h"
#include "Geometry/Geometry2D.h"

// put into namespace, otherwise conflict with fftw3 (and MPI versions)
//namespace Helios_fftw2 {
//};



class FFTSolver_fftw2_mpi : public FFTSolver
{
  cmplxd *data_1, *data_2, *data_3, *data_4, *data_5, *data_6, *data_7, *data_8, *data_9;

  // memory space for X Transforms
  cmplxd *data_X_rIn_kOut, *data_X_kIn_rOut;  
  cmplxd *data_Y_rIn_kOut, *data_Y_kIn_rOut;  
 
  // Do not included fftw2 structs here, otherwise it will collide with fftw3
  
  int perf_flag;
  Array3z AA_xyz;


void copyNyquiest(Array4z f);
int mapAAY(const int y_k) ;


  public:
    FFTSolver_fftw2_mpi(Setup *setup, Parallel *parallel, Geometry<GKC_GEOMETRY> *geo);
    ~FFTSolver_fftw2_mpi();

    int solve(const int FFTtype, const int direction, const int nstacked=1);
    std::string getLibraryName();
    static int    getDecomposition();

    Array3z  multiply(Array3z &A, Array3z &B, Array3z  &R);

    virtual void printOn(ostream &output) const {
         FFTSolver::printOn(output);
         output   << "FFTSolver  |  using FFTW2-MPI";
         
        }
};



#endif // FFTW2_MPI_H

//#endif // FFTW2_MPI


