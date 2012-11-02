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

#include "FFTSolver/FFTSolver.h"


/**
*   @brief Interface for the fftw Fast-Fourier Solver
*
*   Implementes interface to the fftw Fast-Fourier Transform
*   solver, see (<a href="www.fftw.org">fftw homepage</a>).
*   <www.fftw.org>
*
*
*
**/
class FFTSolver_fftw3 : public FFTSolver
{
   std::string wisdom, plan;

   /**
   *   Arrays for complex-to-complex transform for X-direction
   **/ 
   CComplex *data_X_rIn ,     *data_X_rOut;  
   CComplex *data_X_Transp_1, *data_X_Transp_2; 

   CComplex *data_kXOut, *data_kXIn;

   int AA_NkyLD,    ///< use for anti-alias multiplication
       AA_NyLD;     ///< use for anti-alias multiplication

   void transpose    (int Nx, int Ny, int Nz, int Nq, 
                      CComplex In[Nq][Nz][Ny][Nx], CComplex OutT[Nx][Ny][Nz][Nq]);
   void transpose_rev(int Nx, int Ny, int Nz, int Nq, 
                      CComplex In[Nx][Ny][Nz][Nq], CComplex OutT[Nq][Nz][Ny][Nx]);

  public:
   /**
   *   @brief the contstructor
   *
   *   Setup accepts
   *     FFTW.Plan   = { "Measure", "Estimate", "Exhausting" }
   *     FFTW.Wisdom = { "", "link to file" }
   *
   *     @todo allow file link to plan 
   *
   **/
   FFTSolver_fftw3(Setup *setup, Parallel *parallel, Geometry *geo);

   /**
   *   @brief destructor
   *
   *
   **/
   ~FFTSolver_fftw3();

   /**
   *  @brief write some implementation details
   *
   *
   **/
   void solve(const FFT_Type type, const FFT_Sign sign, void *in=nullptr, void *out=nullptr);
   
   void multiply(const CComplex A[NkyLD][NxLD], const CComplex B[NkyLD][NxLD],
                       CComplex R[NkyLD][NxLD]);
   


   virtual void printOn(std::ostream &output) const;

   virtual void initData(FileIO *fileIO) {};
   virtual void writeData(Timing *timing)  {};
   virtual void closeData()  {};
   
   std::string getLibraryName();

};


#endif // __FFTW3_H


#endif // FFTW3


