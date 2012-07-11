/*
 * =====================================================================================
 *
 *       Filename: Fourier.h
 *
 *    Description: Helper functions for Fourier space
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __FOURIER_H
#define __FOURIER_H

#include "Global.h"
#include "Setup.h"
#include "Grid.h"
#include "FFTSolver.h"
#include "Geometry.h"
#include "GeometryShear.h"
#include "Geometry2D.h"

#include "Plasma.h"
#include "SpecialMath.h"

/**
*
*    @brief Fourier3D
*
*
*/ 
class Fourier3D {

    std::vector<int> suppressModeX, suppressModeY, suppressModeZ;
    std::vector<int> convolveModeX, convolveModeY, convolveModeZ;


double epsilon_0, sigma;
  protected:

  FFTSolver *fft;

  public:
  
  Fourier3D(Setup *setup, Grid *grid, FFTSolver *fftsolver);
  virtual  ~Fourier3D();

  std::string getSolverInfo();
  
    /*!  We can suppress various modes, this is set in the setup of the fields.
     *   and sets the Fourier mode to zero. Move to FFT solver.
     */
    int suppressModes(Array4z k2Out, const int field=1);


};


#endif // __Fourier_H
