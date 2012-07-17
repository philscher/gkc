/*
 * =====================================================================================
 *
 *       Filename: FFTSolver.cpp
 *
 *    Description: Interface for various (Fast) Fourier Solvers
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include "FFTSolver.h"


void FFTSolver::setNormalizationConstants() {


	// find out normalization (assume constant normalization accors boundaries)
    // Real -> Complex (Fourier sapce) transform
    rYIn(NxLlD, RyLD, NzLlD, 1) = 1.;
		
    solve(FFT_Y, FFT_FORWARD, NxLD * NzLD);
	Norm_Y_Forward = real(kYOut(NxLlD, 0, NzLlD, 1));		


    // Complex (Fourier sapce) -> Real transform
    kYIn(RxLD, RkyLD, RzLD, 1) = 0.; kYIn(NxLlD, NkyLlD, NzLlD, 1) = 1.;
	solve(FFT_Y, FFT_BACKWARD,  NxLD * NzLD);
	Norm_Y_Backward = rYOut(NxLlD, NyLlD, NzLlD, 1);
        
    // Real -> Complex (Fourier sapce) transform
	rXIn(RxLD, NkyLlD, NzLlD, 1) = 1.;
	solve(FFT_X, FFT_FORWARD, NkyLD * NzLD);
	Norm_X_Forward = (K1xLlD == 0) ? real(kXOut(0, NkyLlD, NzLlD, 1)) : 0.;		


    // Complex (Fourier sapce) -> Real transform
    kXIn(Rk1xL, RkyLD, RzLD, 1) = 0.; if(K1xLlD == 0) kXIn(0, NkyLlD, NzLlD, 1) = 1.;
	solve(FFT_X, FFT_BACKWARD,  NkyLD * NzLD);
	Norm_X_Backward = real(rXOut(NxLlD, NkyLlD, NzLlD, 1));
		
    // broadcast normalization to all nodes
    parallel->barrier();
    parallel->send(Norm_Y_Forward, parallel->Coord(DIR_ALL) == 0); parallel->send(Norm_Y_Backward, parallel->Coord(DIR_ALL) == 0);
    parallel->send(Norm_X_Forward, parallel->Coord(DIR_ALL) == 0); parallel->send(Norm_X_Backward, parallel->Coord(DIR_ALL) == 0);
    parallel->barrier();

};



FFTSolver::FFTSolver(Setup *setup, Parallel *_parallel, Geometry<HELIOS_GEOMETRY> *_geo, double _Norm_XYZ, double _Norm_XY, double _Norm_X, double _Norm_Y) :
      parallel(_parallel), Norm_XYZ(_Norm_XYZ), Norm_XY(_Norm_XY), Norm_X(_Norm_X), Norm_Y(_Norm_Y), geo(_geo)
   {
     // is this really necessary ?
     flags = FFT_X | FFT_Y;
     if(setup->get("Vlasov.useAA" , 0) == 1)   flags |= FFT_AA;
   };


FFTSolver::~FFTSolver() {};

