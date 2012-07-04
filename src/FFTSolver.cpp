/*
 * =====================================================================================
 *
 *       Filename:  FFTSolver.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  02/23/2012 12:51:16 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include "FFTSolver.h"


void FFTSolver::checkNormalization() {


		// find out normalization (assume constant normalization accors boundaries)
        // Real -> Complex (Fourier sapce) transform
		rYIn(NxLlD, RyLD, NzLlD, 1) = 1.;
        //std::cout << " rYIn : " <<  rYIn(NxLlD, RyLD, NzLlD, 1) << std::endl;
		solve(FFT_Y, FFT_FORWARD, NxLD * NzLD);
		Norm_Y_Forward = real(kYOut(NxLlD, 0, NzLlD, 1));		
        //std::cout << " kYOut :  " << kYOut(NxLlD, RkyLD, NzLlD, 1) << std::endl;


        // Complex (Fourier sapce) -> Real transform
        kYIn(RxLD, RkyLD, RzLD, 1) = 0.; kYIn(NxLlD, NkyLlD, NzLlD, 1) = 1.;
        //kYIn = kYOut;     
		solve(FFT_Y, FFT_BACKWARD,  NxLD * NzLD);
        //std::cout << " rYOut : " <<  rYOut(NxLlD, RyLD, NzLlD, 1) << std::endl;
		Norm_Y_Backward = rYOut(NxLlD, NyLlD, NzLlD, 1);
		//std::cout << "Normalization is : " << Norm_Y_Forward << " " << Norm_Y_Backward << " Full : " << Norm_Y << std::endl;     	

        
        
        // Real -> Complex (Fourier sapce) transform
		rXIn(RxLD, NkyLlD, NzLlD, 1) = 1.;
        //std::cout << " rXIn : " <<  rXIn(RxLD, NkyLlD, NzLlD, 1) << std::endl;
		solve(FFT_X, FFT_FORWARD, NkyLD * NzLD);
		Norm_X_Forward = (K1xLlD == 0) ? real(kXOut(0, NkyLlD, NzLlD, 1)) : 0.;		
        //std::cout << " kXOut :  " << kXOut(Rk1xL, NkyLlD, NzLlD, 1) << std::endl;


        // Complex (Fourier sapce) -> Real transform
        kXIn(Rk1xL, RkyLD, RzLD, 1) = 0.; if(K1xLlD == 0) kXIn(0, NkyLlD, NzLlD, 1) = 1.;
        //kXIn = kXOut;     
		solve(FFT_X, FFT_BACKWARD,  NkyLD * NzLD);
        //std::cout << " rXOut : " <<  rXOut(RxLD, NkyLlD, NzLlD, 1) << std::endl;
		Norm_X_Backward = real(rXOut(NxLlD, NkyLlD, NzLlD, 1));
		
        //std::cout << "Normalization is : " << Norm_X_Forward << " " << Norm_X_Backward << " Full : " << Norm_X << std::endl;     	
		//std::cout << "Normalization is : " << Norm_Y_Forward << " " << Norm_Y_Backward << " Full : " << Norm_Y << std::endl;     
        // broadcast normalization to all nodes
        parallel->barrier();
        parallel->send(Norm_Y_Forward, parallel->Coord(DIR_ALL) == 0); parallel->send(Norm_Y_Backward, parallel->Coord(DIR_ALL) == 0);
        parallel->send(Norm_X_Forward, parallel->Coord(DIR_ALL) == 0); parallel->send(Norm_X_Backward, parallel->Coord(DIR_ALL) == 0);
        parallel->barrier();
        //std::cout << parallel->myRank << "  Normalization is X : " << Norm_X_Forward << " " << Norm_X_Backward << " Full : " << Norm_X << std::endl;     	
		//std::cout << parallel->myRank << "  Normalization is Y : " << Norm_Y_Forward << " " << Norm_Y_Backward << " Full : " << Norm_Y << std::endl;     	
        parallel->barrier();

		//check(-1, DMESG("Normalization"));


};




