/*
 * =====================================================================================
 *
 *       Filename:  m_FFTSolver_fftw3_fftw.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/29/2010 12:59:47 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */


#include "FFTSolver_fftw3.h"
#include "Plasma.h"

using namespace Helios_fftw3;

FFTSolver_fftw3::FFTSolver_fftw3(Setup *setup, Parallel *parallel, Geometry<HELIOS_GEOMETRY> *geo) : FFTSolver(setup, parallel, geo, Nx*Nky*Nz, Nx*Nky, Nx, Nky) {
/* 
     // NOTE : We choose fortran ordering 
    //GeneralArrayStorage<3> storage_c; storage_c.ordering() = firstDim, secondDim, thirdDim;
    GeneralArrayStorage<4> storage_c; storage_c.ordering() = firstDim, secondDim, thirdDim, fourthDim;
//    --GeneralArrayStorage<3> storage_c; storage_c.ordering() = thirdDim, secondDim, firstDim;
    
    
    ///////////////////        Set Fourier space domain for 3D transform
    K3zLlD = 0; K3zLuD = Nz-1;
    K3yLlD = 0; K3yLuD = Ny-1;
    K3xLlD = 0; K3xLuD = Nx/2;
    //dir3DHC = DIR_X, 
    
    k3NxL =  K3xLuD - K3xLlD + 1;
    k3NyL =  K3yLuD - K3yLlD + 1;
    k3NzL =  K3zLuD - K3zLlD + 1;
    Rk3xL.setRange(K3xLlD, K3xLuD);
    Rk3yL.setRange(K3yLlD, K3yLuD);
    Rk3zL.setRange(K3zLlD, K3zLuD);

    ///////////////////      Set Fourier space transformation for 2DMany transform 
    K2xLlD = 0; K2xLuD = Nx/2;
    K2yLlD = 0; K2yLuD = Ny-1;
    
    Range Rk2xL; Rk2xL.setRange(K2xLlD, K2xLuD);
    Range Rk2yL; Rk2yL.setRange(K2yLlD, K2yLuD);

    //////////////////      Set Fourier space transformation for 1D transform (not implemented) but possibly
    //                      needed for fluxAvrg. (see Imadera)



#ifdef PARALLEL_OPENMP
    fftw_init_threads();
    fftw_plan_with_nthreads(parallel->numThreads);
#endif



     Array4d r3In_t (RxLD, RyLD, RzLD, Range(1, plasma->nfields), storage_c); r3In .reference(r3In_t );
     Array4d r3Out_t(RxLD, RyLD, RzLD, Range(1, plasma->nfields), storage_c); r3Out.reference(r3Out_t);


    // Configure 2DMany transform
     
     Array4z k2In_t  (Rk2xL, Rk2yL, RzLD, Range(1, plasma->nfields), storage_c); k2In .reference(k2In_t ); 
     Array4z k2Out_t (Rk2xL, Rk2yL, RzLD, Range(1, plasma->nfields), storage_c); k2Out.reference(k2Out_t);
    
     int size[2]   = { Nx, Ny };
     int k_size[2] = { Nx, Ny/2+1};

     plan_Phi2DManyForward  =fftw_plan_many_dft_r2c (2, size, Nz, r3In.data()                   ,NULL, 1, size[0] * size[1], (fftw_complex *) k2Out.data(), NULL, 1, k_size[0]* k_size[1], FFTW_ESTIMATE);  
     plan_Phi2DManyBackward =fftw_plan_many_dft_c2r (2, size, Nz, (fftw_complex * ) k2In.data() ,NULL, 1, k_size[0] * k_size[1], r3Out.data()             , NULL, 1, size[0]  * size[1]  , FFTW_ESTIMATE);  



     if(flags && FFT_X) {
          rXIn.resize(RxLD);
          kXOut.resize(Rk3xL);

          plan_Phi1DAvrgForward  = fftw_plan_dft_1d    (Nx, (fftw_complex *) rXIn.data(), (fftw_complex *) kXOut.data(), FFTW_FORWARD , FFTW_ESTIMATE);
          plan_Phi1DAvrgBackward = fftw_plan_dft_1d    (Nx, (fftw_complex *) kXOut.data(), (fftw_complex *) rXIn.data(), FFTW_BACKWARD, FFTW_ESTIMATE);
    }
    
     if(flags & FFT_XYZ) {
        // Configure 3D-Transform
                K3xLlD = 0;                     K3xLuD = Nx-1;
                K3yLlD = 0;    K3yLuD = Ny/2;
                K3zLlD = 0;                     K3zLuD = Nz/2;
     
                Rk3xL.setRange(K3xLlD, K3xLuD); Rk3yL.setRange(K3yLlD, K3yLuD); Rk3zL.setRange(K3zLlD, K3zLuD);

        GeneralArrayStorage<3> storage_3r; storage_3r.ordering() = thirdDim,  secondDim, firstDim; storage_3r.base() = NxLlD, NyLlD, NzLlD;
        Array3d r33In_t (RxLD, RyLD, RzLD, storage_3r); r33In .reference(r33In_t );
        Array3d r33Out_t(RxLD, RyLD, RzLD, storage_3r); r33Out.reference(r33Out_t);

        GeneralArrayStorage<3> storage_3c; storage_3c.ordering() = thirdDim,  firstDim, secondDim; storage_3c.base() = K3xLlD, K3yLlD, K3zLlD;
        Array3z k3Out_t  (Rk3xL, Rk3yL, Rk3zL, storage_3c); k3Out.reference(k3Out_t ); 
        Array3z k3In_t   (Rk3xL, Rk3yL, Rk3zL, storage_3c); k3In .reference(k3In_t ); 
        
        plan_Phi3DForward      = fftw_plan_dft_r2c_3d(Nz, Ny, Nx, r3In.data(), (fftw_complex *) k3Out.data(), FFTW_ESTIMATE);
        plan_Phi3DBackward     = fftw_plan_dft_c2r_3d(Nz, Ny, Nx, (fftw_complex *) k3In.data(), r3Out.data(), FFTW_ESTIMATE);
     
     }
 * */ 

}


FFTSolver_fftw3::~FFTSolver_fftw3() {

    if(flags && FFT_X) {
        fftw_destroy_plan(plan_Phi1DAvrgForward);
	    fftw_destroy_plan(plan_Phi1DAvrgBackward);
    }
	
    else if(flags && FFT_XYZ) {
        fftw_destroy_plan(plan_Phi3DForward);
	    fftw_destroy_plan(plan_Phi3DBackward);
    }

    else if(flags && FFT_XY) {
    fftw_destroy_plan(plan_Phi2DManyForward);
	fftw_destroy_plan(plan_Phi2DManyBackward);
    }
    #ifdef PARALLEL_OPENMP
     //   fftw_cleanup_threads();
    #endif
  
}

int FFTSolver_fftw3::solve(const int FFTtype, const int direction, const int nstacked) {

        if(FFTtype & FFT_XYZ & flags) {
            if     (direction == FFT_FORWARD ) fftw_execute(plan_Phi3DForward  );
            else if(direction == FFT_BACKWARD) fftw_execute(plan_Phi3DBackward);
            else   check(-1, DMESG("No such FFT direction"));
        } 
        else if(FFTtype & FFT_XY & flags) {
            if      (direction == FFT_FORWARD ) fftw_execute(plan_Phi2DManyForward);
            else if (direction == FFT_BACKWARD) fftw_execute(plan_Phi2DManyBackward);
            else    check(-1, DMESG("No such FFT direction"));
        }
        else if(FFTtype & FFT_X & flags) {
            if      (direction == FFT_FORWARD ) fftw_execute(plan_Phi1DAvrgForward);
            else if (direction == FFT_BACKWARD) fftw_execute(plan_Phi1DAvrgBackward);
            else    check(-1, DMESG("No such FFT direction"));
        } 
        else  check(-1, DMESG("Unknown FFT type or not supported"));

       return HELIOS_SUCCESS;
}


std::string FFTSolver_fftw3::getLibraryName() {
        return std::string(fftw_version);
}
  
int FFTSolver_fftw3::getDecomposition() {
      return DECOMP_NO;
}

