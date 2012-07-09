/*
 * =====================================================================================
 *
 *       Filename:  m_FFTSolver_fftw2_fftw-mpi.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/29/2010 12:57:33 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */


#include "FFTSolver_fftw2.h"

#include "Global.h"

FFTSolver_fftw2::FFTSolver_fftw2(Setup *setup, Parallel *parallel, Geometry<HELIOS_GEOMETRY> *geo) : FFTSolver(setup, parallel, geo, Nx*Ny*Nz, Nx*Ny, Nx) {
       
   check(-1, DMESG("Not working, Please fix"));
   /* 
    int  NxGC = 2 ;  int  NyGC = 2 ; 

    if(parallel->isFFTGroup) {

            plan_Phi2DManyForward   = rfftw2d_create_plan(Nx, Ny, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
            plan_Phi2DManyBackward  = rfftw2d_create_plan(Nx, Ny, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);
//            plan_Phi2DManyForward   = rfftw2d_create_plan(Nx, Ny, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
//            plan_Phi2DManyBackward  = rfftw2d_create_plan(Nx, Ny, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);
            
            plan_Phi3DForward   = rfftw3d_create_plan(Nx, Ny, Nz, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
            plan_Phi3DBackward  = rfftw3d_create_plan(Nx, Ny, Nz, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);


            // check and allocate necessary space, due to pending work buffers the size can be somewhat larger
            
            K2xLlD = 0;       K2xLuD = Nx/2;
            K2yLlD = 0;       K2yLuD = Ny-1;
           
            // BUG, we need mult. 2, otherwise program crashes, why ? 
            data_1 = new double[(Nx+1)*Ny*Nz*2];
            data_2 = new double[(Nx+1)*Ny*Nz*2];
            data_3 = new double[(Nx+1)*Ny*Nz*2];
            data_4 = new double[(Nx+1)*Ny*Nz*2];

            // setup storage space, noth that 
     	    GeneralArrayStorage<4> storage_r; storage_r.ordering() = firstDim, secondDim, thirdDim, fourthDim; storage_r.base() = NxLlD, NyLlD, NzLlD, 1;
            Array3d r4In_t  (data_1, shape(NxLuD - NxLlD+1, NyLuD - NyLlD + 1, NzLuD - NzLlD + 1, Range(1, plasma->nFields)), neverDeleteData, storage_r); r3In  .reference(r3In_t ); 
            Array3d r4Out_t (data_2, shape(NxLuD - NxLlD+1, NyLuD - NyLlD + 1, NzLuD - NzLlD + 1, Range(1, plasma->nFields)), neverDeleteData, storage_r); r3Out .reference(r3Out_t ); 
            
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // 3D - FFT ( is it necessary ?)
            K3xLlD = 0;       K3xLuD = Nx/2;
            K3yLlD = 0;       K3yLuD = Ny-1;
            K3zLlD = 0;       K3zLuD = Nz-1;

            Rk3xL.setRange(K3xLlD, K3xLuD);
            Rk3yL.setRange(K3yLlD, K3yLuD);
            Rk3zL.setRange(K3zLlD, K3zLuD);
            
            
            GeneralArrayStorage<4> storage_3c; storage_3c.ordering() = firstDim, secondDim, thirdDim; storage_3c.base() = K3xLlD, K3yLlD, K3zLlD, 1;
            Array3c k3Out_t  ((cmplxd *) data_1, shape(K3xLuD - K3xLlD+1, K3yLuD - K3yLlD + 1, K3zLuD - K3zLlD + 1), neverDeleteData, storage_3c); k3Out.reference(k3Out_t ); 
            Array3c k3In_t   ((cmplxd *) data_2, shape(K3xLuD - K3xLlD+1, K3yLuD - K3yLlD + 1, K3zLuD - K3zLlD + 1), neverDeleteData, storage_3c); k3In .reference(k3In_t ); 

            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // 2DMany - FFT ( is it necessary ?)

            Range Rk2xL; Rk2xL.setRange(K2xLlD, K2xLuD);
            Range Rk2yL; Rk2yL.setRange(K2yLlD, K2yLuD);
	        GeneralArrayStorage<4> storage_c; storage_c.ordering() = firstDim, secondDim, thirdDim; storage_c.base() = K2xLlD, K2yLlD, NzLlD;
            Array4c k2Out_t  ((cmplxd *) data_3, shape(K2xLuD - K2xLlD + 1, K2yLuD - K2yLlD + 1, NzLD), neverDeleteData, storage_c); k2Out.reference(k2Out_t ); 
            Array4c k2In_t   ((cmplxd *) data_4, shape(K2xLuD - K2xLlD + 1, K2yLuD - K2yLlD + 1, NzLD), neverDeleteData, storage_c); k2In .reference(k2In_t ); 
           

     if(flags & FFT_X) {
  
       ptrdiff_t num_elements, Nk, local_ni, local_istart, Kstart;
       num_elements = fftw_local_size_1d(Nx, parallel->CommFFT1D, FFTW_FORWARD, 0, &local_ni, &local_istart, &Nk, &Kstart);
       check((local_istart + 1 + 2 == RxLD.first()) && (local_ni == RxLD.last() - RxLD.first() + 1), DMESG("FFTW domain for 1D & 3D FFT does not match !"));
      if((Phi1A.numElements() < num_elements) || (Phi1B.numElements() < num_elements)) check(-1, DMESG("r3In/k3Out Array too small"));
      
        plan_Phi1DAvrgForward  = fftw_plan_dft_1d(Nx, (fftw_complex *) Phi1A.data(), (fftw_complex *) Phi1B.data(), parallel->CommFFT1D, FFTW_FORWARD , FFTW_PATIENT);
        plan_Phi1DAvrgBackward = fftw_plan_dft_1d(Nx, (fftw_complex *) Phi1B.data(), (fftw_complex *) Phi1A.data(), parallel->CommFFT1D, FFTW_BACKWARD , FFTW_PATIENT);
    }
    

    }

    */ 

}


FFTSolver_fftw2::~FFTSolver_fftw2() {

/* 
    if(flags & FFT_X) {
        fftw_destroy_plan(plan_Phi1DAvrgForward);
	    fftw_destroy_plan(plan_Phi1DAvrgBackward);
    }
	
    if(flags & FFT_3D) {
        fftw_destroy_plan(plan_Phi3DForward);
	    fftw_destroy_plan(plan_Phi3DBackward);
    }
 * */

    if(flags & FFT_2D) {
        rfftwnd_destroy_plan(plan_Phi2DManyForward);
	    rfftwnd_destroy_plan(plan_Phi2DManyBackward);
    }

        delete[] data_1;
        delete[] data_2;
}

int FFTSolver_fftw2::solve(int FFTtype, int direction) {
        if(FFTtype & FFT_3D) {
          std::cout << " 3DFFT" << std::endl;
            if(direction == FFT_FORWARD)  rfftwnd_one_real_to_complex(plan_Phi3DForward , r3In.data()           , (fftw_complex *)      k3Out.data());
            if(direction == FFT_BACKWARD) rfftwnd_one_complex_to_real(plan_Phi3DBackward, (fftw_complex *) k3In.data(), r3Out.data());
            return HELIOS_SUCCESS;
	}
    //        if     (direction == FFT_FORWARD ) fftw_execute(plan_Phi3DForward  );
    //        else if(direction == FFT_BACKWARD) fftw_execute(plan_Phi3DBackward);
    //        else   check(-1, DMESG("No such FFT direction"));
   /* 
        else if(FFTtype & FFT_X & flags) {

        if(direction == FFT_FORWARD)      fftw_execute(plan_Phi1DAvrgForward);
        
        if(direction == FFT_BACKWARD)     fftw_execute(plan_Phi1DAvrgBackward);
        return HELIOS_SUCCESS;
    }
    * */ 
    else if(FFTtype & FFT_2D ) {
          std::cout << " 2DFFT" << std::endl;
          // Set working space !
            if(direction == FFT_FORWARD)  rfftwnd_real_to_complex(plan_Phi2DManyForward , Nz, r3In.data()           , 1, Nx*Ny, (fftw_complex *) k2Out.data(), 1, (Nx/2+1)*Ny);
            if(direction == FFT_BACKWARD) rfftwnd_complex_to_real(plan_Phi2DManyBackward, Nz, (fftw_complex *) k2In.data(), 1, (Nx/2+1)*Ny, r3Out.data(), 1, Nx*Ny);
            return HELIOS_SUCCESS;
        }
    else  check(-1, DMESG("Unknown FFT type or not supported"));
        
    return HELIOS_FAILED;
}

std::string FFTSolver_fftw2::getLibraryName() {
        return std::string(fftw_version) + std::string("-mpi");
}
  
int FFTSolver_fftw2::getDecomposition() {
 
      return DECOMP_X;
}

       
