/*
 * =====================================================================================
 *
 *       Filename:  m_FFTSolver_p3dfft.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/29/2010 01:09:30 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */


#include "FFTSolver_p3dfft.h"
#include "p3dfft.h"
#include "Plasma.h"


FFTSolver_p3dfft::FFTSolver_p3dfft(Setup *setup, Parallel *parallel) : FFTSolver(setup, parallel, Nx*Ny*Nz, Nx*Ny, Nx) {
  
    if(parallel->Coord(DIR_FFT) == 0) {

        // p3dFFT initialization routines, check if they are equal to ours ...
        int dim[2] = { parallel->decomposition(DIR_X), parallel->decomposition(DIR_Y)};
        p3dfft_setup(dim, Nz, Nx, Ny, true, MPI_Comm_c2f(parallel->Comm[DIR_XYZ]));
        // check if domain decoposition of p3dfft is same as ours 
        int p3d_size[3];
        int conf=1;
        int k3_start[3], k3_end[3]; 
        p3dfft_get_dims(k3_start,k3_end,p3d_size,conf);

        ///////////////////////////////////////////////////////////////////////////////////
        // p3dfft assumes Fortran style arrays, so substract 1
        // check if domains matches
        int NxLlD_fft = k3_start[1] -1 + 3; check((NxLlD_fft == NxLlD) ? 0 : -1, DMESG("FFTSolver Error : Domains are different"));
        int NxLuD_fft = k3_end[1]   -1 + 3; check((NxLuD_fft == NxLuD) ? 0 : -1, DMESG("FFTSolver Error : Domains are different"));
        int NyLlD_fft = k3_start[2] -1 + 3; check((NyLlD_fft == NyLlD) ? 0 : -1, DMESG("FFTSolver Error : Domains are different"));
        int NyLuD_fft = k3_end[2]   -1 + 3; check((NyLuD_fft == NyLuD) ? 0 : -1, DMESG("FFTSolver Error : Domains are different"));
       
      


        // OK, output array for k-space
        conf=2;
        p3dfft_get_dims(k3_start,k3_end,p3d_size,conf);

        // p3dfft assumes Fortran style arrays and X,Y transposed
        K3xLlD = k3_start[1] - 1; K3xLuD = k3_end[1] - 1;
        K3yLlD = k3_start[2] - 1; K3yLuD = k3_end[2] - 1;
        K3zLlD = k3_start[0] - 1; K3zLuD = k3_end[0] - 1;
 
       // kNxL =  K3xLuD - K3xLlD + 1;
      //  kNyL =  K3yLuD - K3yLlD + 1;
       /// kNzL =  K3zLuD - K3zLlD + 1;
       // RkxL.setRange(K3xLlD, K3xLuD);
       // RkyL.setRange(K3yLlD, K3yLuD);
        Rk3zL.setRange(K3zLlD, K3zLuD);
        Rk3xL.setRange(K3xLlD, K3xLuD);
        Rk3yL.setRange(K3yLlD, K3yLuD);
     
        }
     	    
        GeneralArrayStorage<4> storage_3; storage_3.ordering() = fourthDim, thirdDim, firstDim, secondDim;
        Array4d r3In_t  (RxLD, RyLD, RzLD, Range(1, plasma->nfields), storage_3); r3In .reference(r3In_t );
        Array4d r3Out_t (RxLD, RyLD, RzLD, Range(1, plasma->nfields), storage_3); r3Out.reference(r3Out_t);
        // allocate arrays
        //GeneralArrayStorage<3> storage; storage.base() = NxLlD, NyLlD, NzLlD; storage.ordering() = thirdDim, firstDim, secondDim;
        GeneralArrayStorage<3> storage; storage.ordering() = thirdDim, firstDim, secondDim;
        Array3d r33In_t  (RxLD, RyLD, RzLD, storage); r33In .reference(r33In_t );
        Array3d r33Out_t (RxLD, RyLD, RzLD, storage); r33Out.reference(r33Out_t);

        //GeneralArrayStorage<3> storage_k; storage_k.base() = K3xLlD, K3yLlD, K3zLlD; storage_k.ordering() = thirdDim, firstDim, secondDim;
        GeneralArrayStorage<3> storage_k; storage_k.ordering() = thirdDim, firstDim, secondDim;
        Array3c k3In_t  (Rk3xL, Rk3yL, Rk3zL,  storage_k);   k3In.reference(k3In_t);
        Array3c k3Out_t (Rk3xL, Rk3yL, Rk3zL,  storage_k);   k3Out.reference(k3Out_t);

     if(flags & FFT_X) {
       ptrdiff_t num_elements, Nk, local_ni, local_istart, Kstart;
       num_elements = fftw_mpi_local_size_1d(Nx, parallel->Comm[DIR_X], FFTW_FORWARD, 0, &local_ni, &local_istart, &Nk, &Kstart);
       check((local_istart + 1 + 2 == RxLD.first()) && (local_ni == RxLD.last() - RxLD.first() + 1), DMESG("FFTW domain for 1D & 3D FFT does not match !"));
       // BUG really + 1 ?
       rXIn .resize(Range(NxLlD, NxLlD + num_elements + 1)); kXOut.resize(Range(K3xLlD, K3xLlD + num_elements + 1));
       rXOut.resize(Range(NxLlD, NxLlD + num_elements + 1)); kXIn .resize(Range(K3xLlD, K3xLlD + num_elements + 1));
      
       if((rXIn.numElements() < num_elements) || (kXOut.numElements() < num_elements)) check(-1, DMESG("r3In/k3Out Array too small"));
      
        plan_Phi1DAvrgForward  = fftw_mpi_plan_dft_1d(Nx, (fftw_complex *) rXIn.data(), (fftw_complex *) kXOut.data(), parallel->Comm[DIR_X], FFTW_FORWARD , FFTW_ESTIMATE);
        plan_Phi1DAvrgBackward = fftw_mpi_plan_dft_1d(Nx, (fftw_complex *) kXIn.data(), (fftw_complex *) rXOut.data(), parallel->Comm[DIR_X], FFTW_BACKWARD , FFTW_ESTIMATE);
    }
  // why do we have to broadcast these values ?
  //MPI_Bcast(&NxLlD, 1, MPI_INT, parallel->rankFFTMaster, parallel->CommSlaves);
  //MPI_Bcast(&NxLuD, 1, MPI_INT, parallel->rankFFTMaster, parallel->CommSlaves);
  //MPI_Bcast(&NyLlD, 1, MPI_INT, parallel->rankFFTMaster, parallel->CommSlaves);
  //MPI_Bcast(&NyLuD, 1, MPI_INT, parallel->rankFFTMaster, parallel->CommSlaves);


}


FFTSolver_p3dfft::~FFTSolver_p3dfft()
{
        p3dfft_clean();
    if(flags & FFT_X) {
        fftw_destroy_plan(plan_Phi1DAvrgForward);
	    fftw_destroy_plan(plan_Phi1DAvrgBackward);
    }
}

int FFTSolver_p3dfft::solve(const int FFTtype, const int direction, const int nstacked) 
{
        if(FFTtype & FFT_3D & flags) {
            if     (direction == FFT_FORWARD ) p3dfft_ftran_r2c(r33In.data(), (double *) k3Out.data());     
            else if(direction == FFT_BACKWARD) p3dfft_btran_c2r((double *) k3In.data(), r33Out.data());     
            else   check(-1, DMESG("No such FFT direction"));
        } 
        else if(FFTtype & FFT_2D & flags) {
        //    if     (direction == FFT_FORWARD ) p3dfft_ftran_r2c(r33In.data(), (double *) k3Out.data());     
        //    else if(direction == FFT_BACKWARD) p3dfft_btran_c2r((double *) k3In.data(), r33Out.data());     
//            else
   check(-1, DMESG("No such FFT direction"));
        } 
        else if(FFTtype & FFT_X & flags) {

            if     (direction == FFT_FORWARD)      fftw_execute(plan_Phi1DAvrgForward);
            else if(direction == FFT_BACKWARD)     fftw_execute(plan_Phi1DAvrgBackward);
            else   check(-1, DMESG("No such FFT direction"));
    }
    else  check(-1, DMESG("Unknown FFT type or not supported"));
    return HELIOS_FAILED;
    
}



std::string FFTSolver_p3dfft::getLibraryName()
{
        return std::string("P3DFFT 2.2.3");
};
  
int FFTSolver_p3dfft::getDecomposition() {
 
      return DECOMP_XY;
}

