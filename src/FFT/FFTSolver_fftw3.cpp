/*
 * =====================================================================================
 *
 *       Filename:  FFTSolver_fftw3.cpp
 *    Description:  
 *        Version:  0.1
 *         Author:  Paul P. Hilscher (2009, 2011), 
 *
 * =====================================================================================
 */


#include "FFTSolver_fftw3.h"
#include "Plasma.h"

#include <fftw3-mpi.h>
  
fftw_plan plan_YForward, plan_YBackward;
fftw_plan plan_XForward, plan_XBackward;

fftw_plan plan_AA_YForward, plan_AA_YBackward;

FFTSolver_fftw3::FFTSolver_fftw3(Setup *setup, Parallel *parallel, Geometry<HELIOS_GEOMETRY> *geo) : FFTSolver(setup, parallel, geo, Nx*Nky*Nz, Nx*(2*Nky-2), Nx, 2*Nky-2) {
   

   if(parallel->Coord(DIR_FFT) == 0) {

      int perf_flag = (setup->get("FFTW3.Measure", 1) == 1) ? FFTW_MEASURE :  FFTW_ESTIMATE; 

#ifdef PARALLEL_OPENMP
      fftw_init_threads();
      fftw_plan_with_nthreads(parallel->numThreads);
#endif

#ifdef HELIOS_PARALLEL_MPI
      fftw_mpi_init();
#endif
           

      // needs for poissons equation
      if(flags & FFT_X) {

         // set and check bounds 
         long X_NxLD, X_NxLlD, X_NkxL, X_NkxLlD, X_numElements, X_Nx = Nx; 
         X_numElements = fftw_mpi_local_size_1d(Nx, parallel->Comm[DIR_X], FFTW_FORWARD, 0, &X_NxLD, &X_NxLlD, &X_NkxL, &X_NkxLlD);
        
         // allocate arrays
         data_X_kIn  = (cmplxd *) fftw_alloc_complex(X_numElements*NyLD*NzLD*plasma->nfields);
         data_X_kOut = (cmplxd *) fftw_alloc_complex(X_numElements*NyLD*NzLD*plasma->nfields);
         
         data_X_rOut = (cmplxd *) fftw_alloc_complex(X_numElements*NyLD*NzLD*plasma->nfields);
         data_X_rIn  = (cmplxd *) fftw_alloc_complex(X_numElements*NyLD*NzLD*plasma->nfields);
               
         // set and check bounds 
         K1xLlD = X_NkxLlD;       K1xLuD = X_NkxLlD + X_NkxL - 1;
         Rk1xL.setRange(K1xLlD, K1xLuD);
                
         // allocate arrays NOTE : 1D MPI fftw-2 are ALWAYS IN-PLACE TRANSFORMS (but not fftw-3 !)
         GeneralArrayStorage<4> storage_rX; storage_rX.ordering() = fourthDim, thirdDim,  secondDim, firstDim; storage_rX.base() = NxLlD , NkyLlD, NzLlD, 1;
         
         Array4z rXIn_t   ((cmplxd *) data_X_rIn, shape(X_NxLD , NkyLD, NzLD, plasma->nfields), neverDeleteData, storage_rX); rXIn  .reference(rXIn_t ) ; rXIn  = 0.; 
         Array4z rXOut_t  ((cmplxd *) data_X_rOut, shape(X_NxLD , NkyLD, NzLD, plasma->nfields), neverDeleteData, storage_rX); rXOut .reference(rXOut_t) ; rXOut = 0.;
	        
         GeneralArrayStorage<4> storage_cX; storage_cX.ordering() = fourthDim, thirdDim,  secondDim, firstDim; storage_cX.base() = X_NkxLlD, NkyLlD, NzLlD, 1;
         
         Array4z kXOut_t  ((cmplxd *) data_X_kOut, shape(X_NkxL, NkyLD, NzLD, plasma->nfields), neverDeleteData, storage_cX); kXOut.reference(kXOut_t ) ; kXOut = 0.; 
         Array4z kXIn_t   ((cmplxd *) data_X_kIn , shape(X_NkxL, NkyLD, NzLD, plasma->nfields), neverDeleteData, storage_cX); kXIn .reference(kXIn_t  ) ; kXIn  = 0.;
    
         plan_XForward   = fftw_mpi_plan_many_dft(1, &X_Nx, (long) NkyLD,  NxLD, X_NkxL, (fftw_complex *) rXIn.data(), (fftw_complex *) kXOut.data(), parallel->Comm[DIR_X], FFT_FORWARD, perf_flag);
         plan_XBackward  = fftw_mpi_plan_many_dft(1, &X_Nx, (long) NkyLD,  NxLD, X_NkxL, (fftw_complex *) kXIn.data(), (fftw_complex *) rXOut.data(), parallel->Comm[DIR_X], FFT_BACKWARD, perf_flag);
  
      }
           
            
      // Needed to calculate non-linearity in real space
      if(flags & FFT_Y) {
                // array is too big (but for safety, shoukd be 3 but we set 4, maybe there are some roudning errors for AA?!)
                data_Y_kIn  = (cmplxd *) fftw_alloc_complex(4*(NyLD/2+1)*NxLD*NzLD*plasma->nfields);
                data_Y_kOut = (cmplxd *) fftw_alloc_complex(4*(NyLD/2+1)*NxLD*NzLD*plasma->nfields);
                data_Y_rOut = (double *) fftw_alloc_real   (4*NyLD*NxLD*NzLD*plasma->nfields);
                data_Y_rIn  = (double *) fftw_alloc_real   (4*NyLD*NxLD*NzLD*plasma->nfields);
                
                // allocate arrays NOTE : 1D MPI fftw-2 are ALWAYS IN-PLACE TRANSFORMS
                GeneralArrayStorage<4> storage_rY; storage_rY.ordering() = fourthDim, thirdDim,  secondDim, firstDim; storage_rY.base() = NxLlD, NyLlD, NzLlD, 1;
                Array4d rYIn_t   ( data_Y_rIn,  shape(NxLD, NyLD , NzLD, plasma->nfields), neverDeleteData, storage_rY); rYIn .reference(rYIn_t ) ; rYIn  = 0.; 
                Array4d rYOut_t  ( data_Y_rOut, shape(NxLD, NyLD , NzLD, plasma->nfields), neverDeleteData, storage_rY); rYOut.reference(rYOut_t) ; rYOut = 0.;
                
                GeneralArrayStorage<4> storage_cY; storage_cY.ordering() = fourthDim, thirdDim,  secondDim, firstDim; storage_cY.base() = NxLlD, NkyLlD, NzLlD, 1;
                Array4z kYOut_t  (data_Y_kOut, shape(NxLD, NkyLD, NzLD, plasma->nfields), neverDeleteData, storage_cY); kYOut.reference(kYOut_t ) ; kYOut = 0.; 
                Array4z kYIn_t   (data_Y_kIn , shape(NxLD, NkyLD, NzLD, plasma->nfields), neverDeleteData, storage_cY); kYIn .reference(kYIn_t  ) ; kYIn  = 0.;
                   
                //                                                  howmany                                  stride distance,
                plan_YForward  = fftw_plan_many_dft_r2c(1, &NyLD, NxLD*NzLD,                 rYIn.data(), NULL, 1, NyLD ,  (fftw_complex*) kYOut.data(), NULL, 1, NkyLD, perf_flag);
                plan_YBackward = fftw_plan_many_dft_c2r(1, &NyLD, NxLD*NzLD, (fftw_complex*) kYIn.data(), NULL, 1, NkyLD,                  rYOut.data(), NULL, 1, NyLD , perf_flag);


                ////////////////////////   Define Anti-Aliased Arrays /////////////////////////////////////
                AA_NkyLD  = 3*Nky/2 ; AA_NyLD   = 2*AA_NkyLD-2;
                AA_NyLlD  = NyLlD   ; AA_NyLuD  = AA_NyLlD  + AA_NyLD - 1;
                AA_NkyLlD = NkyLlD  ; AA_NkyLuD = AA_NkyLlD + AA_NkyLD - 1;

                Array4d AA_rYIn_t   ( data_Y_rIn , shape(NxLD, AA_NyLD , NzLD, plasma->nfields), neverDeleteData, storage_rY); AA_rYIn .reference(AA_rYIn_t ) ; AA_rYIn  = 0.; 
                Array4d AA_rYOut_t  ( data_Y_rOut, shape(NxLD, AA_NyLD , NzLD, plasma->nfields), neverDeleteData, storage_rY); AA_rYOut.reference(AA_rYOut_t) ; AA_rYOut = 0.;
                
                Array4z AA_kYOut_t  (data_Y_kOut , shape(NxLD, AA_NkyLD, NzLD, plasma->nfields), neverDeleteData, storage_cY); AA_kYOut.reference(AA_kYOut_t ) ; AA_kYOut = 0.; 
                Array4z AA_kYIn_t   (data_Y_kIn  , shape(NxLD, AA_NkyLD, NzLD, plasma->nfields), neverDeleteData, storage_cY); AA_kYIn .reference(AA_kYIn_t  ) ; AA_kYIn  = 0.;
                   
                plan_AA_YForward  = fftw_plan_many_dft_r2c(1, &AA_NyLD, NxLD*NzLD,                 AA_rYIn.data(), NULL, 1, AA_NyLD ,  (fftw_complex*) AA_kYOut.data(), NULL, 1, AA_NkyLD, perf_flag);
                plan_AA_YBackward = fftw_plan_many_dft_c2r(1, &AA_NyLD, NxLD*NzLD, (fftw_complex*) AA_kYIn.data(), NULL, 1, AA_NkyLD,                  AA_rYOut.data(), NULL, 1, AA_NyLD , perf_flag);
     	      

	}

     if(flags & FFT_XY ) check(-1, DMESG("XY does not work for fftw3"));
     if(flags & FFT_XYZ) check(-1, DMESG("3D Fourier transformatio not tested yet"));
   

   checkNormalization();
     
   }

}



int FFTSolver_fftw3::solve(const int FFTtype, const int direction, const int N) {
        
        if     (FFTtype & FFT_XYZ & flags) check(-1, DMESG("No such FFT direction (did you set : FFTSolver.3D = 1 )"));
        else if(FFTtype & FFT_XY  & flags) check(-1, DMESG("No such FFT direction (did you set : FFTSolver.3D = 1 )"));

        else if(FFTtype & FFT_X   & flags) {
             if     (direction == FFT_FORWARD )  fftw_execute(plan_XForward ); 
             else if(direction == FFT_BACKWARD)  fftw_execute(plan_XBackward); 
             else   check(-1, DMESG("No such FFT direction"));
        }
        else if(FFTtype & FFT_Y  & flags) {
            
             if     (direction == FFT_FORWARD )  fftw_execute(plan_YForward);
             else if(direction == FFT_BACKWARD)  fftw_execute(plan_YBackward);
             else   check(-1, DMESG("No such FFT direction"));
        }
        else  check(-1, DMESG("Unknown FFT type or not supported"));

       return HELIOS_SUCCESS;
}

int FFTSolver_fftw3::getDecomposition() {
 
      return DECOMP_X;
}



std::string FFTSolver_fftw3::getLibraryName() {
        return std::string(fftw_version) + std::string("-mpi");
}
  

FFTSolver_fftw3::~FFTSolver_fftw3() {

    //  release fftw-3  
    if(flags & FFT_X) {
        fftw_destroy_plan(plan_XForward);
	    fftw_destroy_plan(plan_XBackward);
    }
    if(flags & FFT_Y) {
        fftw_destroy_plan(plan_YForward);
	    fftw_destroy_plan(plan_YBackward);
        
        fftw_destroy_plan(plan_AA_YForward);
	    fftw_destroy_plan(plan_AA_YBackward);

    }

#ifdef PARALLEL_OPENMP
    fftw_cleanup_threads();
#endif
	
    // clean up arrays 
    fftw_free(data_Y_kIn);
    fftw_free(data_Y_kOut);
    fftw_free(data_Y_rIn);
    fftw_free(data_Y_rOut);
    
    fftw_free(data_X_kIn);
    fftw_free(data_X_rOut);
    fftw_free(data_X_rIn);
    fftw_free(data_X_kOut);

}

       
// Note : We need to take care of aliasing
Array3z FFTSolver_fftw3::multiply(Array3z &A, Array3z &B, Array3z  &R) 
{
   
   // Array A : Copy Values to larger AA-Array (set larger values to zero)
   AA_kYIn = 0.;
   for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) AA_kYIn(RxLD, y_k, RzLD, 1) = A(RxLD, y_k, RzLD); 
   fftw_execute(plan_AA_YBackward);


   AA_rYIn(RxLD, AA_RyLD, RzLD, 1) = AA_rYOut(RxLD, AA_RyLD, RzLD, 1);
   
   // Array B : Copy Values to larger AA-Array (set larger values to zero)
   AA_kYIn = 0.;
   for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) AA_kYIn(RxLD, y_k, RzLD, 1) = B(RxLD, y_k, RzLD); 
   fftw_execute(plan_AA_YBackward);


   //////////////////////// Real Space (multiply values) /////////////////// 
   
   for(int z=NzLlD; z<=NzLuD;z++) { for(int y=AA_NyLlD; y<=AA_NyLuD;y++) { for(int x= NxLlD; x <= NxLuD; x++) {
    AA_rYIn(x,y,z,1) *= AA_rYOut(x,y,z,1) / (Norm_Y * Norm_Y);
   }}}
   
   //////////////////////// End Real Space (multiply values) /////////////////// 
   
   fftw_execute(plan_YForward);
   
   // Array R : Copy Result back to array
   for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) R(RxLD, y_k, RzLD) = AA_kYOut(RxLD, y_k, RzLD); 
   return R;

};

    
void FFTSolver_fftw3::printOn(ostream &output) const {
         FFTSolver::printOn(output);
         output   << "FFTSolver  |  using fftw-3 interface for (" << std::string(fftw_version) << ")" << std::endl;
         
        }
