/*
 * =====================================================================================
 *
 *       Filename:  FFTSolver_fftw2_mpi.cpp
 *
 *    Description:  Wrapper for fftw2-mpi library (www.fftw.org) 
 *
 *        Version:  0.1
 *         Author:  Paul P. Hilscher (2009) 
 *
 * =====================================================================================
 */

#include "FFTSolver_fftw2_mpi.h"
#include "Plasma.h"

#include <rfftw_mpi.h>

// Need to define it here, otherwise it we have name collisions with fftw3.
// Namespace wrapping won't help becauase MACROS do not follow namespace convention.
rfftwnd_mpi_plan  plan_XYZForward;
rfftwnd_mpi_plan  plan_XYZBackward;

rfftwnd_mpi_plan  plan_XYForward;
rfftwnd_mpi_plan  plan_XYBackward;
  
fftw_mpi_plan plan_XForward;
fftw_mpi_plan plan_XBackward;

fftw_mpi_plan plan_YForward;
fftw_mpi_plan plan_YBackward;
            


FFTSolver_fftw2_mpi::FFTSolver_fftw2_mpi(Setup *setup, Parallel *parallel, Geometry<GKC_GEOMETRY> *geo) : FFTSolver(setup, parallel, geo, Nx*Ny*Nz, Nx*Ny, Nx, Ny) {

    if(parallel->Coord(DIR_FFT) == 0) {

            perf_flag = (setup->get("FFTW2.Measure", 0) == 0) ? FFTW_ESTIMATE : FFTW_MEASURE;


#ifdef PARALLEL_OPENMP
   //         fftw_init_threads();
   //         fftw_plan_with_nthreads(parallel->numThreads);
#endif
            // Initialze 2D FFT

            plan_XYForward   = rfftw2d_mpi_create_plan(parallel->Comm[DIR_X], Nx, Ny, FFTW_REAL_TO_COMPLEX , perf_flag );
            plan_XYBackward  = rfftw2d_mpi_create_plan(parallel->Comm[DIR_X], Nx, Ny, FFTW_COMPLEX_TO_REAL , perf_flag );
            
            int fft_NxL, fft_NxL_start, fft_NkyL, fft_NkyL_start, numElem; 
            rfftwnd_mpi_local_sizes(plan_XYForward, &fft_NxL, &fft_NxL_start, &fft_NkyL, &fft_NkyL_start, &numElem);

            // Check if FFT domain and grid domains are equal (add ghost cells)
            int fft_NxLlD =  fft_NxL_start + 3;     int fft_NxLuD = fft_NxL_start + fft_NxL - 1  + 3;
            int fft_NyLlD =  3               ;      int fft_NyLuD = Ny - 1  + 3;
      
            check(((NxLlD == fft_NxLlD) && (NxLuD == fft_NxLuD)) ? 1 : -1, DMESG("FFT Decomposition and Grid decompostion is not unique for X")); 
            check(((NyLlD == fft_NyLlD) && (NyLuD == fft_NyLuD)) ? 1 : -1, DMESG("FFT Decomposition and Grid decompostion is not unique for Y")); 
            
            data_1 = new cmplxd[numElem*NzLD*plasma->nfields];
            data_2 = new cmplxd[numElem*NzLD*plasma->nfields];
            data_3 = new cmplxd[numElem*NzLD*plasma->nfields*4];
            data_8 = new cmplxd[numElem*NzLD*plasma->nfields];
            data_9 = new cmplxd[numElem*NzLD*plasma->nfields];

            // setup storage space, noth that 
     	    GeneralArrayStorage<4> storage_r; storage_r.ordering() = fourthDim, thirdDim,  secondDim, firstDim; storage_r.base() = NxLlD, NkyLlD, NzLlD, 1;
            Array4z r3In_t  (data_1, shape(NxLD, NkyLD, NzLD, plasma->nfields), neverDeleteData, storage_r); r3In  .reference(r3In_t ) ; r3In  = 0.;
            Array4z r3Out_t (data_2, shape(NxLD, NkyLD, NzLD, plasma->nfields), neverDeleteData, storage_r); r3Out .reference(r3Out_t ); r3Out = 0.;
            
            // 2DMany - FFT ( is it necessary ?)
            K2xLlD = 0             ; K2xLuD = Nx-1;
            K2yLlD = fft_NkyL_start; K2yLuD = fft_NkyL_start + fft_NkyL - 1;
           
            Rk2xL.setRange(K2xLlD, K2xLuD);  Rk2yL.setRange(K2yLlD, K2yLuD);
           

            GeneralArrayStorage<4> storage_c; storage_c.ordering() = fourthDim, thirdDim,  firstDim, secondDim; storage_c.base() = K2xLlD, K2yLlD, NzLlD, 1;
            Array4z k2Out_t  ((cmplxd *) data_8, shape(Nx, fft_NkyL, NzLD, plasma->nfields), neverDeleteData, storage_c); k2Out.reference(k2Out_t ); k2Out = 0.;
            Array4z k2In_t   ((cmplxd *) data_9, shape(Nx, fft_NkyL, NzLD, plasma->nfields), neverDeleteData, storage_c); k2In .reference(k2In_t ) ; k2In  = 0.; 


            //
            data_4 = new cmplxd[numElem*4];
            data_5 = new cmplxd[numElem*4];
            data_6 = new cmplxd[numElem*4];
            data_7 = new cmplxd[numElem*4];

            

            // used for e.g. to solve Poisson equation [ A(x,k_y,z) ->A(x_k,y_k,z) ]
            if(flags & FFT_X) {
                 
#define HELIOS_PARALLEL_MPI 
                plan_XForward  = fftw_mpi_create_plan(parallel->Comm[DIR_X], Nx, FFTW_FORWARD , perf_flag );
                plan_XBackward = fftw_mpi_create_plan(parallel->Comm[DIR_X], Nx, FFTW_BACKWARD, perf_flag );
   
                // set and check bounds 
                int X_NxL, X_NxLlD, X_NkxL, X_NkxLlD, X_numElements; 
                fftw_mpi_local_sizes(plan_XForward, &X_NxL, &X_NxLlD, &X_NkxL, &X_NkxLlD, &X_numElements);
                K1xLlD = X_NkxLlD;       K1xLuD = X_NkxLlD + X_NkxL - 1;
                Rk1xL.setRange(K1xLlD, K1xLuD);
            
                data_X_kIn_rOut = new cmplxd[X_numElements*NyLD*NzLD*plasma->nfields];
                data_X_rIn_kOut = new cmplxd[X_numElements*NyLD*NzLD*plasma->nfields];
                
                // allocate arrays NOTE : 1D MPI fftw-2 are ALWAYS IN-PLACE TRANSFORMS
                GeneralArrayStorage<4> storage_rX; storage_rX.ordering() = fourthDim, thirdDim,  secondDim, firstDim; storage_rX.base() = NxLlD , NkyLlD, NzLlD, 1;
                Array4z rXIn_t   ((cmplxd *) data_X_rIn_kOut, shape(X_NxL , NkyLD, NzLD, plasma->nfields), neverDeleteData, storage_rX); rXIn  .reference(rXIn_t ) ; rXIn  = 0.; 
                Array4z rXOut_t  ((cmplxd *) data_X_kIn_rOut, shape(X_NxL , NkyLD, NzLD, plasma->nfields), neverDeleteData, storage_rX); rXOut .reference(rXOut_t) ; rXOut = 0.;
	        
                GeneralArrayStorage<4> storage_cX; storage_cX.ordering() = fourthDim, thirdDim,  secondDim, firstDim; storage_cX.base() = X_NkxLlD, NkyLlD, NzLlD, 1;
                Array4z kXOut_t  ((cmplxd *) data_X_rIn_kOut, shape(X_NkxL, NkyLD, NzLD, plasma->nfields), neverDeleteData, storage_cX); kXOut.reference(kXOut_t ) ; kXOut = 0.; 
                Array4z kXIn_t   ((cmplxd *) data_X_kIn_rOut, shape(X_NkxL, NkyLD, NzLD, plasma->nfields), neverDeleteData, storage_cX); kXIn .reference(kXIn_t  ) ; kXIn  = 0.;
#endif // HELIOS_PARALLEL_MPI     
            }
           
            // Initialize FFT for poloidal direction [ A(x,k_y,z, F) -> A(x,y,z, F) ] . Used e.g. for non-linear calculation
            // Note : y-direction is local (no spatial decomposition in y)
            if(flags & FFT_Y) {
                
              if(flags & FFT_AA) {
                check(-1, DMESG("AA not tested yet"));
                plan_YForward  = fftw_mpi_create_plan(parallel->Comm[DIR_Y], (Ny*3)/2, FFTW_FORWARD , perf_flag );
                plan_YBackward = fftw_mpi_create_plan(parallel->Comm[DIR_Y], (Ny*3)/2, FFTW_BACKWARD, perf_flag );
              } else {
                //plan_YForward  = fftw_mpi_create_plan(parallel->Comm[DIR_Y], Ny, FFTW_FORWARD , perf_flag );
                //plan_YBackward = fftw_mpi_create_plan(parallel->Comm[DIR_Y], Ny, FFTW_BACKWARD, perf_flag );
                plan_YForward  = fftw_mpi_create_plan(parallel->Comm[DIR_Y], Ny, FFTW_FORWARD , perf_flag );
                plan_YBackward = fftw_mpi_create_plan(parallel->Comm[DIR_Y], Ny, FFTW_BACKWARD, perf_flag );
              }
    
            // 
            //
            // s
                int Y_NyL, _Y_NyLlD, Y_NkyL, _Y_kyLlD, Y_numElements; 
                fftw_mpi_local_sizes(plan_YForward, &Y_NyL, &_Y_NyLlD, &Y_NkyL, &_Y_kyLlD, &Y_numElements);

                Y_kyLlD = _Y_kyLlD   ; Y_kyLuD = _Y_kyLlD     + Y_NkyL - 1; 
                Y_NyLlD = _Y_NyLlD+3 ; Y_NyLuD = _Y_NyLlD + 3 + Y_NyL  - 1; 

                Y_RkyL.setRange(Y_kyLlD, Y_kyLuD);
                Y_RyLD.setRange(Y_NyLlD, Y_NyLuD); 
   
                check((NyLlD  == Y_NyLlD) ? 0 : -1, DMESG("FFT - Y , boundaries not equal"));
                check((NyLuD  == Y_NyLuD) ? 0 : -1, DMESG("FFT - Y , boundaries not equal"));
                check((NkyLlD == Y_kyLlD) ? 0 : -1, DMESG("FFT - Y , boundaries not equal"));
                check((NkyLuD == Y_kyLuD) ? 0 : -1, DMESG("FFT - Y , boundaries not equal"));
    
                // need for multiply

                data_Y_kIn_rOut = new cmplxd[Y_numElements*NxLD*NzLD*plasma->nfields];
                data_Y_rIn_kOut = new cmplxd[Y_numElements*NxLD*NzLD*plasma->nfields];
                

                // allocate arrays NOTE : 1D MPI fftw-2 are ALWAYS IN-PLACE TRANSFORMS
                GeneralArrayStorage<4> storage_rY; storage_rY.ordering() = fourthDim, thirdDim,  firstDim, secondDim; storage_rY.base() = NxLlD, NyLlD, NzLlD, 1;
                Array4z rYIn_t   ((cmplxd *) data_Y_rIn_kOut, shape(NxLD, Y_NyL , NzLD, plasma->nfields), neverDeleteData, storage_rY); rYIn  .reference(rYIn_t ) ; rYIn  = 0.; 
                Array4z rYOut_t  ((cmplxd *) data_Y_kIn_rOut, shape(NxLD, Y_NyL , NzLD, plasma->nfields), neverDeleteData, storage_rY); rYOut .reference(rYOut_t) ; rYOut = 0.;
                
                GeneralArrayStorage<4> storage_cY; storage_cY.ordering() = fourthDim, thirdDim,  firstDim, secondDim; storage_cY.base() = NxLlD, NkyLlD, NzLlD, 1;
                Array4z kYOut_t  ((cmplxd *) data_Y_rIn_kOut, shape(NxLD, Y_NkyL, NzLD, plasma->nfields), neverDeleteData, storage_cY); kYOut.reference(kYOut_t ) ; kYOut = 0.; 
                Array4z kYIn_t   ((cmplxd *) data_Y_kIn_rOut, shape(NxLD, Y_NkyL, NzLD, plasma->nfields), neverDeleteData, storage_cY); kYIn .reference(kYIn_t  ) ; kYIn  = 0.;
     	      
            }
            
            
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // 3D - FFT ( is it necessary ?)
            if(flags & FFT_XYZ) {
		
                check(((parallel->decomposition(DIR_Y) > 1) || (parallel->decomposition(DIR_Z) > 1)) ? -1 : 0, DMESG("When 3D FFT is used no decomposition in Y and Z is supported"));    	   
 
                plan_XYZForward   = rfftw3d_mpi_create_plan(parallel->Comm[DIR_X], Nx, Ny, Nz, FFTW_REAL_TO_COMPLEX , perf_flag );
                plan_XYZBackward  = rfftw3d_mpi_create_plan(parallel->Comm[DIR_X], Nx, Ny, Nz, FFTW_COMPLEX_TO_REAL , perf_flag );
             
                rfftwnd_mpi_local_sizes(plan_XYZForward, &fft_NxL, &fft_NxL_start, &fft_NkyL, &fft_NkyL_start, &numElem);
            
		check((NxLlD == (fft_NxL_start+3)) && (NxLuD == (fft_NxL_start+3 + fft_NxL - 1)) ? 1 : -1, DMESG("FFT Decomposition and Grid decompostion is not unique for X")); 
            
//                data_4 = new cmplxd[numElem];
//                data_5 = new cmplxd[numElem];
//                data_6 = new cmplxd[numElem];
                
                K3xLlD = 0;                     K3xLuD = Nx-1;
                K3yLlD = fft_NkyL_start;        K3yLuD = fft_NkyL_start + fft_NkyL - 1;
                K3zLlD = 0;                     K3zLuD = Nz/2;

                Rk3xL.setRange(K3xLlD, K3xLuD); Rk3yL.setRange(K3yLlD, K3yLuD); Rk3zL.setRange(K3zLlD, K3zLuD);
           
                GeneralArrayStorage<3> storage_3r; storage_3r.ordering() = thirdDim,  secondDim, firstDim; storage_3r.base() = NxLlD, NyLlD, NzLlD;
                Array3z r33In_t  (data_4, shape(fft_NxL, NyLD, 2*(Nz/2 + 1)), neverDeleteData, storage_3r); r33In  .reference(r33In_t ); 
                Array3z r33Out_t (data_5, shape(fft_NxL, NyLD, 2*(Nz/2 + 1)), neverDeleteData, storage_3r); r33Out .reference(r33Out_t ); 
            
                GeneralArrayStorage<3> storage_3c; storage_3c.ordering() = thirdDim,  firstDim, secondDim; storage_3c.base() = K3xLlD, K3yLlD, K3zLlD;
                Array3z k3Out_t  ((cmplxd *) data_4, shape(K3xLuD - K3xLlD+1, K3yLuD - K3yLlD + 1, K3zLuD - K3zLlD + 1), neverDeleteData, storage_3c); k3Out.reference(k3Out_t ); 
                Array3z k3In_t   ((cmplxd *) data_5, shape(K3xLuD - K3xLlD+1, K3yLuD - K3yLlD + 1, K3zLuD - K3zLlD + 1), neverDeleteData, storage_3c); k3In .reference(k3In_t ); 
            
 	      
            }
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    

    }


}


FFTSolver_fftw2_mpi::~FFTSolver_fftw2_mpi() {

    if(flags & FFT_X) {
        fftw_mpi_destroy_plan(plan_XForward);
	    fftw_mpi_destroy_plan(plan_XBackward);
    }
	
    if(flags & FFT_XYZ) {
        rfftwnd_mpi_destroy_plan(plan_XYZForward);
      	rfftwnd_mpi_destroy_plan(plan_XYZBackward);
        
        delete[] data_4; delete[] data_5; delete[] data_6;
    }


    if(flags & FFT_XY) {
            rfftwnd_mpi_destroy_plan(plan_XYForward);
      	    rfftwnd_mpi_destroy_plan(plan_XYBackward);
    }

        delete[] data_1;
        delete[] data_2;
        delete[] data_3;
}

int FFTSolver_fftw2_mpi::solve(const int FFTtype, const int direction, const int N) {
        
        if(FFTtype & FFT_XYZ & flags) {
//             if     (direction == FFT_FORWARD )  rfftwnd_mpi(plan_XYZForward , 1, (fftw_real *) r33In.data(), data_6, FFTW_TRANSPOSED_ORDER);
//             else if(direction == FFT_BACKWARD)  rfftwnd_mpi(plan_XYZBackward, 1, (fftw_real *) k3In.data(),  data_6, FFTW_TRANSPOSED_ORDER);
  //           else   check(-1, DMESG("No such FFT direction (did you set : FFTSolver.3D = 1 )"));
        } 
	
        else if(FFTtype & FFT_XY & flags) {
//             if     (direction == FFT_FORWARD )  rfftwnd_mpi(plan_XYForward , nstacked*NzLD, (fftw_real *) r3In.data(), data_3, FFTW_TRANSPOSED_ORDER);
//             else if(direction == FFT_BACKWARD)  rfftwnd_mpi(plan_XYBackward, nstacked*NzLD, (fftw_real *) k2In.data(), data_3, FFTW_TRANSPOSED_ORDER);
//             else   check(-1, DMESG("No such FFT direction"));
        } 
        else if(FFTtype & FFT_X & flags) {
             if     (direction == FFT_FORWARD )  fftw_mpi(plan_XForward  , N, (fftw_complex *) rXIn.data(), (fftw_complex *) data_3);
             else if(direction == FFT_BACKWARD)  fftw_mpi(plan_XBackward , N, (fftw_complex *) kXIn.data(), (fftw_complex *) data_3);
             else   check(-1, DMESG("No such FFT direction"));
        }
        else if(FFTtype & FFT_Y & flags) {
             if     (direction == FFT_FORWARD )  fftw_mpi(plan_YForward  , N, (fftw_complex *) rYIn.data(), (fftw_complex *) data_3);
             else if(direction == FFT_BACKWARD)  fftw_mpi(plan_YBackward , N, (fftw_complex *) kYIn.data(), (fftw_complex *) data_3);
             else   check(-1, DMESG("No such FFT direction"));
        }
        else  check(-1, DMESG("Unknown FFT type or not supported"));

       return GKC_SUCCESS;
}

std::string FFTSolver_fftw2_mpi::getLibraryName() {
        return std::string(fftw_version) + std::string("-mpi");
}
  
int FFTSolver_fftw2_mpi::getDecomposition() {
 
      return DECOMP_X;
}

       


void FFTSolver_fftw2_mpi::copyNyquiest(Array4z f) {
  int N = (flags & FFT_AA) ? (Ny*3)/2 : Ny;
  f(RxLD, N-Ny/2, RzLD, 1) = -f(RxLD, Ny/2, RzLD, 1); 
};

int FFTSolver_fftw2_mpi::mapAAY(const int y_k) {
    int N = (flags & FFT_AA) ? (Ny*3)/2 : Ny;
    return (y_k <= Ny/2) ? y_k : N - (Ny-y_k);
};




Array3z FFTSolver_fftw2_mpi::multiply(Array3z &A, Array3z &B, Array3z  &R) {

   // back transform dphi_dx(x,ky,z,m,s) -> dphi_dx(x,y,z,m,s) 
   kYIn = 0.;
   for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) kYIn(RxLD, mapAAY(y_k), RzLD, 1) = A(RxLD, y_k, RzLD); 
   if(flags & FFT_AA) copyNyquiest(kYIn);
   solve(FFT_Y, FFT_BACKWARD, NxLD * NzLD);
  
   // save  
   AA_xyz(RxLD, Y_RyLD, RzLD) = rYOut(RxLD, Y_RyLD, RzLD, 1);

   // second Term : dphi_dy
   kYIn = 0.;
   //#-pragma omp parallel for
   for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) kYIn(RxLD, mapAAY(y_k), RzLD, 1) = B(RxLD, y_k, RzLD);
   if(flags & FFT_AA) copyNyquiest(kYIn);
   solve(FFT_Y, FFT_BACKWARD, NxLD * NzLD);


   // multiply-values 
   
   //AA_xyz(RxLD, Y_RyLD, RzLD) *= rYOut(RxLD, Y_RyLD, RzLD, 1)/Norm_Y;
   for(int z=NzLlD; z<=NzLuD;z++) { for(int y=Y_NyLlD; y<=Y_NyLuD;y++) { for(int x= NxLlD; x <= NxLuD; x++) {
    AA_xyz(x,y,z) = AA_xyz(x,y,z) * rYOut(x,y,z,1);
   }}}



   // forward transform
   rYIn(RxLD,Y_RyLD, RzLD, 1) = AA_xyz(RxLD, Y_RyLD, RzLD); 
   solve(FFT_Y, FFT_FORWARD, NxLD * NzLD);

   return kYOut(RxLD, RkyLD, RzLD, 1);

    
    
};
