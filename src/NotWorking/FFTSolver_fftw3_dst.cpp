/*
 * =====================================================================================
 *
 *       Filename:  m_FFTSolver_fftw3_dst_fftw.cpp
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


#include "FFTSolver_fftw3_dst.h"
#include "Plasma.h"
#include <complex.h>

using namespace Helios_fftw3;


FFTSolver_fftw3_dst::FFTSolver_fftw3_dst(Setup *setup, Parallel *parallel, Geometry<HELIOS_GEOMETRY> *geo) : FFTSolver(setup, parallel, geo, Nx*Nky*Nz, 2*Nx*Nky, Nx, Nky, false) 
{

/* 
#ifdef PARALLEL_OPENMP
    fftw_init_threads();
    fftw_plan_with_nthreads(parallel->numThreads);
#endif
    DTtype = HELIOS_FFT_DST;

    // NOTE : We choose fortran ordering 
    //GeneralArrayStorage<4> storage_c; storage_c.ordering() = firstDim, secondDim, thirdDim, fourthDim;
    
    
    //GeneralArrayStorage<4> storage_c; storage_c.ordering() = firstDim, secondDim, thirdDim, fourthDim;
    GeneralArrayStorage<4> storage_c; storage_c.ordering() = secondDim, firstDim, thirdDim, fourthDim;
    

    ///////////////////      Set Fourier space transformation for 2DMany transform 
    // DST (Discrete Sin Transformation in x-direction)
    
    // Note : We sum sin from 0, to N-1, thus domain is different
    Range RxLD2; RxLD2.setRange(NxLlD, NxLuD+Nx);
    Array4d r3In_t (RxLD2, RyLD, RzLD, Range(1, plasma->nfields), storage_c); r3In .reference(r3In_t );
    Array4d r3Out_t(RxLD2, RyLD, RzLD, Range(1, plasma->nfields), storage_c); r3Out.reference(r3Out_t);


    K2xLlD = 0; K2xLuD = 2*Nx-1;
    //K2xLlD = 0; K2xLuD = Nx-1;
    K2yLlD = 0; K2yLuD = Ny/2  ;
    
    Rk2xL.setRange(K2xLlD, K2xLuD);
    Rk2yL.setRange(K2yLlD, K2yLuD);
    
    Rk2xL2.setRange(K2xLlD, 2*Nx-1);
     
    // Configure 2DMany transform
    Array4z k2In_t  (Rk2xL2, Rk2yL, RzLD, Range(1, plasma->nfields), storage_c); k2In .reference(k2In_t ); 
    Array4z k2Out_t (Rk2xL2, Rk2yL, RzLD, Range(1, plasma->nfields), storage_c); k2Out.reference(k2Out_t); 
    
    
    int size[2]   = { 2*Nx, Ny     };
    int k_size[2] = { 2*Nx, Ny/2+1 };


    plan_Phi2DManyForward  =fftw_plan_many_dft_r2c (2, size, Nz, r3In.data()                   ,NULL, 1, size[0] * size[1], (fftw_complex *) k2Out.data(), NULL, 1, k_size[0]*k_size[1], FFTW_ESTIMATE);  
    plan_Phi2DManyBackward =fftw_plan_many_dft_c2r (2, size, Nz, (fftw_complex * ) k2In.data() ,NULL, 1, k_size[0] * k_size[1], r3Out.data()             , NULL, 1, size[0]*size[1]    , FFTW_ESTIMATE);  

//     const fftw_r2r_kind fftw_kind_forward[] = { FFTW_RODFT10, FFTW_R2HC };
 //    const fftw_r2r_kind fftw_kind_forward[] = { FFTW_R2HC, FFTW_R2HC };
//     plan_Phi2DManyForward  =fftw_plan_many_r2r (2, size, Nz, r3In.data()                   ,NULL, 1, size[0] * size[1], (double *) k2Out.data(), NULL, 1, k_size[0]*k_size[1], fftw_kind_forward,  FFTW_ESTIMATE);  
//     const fftw_r2r_kind fftw_kind_backward[] = { FFTW_RODFT10, FFTW_HC2R };
 //    const fftw_r2r_kind fftw_kind_backward[] = { FFTW_HC2R, FFTW_HC2R };
  //   plan_Phi2DManyBackward =fftw_plan_many_r2r (2, size, Nz, (double * ) k2In.data() ,NULL, 1, k_size[0] * k_size[1], r3Out.data()             , NULL, 1, size[0]*size[1]    , fftw_kind_backward, FFTW_ESTIMATE);  

 * */

    //////////////////      Set Fourier space transformation for 1D transform (not implemented) but possibly
     if(flags & FFT_X) {
          
          rXIn.resize(RxLD);
          kXOut.resize(Rk2xL);

          plan_Phi1DAvrgForward  = fftw_plan_dft_1d    (Nx, (fftw_complex *) rXIn.data() , (fftw_complex *) kXOut.data(), FFTW_FORWARD , FFTW_ESTIMATE);
          plan_Phi1DAvrgBackward = fftw_plan_dft_1d    (Nx, (fftw_complex *) kXOut.data(), (fftw_complex *) rXIn.data() , FFTW_BACKWARD, FFTW_ESTIMATE);
    }
   /* 
     if(flags & FFT_XYZ) {
    ///////////////////        Set Fourier space domain for 3D transform
        
        K3xLlD = 0;  K3xLuD = Nx-1;
        K3yLlD = 0;  K3yLuD = Ny/2;
        K3zLlD = 0;  K3zLuD = Nz/2;
        
        Rk3xL.setRange(K3xLlD, K3xLuD);
        Rk3yL.setRange(K3yLlD, K3yLuD);
        Rk3zL.setRange(K3zLlD, K3zLuD);
        
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


FFTSolver_fftw3_dst::~FFTSolver_fftw3_dst() {

    if(flags & FFT_X) {
        fftw_destroy_plan(plan_Phi1DAvrgForward);
	    fftw_destroy_plan(plan_Phi1DAvrgBackward);
    }
	
    else if(flags & FFT_XYZ) {
        fftw_destroy_plan(plan_Phi3DForward);
	    fftw_destroy_plan(plan_Phi3DBackward);
    }

    else if(flags & FFT_XY) {
    fftw_destroy_plan(plan_Phi2DManyForward);
	fftw_destroy_plan(plan_Phi2DManyBackward);
    }
    #ifdef PARALLEL_OPENMP
//        fftw_cleanup_threads();
    #endif
  
}

int FFTSolver_fftw3_dst::solve(const int FFTtype, const int direction, const int nstacked) {

/* 
 * */
        Range RnS(1,nstacked);

        if(FFTtype & FFT_XYZ & flags) {
            if     (direction == FFT_FORWARD ) fftw_execute(plan_Phi3DForward  );
            else if(direction == FFT_BACKWARD) fftw_execute(plan_Phi3DBackward);
            else   check(-1, DMESG("No such FFT direction"));
        } 
        else if(FFTtype & FFT_XY & flags) {

            if (direction == FFT_FORWARD ) {

              // for sin : [1 2 3 4 -1 -2 -3 -4]   for cos : [1 2 3 4 1 2 3 4
              if      (DTtype & HELIOS_FFT_DST)  for(int x = NxLlD; x <= NxLuD; x++) r3In(NxLlD+NxLuD+Nx-x, RyLD, RzLD, RnS) = -r3In(x,RyLD, RzLD, RnS);
              else if (DTtype & HELIOS_FFT_DCT)  for(int x = NxLlD; x <= NxLuD; x++) r3In(NxLlD+NxLuD+Nx-x, RyLD, RzLD, RnS) =  r3In(x,RyLD, RzLD, RnS);
              else if (DTtype & HELIOS_FFT_DFT)  for(int x = NxLlD; x <= NxLuD; x++) r3In(NxLlD+NxLuD+Nx-x, RyLD, RzLD, RnS) =  0.;
              else check(-1, DMESG("No such FFT Type"));
             

              // FFT FORWARD
              fftw_execute(plan_Phi2DManyForward);
                  

             // for DST-II and DCT-II, we have a frequncy shift of 1/2 bin, for which we need to correct.
             if((DTtype & HELIOS_FFT_DST) || (DTtype & HELIOS_FFT_DCT)) {

   //             for(int x_k = K2xLlD; x_k <= K2xLuD; x_k++)  k2Out(x_k, Rk2yL, RzLD, RnS) *= cmplxd(0.,1.)/exp(-2.*M_PI*cmplxd(0.,1.)*(1.*((double)Nx)-0.5)*((double) x_k)/(1.*Nx));
                //for(int x_k = K2xLlD; x_k <= K2xLuD; x_k++)  k2Out(x_k, Rk2yL, RzLD, RnS) /= exp(-2.*M_PI*cmplxd(0.,1.)*((2.*(double)Nx)-0.5)*((double) x_k)/(2.*Nx));
                for(int x_k = K2xLlD; x_k <= K2xLuD; x_k++)  k2Out(x_k, Rk2yL, RzLD, RnS) *= 1./exp(-2.*M_PI*cmplxd(0.,1.)*(1.*((double)Nx)-0.5)*((double) x_k)/(1.*Nx));

             }

             // We need to modift real and imaginary values for DCT and DST. For DFT ? Dont know, strides ?
             if(DTtype & HELIOS_FFT_DFT) {
            
             }
             else if(DTtype & HELIOS_FFT_DST) {
/* 
                 real(k2Out(Rk2xL , Rk2yL, RzLD, RnS))  =  0.;
                 imag(k2Out(Rk2xL , Rk2yL, RzLD, RnS)) *= -2.;
                 
                 imag(k2Out(K2xLlD, Rk2yL, RzLD, RnS))  =  0.;
                 imag(k2Out(K2xLuD, Rk2yL, RzLD, RnS))  =  0.;
 * */
//                 real(k2Out(Rk2xL , 0, RzLD, RnS))  =  0.;
//                 imag(k2Out(Rk2xL , 0, RzLD, RnS)) *= -2.;
                 
                  k2Out(0, Rk2yL, RzLD, RnS) =  0.;
        //          k2Out(2*Nx-1, Rk2yL, RzLD, RnS) =  0.;
                 //imag(k2Out(K2xLuD, 0, RzLD, RnS))  =  0.;
//                 k2Out(Nx-1, 0, RzLD, RnS)   =  0.;


            } 
             else if(DTtype & HELIOS_FFT_DCT) {
                 
                 real(k2Out(Rk2xL , Rk2yL, RzLD, RnS)) *=  2.;
                 imag(k2Out(Nx-1  , Rk2yL, RzLD, RnS))  =  0.;
                 imag(k2Out(Rk2xL , Rk2yL, RzLD, RnS))  =  0.;
                 
                 real(k2Out(0   , Rk2yL, RzLD, RnS)) /=  2.;
                 real(k2Out(Nx-1, Rk2yL, RzLD, RnS)) /=  2.;

            }
             else check(-1, DMESG("No such option"));
             
            }
            else if (direction == FFT_BACKWARD) {




              if      (DTtype & HELIOS_FFT_DST)  {
                    //real(k2In(Range(0, 2*Nx-1), Rk2yL, RzLD, RnS)) = 0.;
                    //for(int x_k = K2xLlD+1; x_k <= K2xLuD; x_k++)  imag(k2In(2*Nx-x_k, Rk2yL, RzLD, RnS)) = -imag(k2In(x_k, Rk2yL, RzLD, RnS)); 
//                    real(k2In(Range(0, 2*Nx-1), 0, RzLD, RnS)) = 0.;
        //            for(int x_k = K2xLlD+1; x_k <= K2xLuD; x_k++)  imag(k2In(2*Nx-x_k, 0, RzLD, RnS)) = -imag(k2In(x_k, 0, RzLD, RnS)); 
              }
              else if (DTtype & HELIOS_FFT_DCT)  ;
              else if (DTtype & HELIOS_FFT_DFT)  ;
              else check(-1, DMESG("No such FFT Type"));
            
              fftw_execute(plan_Phi2DManyBackward);
            




            if(DTtype & HELIOS_FFT_DFT) ;
            else if(DTtype & HELIOS_FFT_DST) {

                 // For me its a bit cheating, but it turns out that for lower grid numbers the
                 // upper boundary is less accurate zero. Or is it a BUG ?!
//                 r3Out(NxLlD, RyLD, RzLD, RnS) = 0.;
//                 r3Out(NxLuD, RyLD, RzLD, RnS) = 0.;

            } 
             else if(DTtype & HELIOS_FFT_DCT) {
            }
            else    check(-1, DMESG("No such FFT direction"));
            }
        }
        else if(FFTtype & FFT_X & flags) {
            if      (direction == FFT_FORWARD ) fftw_execute(plan_Phi1DAvrgForward);
            else if (direction == FFT_BACKWARD) fftw_execute(plan_Phi1DAvrgBackward);
            else    check(-1, DMESG("No such FFT direction"));
        } 
        else  check(-1, DMESG("Unknown FFT type or not supported"));

       return HELIOS_SUCCESS;
}


std::string FFTSolver_fftw3_dst::getLibraryName() {
        return std::string("FFTW3 using Dirichlet Boundary Conditions in X (DST)");
}
  
int FFTSolver_fftw3_dst::getDecomposition() {
      return DECOMP_NO;
}

