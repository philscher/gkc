/*
 * =====================================================================================
 *
 *       Filename: Vlasov_Optim.cpp
 *
 *    Description: Implementation of optimized algorithms for Vlasov 
 *
 *         Author: Paul P. Hilscher (2012-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */
//#define INSTRSET 7            // SSE2 required
//#include "external/VectorClass/vectorclass.h"
//#define INSTRSET 7            // SSE2 required
//#include "external/VectorClass/special/complexvec.h"
//#define INSTRSET 7            // SSE2 required

#include "Vlasov_Optim.h"
#include "Geometry/Geometry2D.h"




VlasovOptim::VlasovOptim(Grid *_grid, Parallel *_parallel, Setup *_setup, FileIO *fileIO, Geometry *_geo, FFTSolver *fft, Benchmark *_bench, Collisions *_coll)    
: Vlasov(_grid, _parallel, _setup, fileIO, _geo, fft, _bench, _coll)
{
   // no need to call : Is done is base constructor -- or ?    Vlasov::initData(fileIO);    
}

void VlasovOptim::solve(std::string equation_type, Fields *fields, CComplex *_fs, CComplex *_fss,
                        double dt, int rk_step, const double rk[3]) 
{
 
  if((equation_type == "VlasovAux_ES")) {

      Vlasov_2D((A6sz) _fs, (A6sz) _fss, (A6sz) f0, 
                (A6sz) f, (A6sz) ft, (A6sz) Coll, (A6sz) fields->Field,
                (A4sz) nonLinearTerm, X, V, M, dt, rk_step, rk);
  }

  return;
}

                           

 void    VlasovOptim::Vlasov_2D(
                           const cmplx16 fs        [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           cmplx16 fss             [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const cmplx16 f0        [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const cmplx16 f1        [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           cmplx16 ft              [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const cmplx16 Coll      [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const cmplx16 Fields[Nq][NsLD][NmLD][NzLB][Nky][NxLB+4]      ,
                           cmplx16 nonLinearTerm               [NzLD][Nky][NxLD][NvLD]  ,
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           const double dt, const int rk_step, const double rk[3])
{ 
 
  const Geometry2D *geo = static_cast<Geometry2D*>(Vlasov::geo);


  const int BlockSize_X=bench->BlockSize_X ;
  const int BlockSize_V=bench->BlockSize_V ;

  // MAXIMUM_ALIGN preprocessor directive !
   __assume_aligned(fs    , 64);
   __assume_aligned(fss   , 64);
   __assume_aligned(f1    , 64);
   __assume_aligned(f0    , 64);
   __assume_aligned(ft    , 64);
   __assume_aligned(Fields, 64);
   __assume(BlockSize_V % 4 == 0);
   __assume(BlockSize_X % 2 == 0);
   
   for(int s = NsLlD; s <= NsLuD; s++) {
        
      const double w_n   = species[s].w_n;
      const double w_T   = species[s].w_T;
      const double alpha = species[s].alpha;
      const double sigma = species[s].sigma;
      const double kw_T  = 1./species[s].T0;
    
      const double sub = (species[s].doGyro) ? 3./2. : 1./2.;
        
      for(int m=NmLlD; m<= NmLuD;m++) { 
  
       for(int z=NzLlD; z<= NzLuD;z++) { for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { 
             
          const double ky_im = fft->ky(y_k);
        
       for(int x=NxLlD; x<= NxLuD; x+=BlockSize_X) { 
       for(int v=NvLlD; v<= NvLuD; v+=BlockSize_V) { 


       // Block Loop
       for(int xx=x; xx<= x+BlockSize_X;xx++) {

             //const double Krook_nu = krook_nu * ( (X[x] > 0.8 * Lx/2.) || (X[x] < -0.8 * Lx/2.))  ?  0.1 * pow2(abs(X[x]) - 0.8 * Lx/2.): 0.;

             //const Complex dphi_dx  = (8.*(Fields[Field::phi][s][m][z][y_k][xx+1] - Fields[Field::phi][s][m][z][y_k][xx-1]) - (Fields[Field::phi][s][m][z][y_k][xx+2] - Fields[Field::phi][s][m][z][y_k][xx-2]))/(12.*dx)  ;  
              const double phi_re     = Fields[Field::phi][s][m][z][y_k][xx].re;
              const double phi_im     = Fields[Field::phi][s][m][z][y_k][xx].im;
             
             // BUG : does it needs a mutex ?
              const double kp = 0.4 * X[x] * (2.*M_PI)/Ly * y_k; //geo->get_kp(xx, ky, z);
             //updateCFL(dphi_dx, ky*phi_, 0.);


           //  Complex4d Phi; Phi.load_a((double*)&Fields[Field::phi][s][m][z][y_k][xx]);
            //for(int v=NvLlD; v<= NvLuD-4;v+=4) {
         
            //#pragma unroll(8)
            //#pragma unroll
            #pragma nofusion
            #pragma ivdep
            #pragma  vector nontemporal(fss)
            simd_for(int vv=v; vv<= v+BlockSize_V;vv++) {
  
       //        const  Complex g    = fs [s][m][z][y_k][xx][vv];
       //        const  Complex F0   = f0 [s][m][z][y_k][xx][vv];

       const double dfs_dv_re   = (8.  *(fs[s][m][z][y_k][xx][vv+1].re - fs[s][m][z][y_k][xx][vv-1].re) - (fs[s][m][z][y_k][xx][vv+2].re - fs[s][m][z][y_k][xx][vv-2].re))*_kw_12_dv;
       const double ddfs_dvv_re = (16. *(fs[s][m][z][y_k][xx][vv+1].re + fs[s][m][z][y_k][xx][vv-1].re) - (fs[s][m][z][y_k][xx][vv+2].re + fs[s][m][z][y_k][xx][vv-2].re) - 30.*fs[s][m][z][y_k][xx][vv].re) * _kw_12_dv_dv;
     
       /////////////// Finally the Vlasov equation calculate the time derivative      //////////////////////
       
       // Note, we split here real and imaginary part to enhance vectorization
       const double dg_dt_re = 
            ky_im* (-(w_n + w_T * (((V[vv]*V[vv])+ M[m])/kw_T  - 3./2.)) * f0 [s][m][z][y_k][xx][vv].re * phi_im )
             - alpha  * V[vv]* kp  * ( fs[s][m][z][y_k][xx][vv].im + sigma * phi_im * f0 [s][m][z][y_k][xx][vv].re)
             + Coll[s][m][z][y_k][xx][vv].re;
        
        ft [s][m][z][y_k][xx][vv].re = rk[0] * ft[s][m][z][y_k][xx][vv].re + rk[1] * dg_dt_re                                ;
        fss[s][m][z][y_k][xx][vv].re = f1[s][m][z][y_k][xx][vv].re         + (rk[2] * ft[s][m][z][y_k][xx][vv].re + dg_dt_re) * dt;
       
       const double ddfs_dvv_im = (16. *(fs[s][m][z][y_k][xx][vv+1].im + fs[s][m][z][y_k][xx][vv-1].im) 
                                      - (fs[s][m][z][y_k][xx][vv+2].im + fs[s][m][z][y_k][xx][vv-2].im) 
                                   - 30.*fs[s][m][z][y_k][xx][vv].im) * _kw_12_dv_dv;

       const double dfs_dv_im   = (8.  *(fs[s][m][z][y_k][xx][vv+1].im - fs[s][m][z][y_k][xx][vv-1].im) 
                                      - (fs[s][m][z][y_k][xx][vv+2].im - fs[s][m][z][y_k][xx][vv-2].im))
                                  * _kw_12_dv;
       
        const  double dg_dt_im = 
            ky_im* (-(w_n + w_T * (((V[vv]*V[vv])+ M[m])/kw_T  - 3./2.)) * f0 [s][m][z][y_k][xx][vv].re * phi_re )
             - alpha  * V[vv]* kp  * ( fs[s][m][z][y_k][xx][vv].re + sigma * phi_re * f0 [s][m][z][y_k][xx][vv].re)
             + Coll[s][m][z][y_k][xx][vv].im;
        
        ft [s][m][z][y_k][xx][vv].im = rk[0] * ft[s][m][z][y_k][xx][vv].im + rk[1] * dg_dt_im                                ;
        fss[s][m][z][y_k][xx][vv].im = f1[s][m][z][y_k][xx][vv].im         + (rk[2] * ft[s][m][z][y_k][xx][vv].im + dg_dt_im) * dt;

     
       } } // Block-Loop
         
       } } 
         
       } } 
       } }
}


