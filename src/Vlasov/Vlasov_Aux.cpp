/*
 * =====================================================================================
 *
 *       Filename:  Vlasov_2D.cpp
 *
 *    Description: Vlasov Solver Implementation for 2D Geometry
 *                 or other special types.
 *
 *         Author: Paul P. Hilscher (2009-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include "Vlasov/Vlasov_Aux.h"


  
VlasovAux::VlasovAux(Grid *_grid, Parallel *_parallel, Setup *_setup, FileIO *_fileIO, Geometry *_geo, FFTSolver *_fft, Benchmark *_bench, Collisions *_coll)    
        : VlasovCilk(_grid, _parallel, _setup, _fileIO, _geo, _fft, _bench, _coll) 
          
{

  // Cast to two dimensional geometry module (to access k_\parallel)
  geo = static_cast<Geometry2D*>(Vlasov::geo);

}


void VlasovAux::solve(std::string equation_type, Fields *fields, CComplex *f_in, CComplex *f_out, double dt, int rk_step, const double rk[3]) 
{

  // do I need both, we can stick to e-m ? Speed penality ?
  
  if(equation_type == "ES")

      Vlasov_ES   ((A6zz) f_in, (A6zz) f_out     , (A6zz) f0, (A6zz) f, 
                   (A6zz) ft , (A6zz) Coll, (A6zz) fields->Field, 
                   (A3zz) nonLinearTerm, X, V, M, dt, rk_step, rk);

  else if(equation_type == "EM")

      Vlasov_EM   ((A6zz) f_in, (A6zz) f_out, (A6zz) f0, (A6zz) f,
                   (A6zz) ft, (A6zz) Coll, (A6zz) fields->Field, (A3zz) nonLinearTerm,
                   (A4zz) Xi, (A4zz) G, dt, rk_step, rk);

  else if(equation_type == "Landau_Damping")
    
      Landau_Damping((A6zz) f_in, (A6zz) f_out, (A6zz) f0, (A6zz) f, 
                     (A6zz) ft , (A6zz) fields->Field, 
                      X, V, M, dt, rk_step, rk);
  
  else   check(-1, DMESG("No Such Equation"));

  return;
}



void VlasovAux::Vlasov_ES(
                           const CComplex fs        [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex fss             [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f0        [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f1        [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex ft              [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Coll      [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Fields[Nq][NsLD][NmLB][NzLB][Nky][NxLB+4],
                           CComplex       nonLinearTerm               [Nky][NxLD  ][NvLD],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           const double dt, const int rk_step, const double rk[3])
{ 


  //#pragma ivdep 
  for(int s = NsLlD; s <= NsLuD; s++) {
        
      // some abbrevations
      const double *w_n  = species[s].w_n;
      const double *w_T  = species[s].w_T;
      const double alpha = species[s].alpha;
      const double sigma = species[s].sigma;
    
      const double sub = (species[s].doGyro) ? 3./2. : 1./2.;
        
      const bool doGyro  = (species[s].doGyro);
      const bool isGyro1 = (species[s].gyroModel == "Gyro-1");
    
      const double rho_t2 = species[s].T0 * species[s].m / (pow2(species[s].q) * plasma->B0); 
     
      // Cannot do omp parallelization here (due to nonlinear term ..)
      for(int m = NmLlD; m <= NmLuD; m++) { for(int z = NzLlD; z <= NzLuD; z++) { 
      
      // calculate non-linear term (rk_step == 0 for eigenvalue calculations)
     // #pragma omp single
     //   {
      if(doNonLinear         && (rk_step != 0)) calculatePoissonBracket(nullptr, nullptr, fs, Fields, z, m, s, nonLinearTerm, Xi_max, false); 
      if(doNonLinearParallel && (rk_step != 0)) calculateParallelNonLinearity2(fs, Fields, z, m, s, nonLinearTerm);
      //  }
      //#pragma omp flush

      #pragma omp for
      for(int y_k=NkyLlD; y_k <= NkyLuD; y_k++) { 
             
         const CComplex ky = _imag  * fft->ky(y_k);
             
         for(int x=NxLlD; x<= NxLuD; x++) { 
      
           const double kw_T  = 1./species[s].T[x];
         
           const CComplex phi_   = Fields[Field::phi][s][m][z][y_k][x];
       
           CComplex half_eta_kperp2_phi = 0;
           if(isGyro1) { // first order approximation for gyro-kinetics
             const CComplex ddphi_dx_dx = (16. *(Fields[Field::phi][s][m][z][y_k][x+1] + Fields[Field::phi][s][m][z][y_k][x-1]) -
                                                (Fields[Field::phi][s][m][z][y_k][x+2] + Fields[Field::phi][s][m][z][y_k][x-2]) 
                                         - 30.*phi_) * _kw_12_dx_dx;
             half_eta_kperp2_phi     = rho_t2 * 0.5 * w_T[x]  * ( (ky*ky) * phi_ + ddphi_dx_dx ) ; 
           }
             
           // Sign has no influence on result ...
           const CComplex kp = geo->get_kp(x, ky, z);

         //#pragma unroll(8)
         //#pragma unroll
         //#pragma vector aligned 
         //#pragma vector nontemporal(fss)
         simd_for(int v=NvLlD; v<= NvLuD; v++) { 
           
            const CComplex   dg_dv = (8.*(fs[s][m][z][y_k][x][v+1] - fs[s][m][z][y_k][x][v-1]) 
                                       - (fs[s][m][z][y_k][x][v+2] - fs[s][m][z][y_k][x][v-2]))/(12.*dv);

            const  CComplex f0_    = f0 [s][m][z][y_k][x][v];

            const  CComplex g      = fs [s][m][z][y_k][x][v];
        
            ///////////////   The time derivative of the Vlasov equation      //////////////////////
       
            const CComplex dg_dt =

             -  nonLinearTerm[y_k][x][v]                                                   // Non-linear ( array is zero for linear simulations) 
             +  ky* (-(w_n[x] + w_T[x] * (((V[v]*V[v])+ (doGyro ? M[m] : 0.))*kw_T  - sub)) * f0_ * phi_    // Driving term (Temperature/Density gradient)
             -  half_eta_kperp2_phi * f0_)                                            // Contributions from gyro-1 (0 if not neq Gyro-1)
             -  alpha  * V[v]* kp  * ( g + sigma * phi_ * f0_)                        // Linear Landau damping
             +  Coll[s][m][z][y_k][x][v]  ;                                           // Collisional operator
         
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        
        ft [s][m][z][y_k][x][v] = rk[0] * ft[s][m][z][y_k][x][v] + rk[1] * dg_dt                                ;
        fss[s][m][z][y_k][x][v] = f1[s][m][z][y_k][x][v]         + (rk[2] * ft[s][m][z][y_k][x][v] + dg_dt) * dt;
     
       }
         
         }} }
      
      }
   }
}


void VlasovAux::Vlasov_EM(
                           CComplex fs       [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex fss      [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f0 [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f1 [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex ft       [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Coll      [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Fields[Nq][NsLD][NmLD][NzLB][Nky][NxLB+4]      ,
                           CComplex    nonLinearTerm               [Nky][NxLD  ][NvLD],
                           CComplex Xi       [NzLB][Nky][NxLB+4][NvLB],
                           CComplex G        [NzLB][Nky][NxLB][NvLB],
                           const double dt, const int rk_step, const double rk[3])
{ 

   for(int s = NsLlD; s <= NsLuD; s++) {
        
      // abbrevations
      const double *w_n   = species[s].w_n;
      const double *w_T   = species[s].w_T;
      const double alpha = species[s].alpha;
      const double sigma = species[s].sigma;
    
      const double sub   = (species[s].doGyro) ? 3./2. : 1./2.;
      
      bool isGyro1 = (species[s].gyroModel == "Gyro-1");
      

    for(int m = NmLlD; m <= NmLuD; m++) { 
 
          setupXiAndG(fs, f0 , Fields, Xi, G, m , s);
       
    for(int z = NzLlD; z <= NzLuD; z++) { 
      
          if(doNonLinear && (rk_step != 0)) calculatePoissonBracket(G, Xi, nullptr, nullptr, z, m, s, nonLinearTerm, Xi_max, true); 
      
    for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int x = NxLlD; x <= NxLuD; x++) { 
      
      const double _kw_T  = 1./species[s].T[x];

      const CComplex phi_ = Fields[Field::phi][s][m][z][y_k][x];
               
      const CComplex dphi_dx = (8.*(Fields[Field::phi][s][m][z][y_k][x+1] - Fields[Field::phi][s][m][z][y_k][x-1]) 
                                 - (Fields[Field::phi][s][m][z][y_k][x+2] - Fields[Field::phi][s][m][z][y_k][x-2])) * _kw_12_dx;

      const CComplex ky = _imag *  fft->ky(y_k);
      const CComplex kp = geo->get_kp(x, ky, z);

    #pragma ivdep
    #pragma vector always 
    for(int v = NvLlD; v <= NvLuD; v++) {
           
      // sign has no influence on result ...
      const CComplex kp = geo->get_kp(x, ky, z);

      const CComplex g    = fs[s][m][z][y_k][x][v];
      const CComplex F0   = f0[s][m][z][y_k][x][v];

      const CComplex G_   =  G[z][y_k][x][v];
      const CComplex Xi_  = Xi[z][y_k][x][v];

    /////////////// Finally the Vlasov equation calculate the time derivatve      //////////////////////
   
            
    const CComplex dg_dt = 
    
    +  nonLinearTerm[y_k][x][v]                                             // Non-linear ( array is zero for linear simulations) 
    -  ky* (w_n[x] + w_T[x] * (((V[v]*V[v])+ M[m]) * _kw_T  - sub)) * F0 * Xi_     // Driving term (Temperature/Density gradient)
    -  alpha  * V[v]* kp  * G_                                              // Linear Landau damping
    +  Coll[s][m][z][y_k][x][v]  ;                                          // Collisional operator
         
        
    //////////////////////////// Vlasov End ////////////////////////////
  
    //  time-integrate the distribution function    
    ft [s][m][z][y_k][x][v] = rk[0] * ft[s][m][z][y_k][x][v] +  rk[1] * dg_dt             ;
    fss[s][m][z][y_k][x][v] =         f1[s][m][z][y_k][x][v] + (rk[2] * ft[s][m][z][y_k][x][v] + dg_dt) * dt;



   } } } 
   } }
   }
}

void VlasovAux::Landau_Damping(
                           CComplex fs       [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex fss      [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f0 [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f1 [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex ft       [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Fields[Field::phi][NsLD][NmLD][NzLB][Nky][NxLB+4],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           const double dt, const int rk_step, const double rk[3])
{ 

   for(int s = NsLlD; s <= NsLuD; s++) { for(int m=NmLlD; m<= NmLuD;m++) { 
      
     const double alpha = species[s].alpha;
  
         for(int z=NzLlD; z<= NzLuD;z++) {       for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) {
         for(int x=NxLlD; x<= NxLuD;x++) {  simd_for(int v=NvLlD; v<= NvLuD;v++)         {
       
           const CComplex f0_     = f0 [s][m][z][y_k][x][v];

           const CComplex dphi_dx = (8.*(Fields[Field::phi][s][m][z][y_k][x+1] - Fields[Field::phi][s][m][z][y_k][x-1]) 
                                      - (Fields[Field::phi][s][m][z][y_k][x+2] - Fields[Field::phi][s][m][z][y_k][x-2]))/(12.*dx)  ;  
           
           const CComplex   df_dv = (8.*(fs[s][m][z][y_k][x][v+1] - fs[s][m][z][y_k][x][v-1]) 
                                      - (fs[s][m][z][y_k][x][v+2] - fs[s][m][z][y_k][x][v-2]))/(12.*dv);

           const CComplex   df_dx = (8. *(fs[s][m][z][y_k][x+1][v] - fs[s][m][z][y_k][x-1][v])  
                                       - (fs[s][m][z][y_k][x+2][v] - fs[s][m][z][y_k][x-2][v]))/(12.*dx) ;
           
           /////////////// 2D (x,v) Landau Damping Test       //////////////////////
        
           const CComplex dg_dt = - alpha * V[v] * (df_dx +  f0_ * dphi_dx)   // linear Landau damping term
                                  + dphi_dx  * df_dv;                         // non-linear Landau term

           //////////////////////////// Landau End ////////////////////////////

           //  time-integrate the distribution function    
           ft [s][m][z][y_k][x][v] = rk[0] * ft[s][m][z][y_k][x][v] + rk[1] * dg_dt             ;
           fss[s][m][z][y_k][x][v] = f1[s][m][z][y_k][x][v]         + (rk[2] * ft[s][m][z][y_k][x][v] + dg_dt) * dt;

      }}} }}
   }
}



void VlasovAux::calculateParallelNonLinearity(
                                const CComplex f          [NsLD][NmLD][NzLB][Nky][NxLB   ][NvLB],
                                const CComplex Fields [Nq][NsLD][NmLD ][NzLB][Nky][NxLB+4], // in case of e-s
                                const int z, const int m, const int s                     ,
                                CComplex nonLinearTerm[Nky][NxLD][NvLD])
{
  
  const double half_vth_kw_T = 0.5 * species[s].q * species[s].v_th / species[s].T0; 

  #pragma omp for
  for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { 
  
    const CComplex ky   = _imag  * fft->ky(y_k);
          
  for(int x = NxLlD; x <= NxLuD; x++) { 
 
    const CComplex kp   = geo->get_kp(x, ky, z);
    const CComplex phi_ = Fields[Field::phi][s][m][z][y_k][x];

  simd_for(int v = NvLlD; v <= NvLuD; v++) { 

    const CComplex  dg_dv = (8.*(f[s][m][z][y_k][x][v+1] - f[s][m][z][y_k][x][v-1]) 
                              - (f[s][m][z][y_k][x][v+2] - f[s][m][z][y_k][x][v-2])) * _kw_12_dv;

    //  Kinetic energy can only be contain in the m=0 mode, otherwise it will
    //  cancel out. Thus we direclty only calculte the m=0 mode which is given
    //  by  (m=0) = (m') + (m'') = 0 -> m' = -m'' = m''*.
    //  m = 0 needs to be per definition real
    //
    //  m_sigma is absorbed in dg_dv
    //
    nonLinearTerm[0][x][v] -= half_vth_kw_T * creal((kp * phi_) * conj(dg_dv));
               
   
  } } }


};

void VlasovAux::calculateParallelNonLinearity2(
                                const CComplex f          [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                                const CComplex Fields [Nq][NsLD][NmLD][NzLB][Nky][NxLB+4], 
                                const int z, const int m, const int s                    ,
                                CComplex nonLinearTerm[Nky][NxLD][NvLD])
{
  nonLinearTerm[:][NxLlD:NxLD][NvLlD:NvLD] = 0.; 
  
  const double _kw_fft_Norm  = 1./(fft->Norm_Y_Backward * fft->Norm_Y_Backward * fft->Norm_Y_Forward);
  const double half_vth_kw_T = 0.5 * species[s].q *  pow2(species[s].v_th) / species[s].T0; 
 

  // Note : we do not use any anti-aliasing method here!

  CComplex Arr_NkyNx[Nky][NxLD], NL_NxNky[Nky][NxLD];
  
  doubleAA xy_phi_kp[NyLD][NxLD], xy_dg_dv[NyLD][NxLD], NL_NxNy[NyLD][NxLD];

  for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) {  for(int x = NxLlD; x <= NxLuD; x++) { 
 
    const CComplex ky   = _imag  * fft->ky(y_k);
    const CComplex kp   = geo->get_kp(x, ky, z);
    const CComplex phi_ = Fields[Field::phi][s][m][z][y_k][x];

    Arr_NkyNx[y_k][x-NxLlD] = kp * phi_ ;

  } }

  fft->solve(FFT_Type::Y_NL, FFT_Sign::Backward, (CComplex *) Arr_NkyNx, &xy_phi_kp);

  for(int v = NvLlD; v <= NvLuD; v++) {

  for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { for(int x = NxLlD; x <= NxLuD; x++) { 
  
    Arr_NkyNx[y_k][x-NxLlD] = (8.*(f[s][m][z][y_k][x][v+1] - f[s][m][z][y_k][x][v-1]) 
                                - (f[s][m][z][y_k][x][v+2] - f[s][m][z][y_k][x][v-2])) * _kw_12_dv;
 
  } } // y_k, x
  
  fft->solve(FFT_Type::Y_NL, FFT_Sign::Backward, (CComplex *) Arr_NkyNx, &xy_dg_dv);

  // multiply non-linear terms
  NL_NxNy[:][:] = _kw_fft_Norm * half_vth_kw_T * xy_phi_kp[:][:] * xy_dg_dv[:][:] ;
   
  fft->solve(FFT_Type::Y_NL, FFT_Sign::Forward, NL_NxNy, (CComplex *) NL_NxNky);

  nonLinearTerm[:][NxLlD:NxLD][v] -= NL_NxNky[:][:]; 
  
  } // v


};

