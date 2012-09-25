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


  
VlasovAux::VlasovAux(Grid *_grid, Parallel *_parallel, Setup *_setup, FileIO *_fileIO, Geometry *_geo, FFTSolver *_fft, Benchmark *_bench)    
        : VlasovCilk(_grid, _parallel, _setup, _fileIO, _geo, _fft, _bench) 
{

  // Cast to two dimensional geometry module (to access k_\parallel)
  geo = static_cast<Geometry2D*>(Vlasov::geo);

}


int VlasovAux::solve(std::string equation_type, Fields *fields, Array6C f_in, Array6C f_out, double dt, int rk_step, const double rk[3]) 
{

  // do I need both, we can stick to e-m ? Speed penality ?
  
  if(equation_type == "ES")

      Vlasov_ES   ((A6zz) f_in.dataZero(), (A6zz) f_out.dataZero()      , (A6zz) f0.dataZero(), (A6zz) f.dataZero(), 
                   (A6zz) ft.dataZero()  , (A6zz) fields->Field.dataZero(), 
                   (A3zz) nonLinearTerms, X, V, M, dt, rk_step, rk);

  else if(equation_type == "EM")

      Vlasov_EM    ((A6zz) f_in.dataZero(), (A6zz) f_out.dataZero(), (A6zz) f0.dataZero(), (A6zz) f.dataZero(),
                   (A6zz) ft.dataZero(), (A6zz) fields->Field.dataZero(), (A3zz) nonLinearTerms,
                   (A4zz) Xi, (A4zz) G, X, V, M, dt, rk_step, rk);

  else if(equation_type == "2DLandauDamping")
    
      Landau_Damping((A6zz) f_in.dataZero(), (A6zz) f_out.dataZero(), (A6zz) f0.dataZero(), (A6zz) f.dataZero(), 
                     (A6zz) ft.dataZero()  , (A6zz) fields->Field.dataZero(), 
                      X, V, M, dt, rk_step, rk);
  
  else   check(-1, DMESG("No Such Equation"));

  return GKC_SUCCESS;
}



void VlasovAux::Vlasov_ES(
                           const CComplex fs        [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           CComplex fss             [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex f0        [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex f1        [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           CComplex ft              [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex Fields[Nq][NsLD][NmLB][NzLB][NkyLD][NxLB+4],
                           CComplex       NL                          [NkyLD][NxLD  ][NvLD],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           const double dt, const int rk_step, const double rk[3])
{ 


  //#pragma ivdep 
  for(int s = NsLlD; s <= NsLuD; s++) {
        
      // some abbrevations
      const double w_n   = plasma->species[s].w_n;
      const double w_T   = plasma->species[s].w_T;
      const double alpha = plasma->species[s].alpha;
      const double sigma = plasma->species[s].sigma;
      const double kw_T  = 1./plasma->species[s].T0;
    
      const double sub = (plasma->species[s].doGyro) ? 3./2. : 1./2.;
        
      const double v2_rms   = 1.;//pow2(alpha);
 
      bool isGyro1 = (plasma->species[s].gyroModel == "Gyro-1");
     
      // Cannot do omp parallelization here (due to nonlinear term ..)
      for(int m=NmLlD; m<= NmLuD;m++) { for(int z=NzLlD; z<= NzLuD;z++) { 
      
      // calculate non-linear term (rk_step == 0 for eigenvalue calculations)
      if(nonLinear && (rk_step != 0)) calculatePoissonBracket(nullptr, nullptr, fs, Fields, z, m, s, NL, Xi_max, false); 
       

      omp_for(int y_k=NkyLlD; y_k <= NkyLuD; y_k++) { 
             
         const CComplex ky = ((CComplex) (0. + 1.j))  * fft->ky(y_k);
             
         for(int x=NxLlD; x<= NxLuD; x++) { 
         
           const CComplex phi_   = Fields[Field::phi][s][m][z][y_k][x];
       
           CComplex half_eta_kperp2_phi = 0;
           if(isGyro1) { // first order approximation for gyro-kinetics
             const CComplex ddphi_dx_dx = (16. *(Fields[Field::phi][s][m][z][y_k][x+1] + Fields[Field::phi][s][m][z][y_k][x-1]) - (Fields[Field::phi][s][m][z][y_k][x+2] + Fields[Field::phi][s][m][z][y_k][x-2]) - 30.*phi_) * _kw_12_dx_dx;
             half_eta_kperp2_phi     = 0.5 * w_T  * ( (ky*ky) * phi_ + ddphi_dx_dx ) ; 
           }
             
           // Sign has no influence on result ...
           const CComplex kp = geo->get_kp(x, ky, z);

         //#pragma unroll(8)
         //#pragma unroll
         //#pragma vector aligned 
         //#pragma vector nontemporal(fss)
         simd_for(int v=NvLlD; v<= NvLuD; v++) { 

            const  CComplex g      = fs [s][m][z][y_k][x][v];
            const  CComplex f0_    = f0 [s][m][z][y_k][x][v];

            const CComplex dfs_dv  = (8.  *(fs[s][m][z][y_k][x][v+1] - fs[s][m][z][y_k][x][v-1]) - (fs[s][m][z][y_k][x][v+2] - fs[s][m][z][y_k][x][v-2]))*_kw_12_dv;
            const CComplex ddfs_dv = (16. *(fs[s][m][z][y_k][x][v+1] + fs[s][m][z][y_k][x][v-1]) - (fs[s][m][z][y_k][x][v+2] + fs[s][m][z][y_k][x][v-2]) - 30.*g) * _kw_12_dv_dv;
        
            // hyperdiffusion terms
            // const Complex d4fs_d4x = (-(fs[s][m][z][y_k][x-2][v] + fs[s][m][z][y_k][x+2][v]) + 4. * (fs[s][m][z][y_k][x-1][v] + fs[s][m][z][y_k][x+1][v]) - 6.*fs[s][m][z][y_k][x][v])/_16_dx4;
       
       
            ///////////////   The time derivative of the Vlasov equation      //////////////////////
       
            const CComplex dg_dt = 
             +  NL[y_k][x][v]                                                         // Non-linear ( array is zero for linear simulations) 
             +  ky* (-(w_n + w_T * (((V[v]*V[v])+ M[m])*kw_T  - sub)) * f0_ * phi_    // Driving term (Temperature/Density gradient)
             -  half_eta_kperp2_phi * f0_)                                            // Contributions from gyro-1 (0 if not neq Gyro-1)
             -  alpha  * V[v]* kp  * ( g + sigma * phi_ * f0_)                        // Linear Landau damping
             +  collisionBeta  * (g  + alpha * V[v] * dfs_dv + 2. * ddfs_dv);         // Lennard-Bernstein Collision term
         
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        
        ft [s][m][z][y_k][x][v] = rk[0] * ft[s][m][z][y_k][x][v] + rk[1] * dg_dt                                ;
        fss[s][m][z][y_k][x][v] = f1[s][m][z][y_k][x][v]         + (rk[2] * ft[s][m][z][y_k][x][v] + dg_dt) * dt;
     
       }
         
         }} }
      
      }
   }
}


void VlasovAux::Vlasov_EM(
                           CComplex fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           CComplex fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex f0 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           CComplex ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex Fields[Nq][NsLD][NmLD][NzLB][NkyLD][NxLB+4]      ,
                           CComplex    nonLinear               [NkyLD][NxLD  ][NvLD],
                           CComplex Xi       [NzLB][NkyLD][NxLB+4][NvLB],
                           CComplex G        [NzLB][NkyLD][NxLB][NvLB],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           const double dt, const int rk_step, const double rk[3])
{ 

   for(int s = NsLlD; s <= NsLuD; s++) {
        
      // abbrevations
      const double w_n   = plasma->species[s].w_n;
      const double w_T   = plasma->species[s].w_T;
      const double alpha = plasma->species[s].alpha;
      const double sigma = plasma->species[s].sigma;
      const double Temp  = plasma->species[s].T0;
    
      const double sub   = (plasma->species[s].doGyro) ? 3./2. : 1./2.;
      
      bool isGyro1 = (plasma->species[s].gyroModel == "Gyro-1");
      

    for(int m=NmLlD; m<= NmLuD;m++) { 
 
//         setupXiAndG(fs, f0 , phi, Ap, Bp, Xi, G, V, M, m , s);
//         if(calculate_nonLinear && (rk_step != 0)) calculatePoissonBracket(Vlasov::Xi, _fs,m, s);
       
         
    for(int z=NzLlD; z<= NzLuD;z++) { omp_for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int x=NxLlD; x<= NxLuD;x++) { 
       
          const CComplex phi_ = (__extension__ 1.0i) * Fields[Field::phi][s][m][z][y_k][x];
               
          const CComplex dphi_dx = (8.*(Fields[Field::phi][s][m][z][y_k][x+1] - Fields[Field::phi][s][m][z][y_k][x-1]) 
                                     - (Fields[Field::phi][s][m][z][y_k][x+2] - Fields[Field::phi][s][m][z][y_k][x-2]))/(12.*dx)  ;  

          const CComplex ky = ((CComplex) (0. + 1.i)) *  fft->ky(y_k);
          const CComplex kp = geo->get_kp(x, ky, z);

   #pragma ivdep
   #pragma vector always 
   for(int v=NvLlD; v<= NvLuD;v++) {
           
           // Sign has no influence on result ...
           const CComplex kp = geo->get_kp(x, ky, z);
        

      const CComplex g    = fs[s][m][z][y_k][x][v];
      const CComplex F0   = f0 [s][m][z][y_k][x][v];
      const CComplex G_   = G[z][y_k][x][v];
      const CComplex Xi_  = Xi[z][y_k][x][v];

      // Velocity derivaties for Lennard-Bernstein Collisional Model
      const CComplex dfs_dv   = (8. *(fs[s][m][z][y_k][x][v+1] - fs[s][m][z][y_k][x][v-1]) - (fs[s][m][z][y_k][x][v+2] - fs[s][m][z][y_k][x][v-2]))/(12.*dv);
      const CComplex ddfs_dvv = (16.*(fs[s][m][z][y_k][x][v+1] + fs[s][m][z][y_k][x][v-1]) - (fs[s][m][z][y_k][x][v+2] + fs[s][m][z][y_k][x][v-2]) - 30.*fs[s][m][z][y_k][x][v])/(12.*dv*dv);
      const double v2_rms = 1.;//pow2(alpha)
     
     
      // calculate first order Average
      CComplex half_eta_kperp2_Xi = 0;
      if(isGyro1) { // first order approximation for gyro-kinetics
             const CComplex ddXi_dx_dx = (16. *(Xi[z][y_k][x+1][v] + Xi[z][y_k][x-1][v])
                                             - (Xi[z][y_k][x+2][v] + Xi[z][y_k][x-2][v]) - 30.*Xi_) * _kw_12_dx_dx;
             half_eta_kperp2_Xi     = 0.5 * w_T  * ( (ky*ky) * Xi_ + ddXi_dx_dx ) ; 
           }
             
    
     /////////////// Finally the Vlasov equation calculate the time derivatve      //////////////////////
   
            
     const CComplex dg_dt = 
             +  nonLinear[y_k][x][v]                                                 // Non-linear ( array is zero for linear simulations) 
             +  ky* (-(w_n + w_T * (((V[v]*V[v])+ M[m])/Temp  - sub)) * F0 * Xi_     // Driving term (Temperature/Density gradient)
             -  half_eta_kperp2_Xi * F0)                                             // Contributions from gyro-1 (0 if not neq Gyro-1)
             -  alpha  * V[v]* kp  * ( g + sigma * Xi_ * F0)                         // Linear Landau damping
             +  collisionBeta  * (g  + alpha * V[v] * dfs_dv + 2. * ddfs_dvv);       // Lennard-Bernstein Collision term
         
        
        //////////////////////////// Vlasov End ////////////////////////////
  
        //  time-integrate the distribution function    
        ft [s][m][z][y_k][x][v] = rk[0] * ft[s][m][z][y_k][x][v] + rk[1] * dg_dt             ;
        fss[s][m][z][y_k][x][v] = f1[s][m][z][y_k][x][v]         + (rk[2] * ft[s][m][z][y_k][x][v] + dg_dt) * dt;



      }}} }}
   }
}

void VlasovAux::Landau_Damping(
                           CComplex fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           CComplex fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex f0 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           CComplex ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex Fields[Field::phi][NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           const double dt, const int rk_step, const double rk[3])
{ 

   for(int s = NsLlD; s <= NsLuD; s++) { for(int m=NmLlD; m<= NmLuD;m++) { 
      
     const double alpha = plasma->species[s].alpha;
  
         for(int z=NzLlD; z<= NzLuD;z++) {       for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) {
     omp_for(int x=NxLlD; x<= NxLuD;x++) {  simd_for(int v=NvLlD; v<= NvLuD;v++)         {
       
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


void    VlasovAux::Vlasov_2D_Fullf(
                           CComplex fs               [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           CComplex fss              [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex f0         [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex f1         [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           CComplex ft               [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const  CComplex Fields[Nq][NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           const double dt, const int rk_step, const double rk[3])
{ 

   for(int s = NsLlD; s <= NsLuD; s++) {
        
      // small abbrevations
      const double w_n   = plasma->species[s].w_n;
      const double w_T   = plasma->species[s].w_T;
      const double alpha = plasma->species[s].alpha;
      const double sigma = plasma->species[s].sigma;
      const double Temp  = plasma->species[s].T0;
    
      const double sub = (plasma->species[s].doGyro) ? 3./2. : 1./2.;
      
      for(int m=NmLlD; m<= NmLuD;m++) { 
  
       // gyro-fluid model
 //      if(calculate_nonLinear && (rk_step != 0)) calculatePoissonBracket(fields->phi, _fs,m, s);
       //if(nonLinear)   calculatePoissonBracket(fields->phi(RxLD, RkyLD, RzLD, m, s), _fs,m, s);
       
       // calculate for estimation of CFL condition
       for(int z=NzLlD; z<= NzLuD;z++) { omp_for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int x=NxLlD; x<= NxLuD;x++) { 
       
             const CComplex dphi_dx  = (8.*(Fields[Field::phi][s][m][z][y_k][x+1] - Fields[Field::phi][s][m][z][y_k][x-1]) - (Fields[Field::phi][s][m][z][y_k][x+2] - Fields[Field::phi][s][m][z][y_k][x-2]))/(12.*dx)  ;  

             const CComplex ky = ((CComplex) (0. + 1.i)) *  fft->ky(y_k);
             const CComplex kp = geo->get_kp(x, ky, z);

            for(int v=NvLlD; v<= NvLuD;v++) {
        



        // abbreviations to not clutter the equations
        const CComplex g_   = fs [s][m][z][y_k][x][v];
        const CComplex f0_  = f0 [s][m][z][y_k][x][v];
        const CComplex phi_ = Fields[Field::phi][s][m][z][y_k][x];

        // Hyper diffusion terms
        //const Complex d4g_dv    =  0.;
//        const Complex d4g_dv    =  -1.e-3 * (-39. *(fs[s][m][z][y_k][x][v+1] - fs[s][m][z][y_k][x][v-1])  + 12. *(fs[s][m][z][y_k][x][v+2] - fs[s][m][z][y_k][x][v-2]) + 56. * fs[s][m][z][y_k][x][v]);///pow4(dv);

        const CComplex dfs_dv    = (8.  *(fs[s][m][z][y_k][x][v+1] - fs[s][m][z][y_k][x][v-1]) - (fs[s][m][z][y_k][x][v+2] - fs[s][m][z][y_k][x][v-2]))/(12.*dv);
        const CComplex ddfs_dvv  = (16. *(fs[s][m][z][y_k][x][v+1] + fs[s][m][z][y_k][x][v-1]) - (fs[s][m][z][y_k][x][v+2] + fs[s][m][z][y_k][x][v-2]) - 30.*fs[s][m][z][y_k][x][v])/(12.*pow2(dv));
        const double v2_rms = 1.;//pow2(alpha)
        ;
        /////////////// Finally the Vlasov equation calculate the time derivatve      //////////////////////
        CComplex dg_dt = 


             
             // driving term (use dphi_dy instead of dXi_dy, because v * A does not vanish due to numerical errors)
             //ky * (-(w_n + w_T * ((pow2(V[v])+ M[m])/Temp  - sub)) * F0 * phi_
             ky * (-(w_T * pow2(V[v]) )* f0_ * phi_
                    -(w_n + w_T * (M[m]/Temp  - sub)) * f0_ * phi_)

             // add first order gyro-average term (zero when full-gyro)
//             + 0.5 * w_T  * k2_Fields[Field::phi][z][y_k][x] * F0

              // Landau Damping term and parallel ... ? - alpha  * V[v]* geo->get_kp(x)  * ( g + sigma * phi * F0)) 
           - alpha  * V[v]* kp * ( g_ + sigma * phi_ * f0_ )
           // Collisional terms 
             + collisionBeta * (g_  + V[v] * dfs_dv + v2_rms * ddfs_dvv)
          ;

            // Energy evolution term
//          + rhoOverLn * 1./mass*(shear * X[x] + theta) * dXi_dy *df1_dv;


        //////////////////////////// Vlasov End ////////////////////////////
        //  time-integrate the distribution function    
        ft [s][m][z][y_k][x][v] = rk[0] * ft[s][m][z][y_k][x][v] + rk[1] * dg_dt             ;
        fss[s][m][z][y_k][x][v] = f1[s][m][z][y_k][x][v]         + (rk[2] * ft[s][m][z][y_k][x][v] + dg_dt) * dt;

      }}} }}
   }
}




