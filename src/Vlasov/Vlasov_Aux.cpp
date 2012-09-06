/*
 * =====================================================================================
 *
 *       Filename:  Vlasov_2D.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  09/01/2012 04:06:01 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include "Vlasov/Vlasov_Aux.h"


// Replace asap with C99 function or so
extern CComplex conj(CComplex c) {
   
     struct _Z{ double re; double im; };

     _Z *a = (_Z*) &c ;
     a->im = - a->im;

    return c; 
};
 
  
VlasovAux::VlasovAux(Grid *_grid, Parallel *_parallel, Setup *_setup, FileIO *fileIO, Geometry *_geo, FFTSolver *fft, Benchmark *_bench)    
: VlasovCilk(_grid, _parallel, _setup, fileIO, _geo, fft, _bench) 
{

  // Cast to two dimwnsional geometry module (to access k_\parallel)
  geo = static_cast<Geometry2D*>(Vlasov::geo);

}


int VlasovAux::solve(std::string equation_type, Fields *fields, Array6C f_in, Array6C f_out, double dt, int rk_step, const double rk[3]) 
{

  // do I need both, we can stick to e-m ? Speed penality ?
  if(equation_type == "ES")

      Vlasov_ES   ((A6zz) f_in.dataZero(), (A6zz) f_out.dataZero()      , (A6zz) f0.dataZero(), (A6zz) f.dataZero(), 
                   (A6zz) ft.dataZero() , (A5zz) fields->phi.dataZero(), 
                   (A3zz) nonLinearTerms.dataZero(), X.dataZero(), V.dataZero(), M.dataZero(), dt, rk_step, rk);

  else if(equation_type == "EM")

      Vlasov_EM    ((A6zz) f_in.dataZero(), (A6zz) f_out.dataZero(), (A6zz) f0.dataZero(), (A6zz) f.dataZero(),
                   (A6zz) ft.dataZero(), (A5zz) fields->phi.dataZero(), (A5zz) fields->Ap.dataZero(),
                   (A5zz) fields->Bp.dataZero(), (A3zz) nonLinearTerms.dataZero(),
                   (A4zz) Xi.dataZero(), (A4zz) G.dataZero(), X.dataZero(), V.dataZero(), M.dataZero(), dt, rk_step, rk);

  else if(equation_type == "2D_Island") 
    
      Vlasov_2D_Island((A6zz) f_in.dataZero(), (A6zz) f_out.dataZero(), (A6zz) f0.dataZero(), (A6zz) f.dataZero(), 
                    (A6zz) ft.dataZero(), (A5zz) fields->phi.dataZero(), (A3zz) nonLinearTerms.dataZero(),
                    X.dataZero(), V.dataZero(), M.dataZero(), dt, rk_step, rk);

//  else if(equation_type == "2DLandauDamping") Landau_Damping((A6z) f_in.dataZero(), (A6z) f_out.dataZero(), (A6z) f0.dataZero(), (A6z) f.dataZero(), (A6z) ft.dataZero(), (A5z) fields->phi.dataZero(), (A4z) k2p_phi.dataZero(), (A3z) dphi_dx.dataZero(), X.dataZero(), V.dataZero(), M.dataZero(), fields, dt, rk_step, rk);
  else   check(-1, DMESG("No Such Equation"));

  return GKC_SUCCESS;
}



void VlasovAux::Vlasov_ES(
                           const CComplex fs [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           CComplex fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex f0 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           CComplex ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex phi[NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           CComplex       NL                   [NkyLD][NxLD  ][NvLD],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           const double dt, const int rk_step, const double rk[3])
{ 

  // KW is Kehrwert (German for Multiplicative Inverse)
  // increases speed by 5%, set directly in Vlasov ? Arrangment ?  
  const double _kw_12_dv    = 1./(12.*dv);
  const double _kw_12_dv_dv = 1./(12.*dv*dv);
  const double _kw_16_dx4   = 1./(16.*pow4(dx));
  const double _kw_12_dx_dx = 1./(12.*dx*dx);

  #pragma ivdep 
  for(int s = NsLlD; s <= NsLuD; s++) {
        
      // some abbrevations
      const double w_n   = plasma->species(s).w_n;
      const double w_T   = plasma->species(s).w_T;
      const double alpha = plasma->species(s).alpha;
      const double sigma = plasma->species(s).sigma;
      const double kw_T  = 1./plasma->species(s).T0;
    
      const double sub = (plasma->species(s).doGyro) ? 3./2. : 1./2.;
        
      const double v2_rms   = 1.;//pow2(alpha);
  
      bool isGyro1 = (plasma->species(s).gyroModel == "Gyro-1");

      
      for(int m=NmLlD; m<= NmLuD;m++) { for(int z=NzLlD; z<= NzLuD;z++) { 
      
         // calculate non-linear term (rk_step == 0 for eigenvalue calculations)
         if(nonLinear && (rk_step != 0)) calculatePoissonBracket(nullptr, nullptr, fs, phi, z, m, s, NL, Xi_max, false); 
       

      omp_for(int y_k=NkyLlD; y_k <= NkyLuD; y_k++) { 
             
         //const CComplex ky = _Imaginary * fft->ky(y_k);
         const CComplex ky = ((CComplex) (0. + 1.j))  * fft->ky(y_k);
             
         for(int x=NxLlD; x<= NxLuD; x++) { 
       
           CComplex eta_kp2_phi = 0;
           if(isGyro1) { // first order approximation for gyro-kinetics
             const CComplex ddphi_dx_dx = (16. *(phi[s][m][z][y_k][x+1] + phi[s][m][z][y_k][x-1]) - (phi[s][m][z][y_k][x+2] + phi[s][m][z][y_k][x-2]) - 30.*phi[s][m][z][y_k][x]) * _kw_12_dx_dx;
             eta_kp2_phi                = 0.5 * w_T  * ( ky*ky * phi[s][m][z][y_k][x] + ddphi_dx_dx ); 
           }
             
           const CComplex kp = geo->get_kp(x, ky, z);

         //#pragma unroll(8)
         //#pragma unroll
         #pragma  vector nontemporal(fss)
         simd_for(int v=NvLlD; v<= NvLuD; v++) { 



             const CComplex phi_     = phi[s][m][z][y_k][x];
             
             //updateCFL(dphi_dx, ky*phi_, 0.);
       

            const  CComplex g    = fs [s][m][z][y_k][x][v];
            const  CComplex f0_  = f0 [s][m][z][y_k][x][v];

            const CComplex dfs_dv   = (8.  *(fs[s][m][z][y_k][x][v+1] - fs[s][m][z][y_k][x][v-1]) - (fs[s][m][z][y_k][x][v+2] - fs[s][m][z][y_k][x][v-2]))*_kw_12_dv;
            const CComplex ddfs_dv = (16. *(fs[s][m][z][y_k][x][v+1] + fs[s][m][z][y_k][x][v-1]) - (fs[s][m][z][y_k][x][v+2] + fs[s][m][z][y_k][x][v-2]) - 30.*fs[s][m][z][y_k][x][v]) * _kw_12_dv_dv;
        
            // hyperdiffusion terms
            //     const Complex d4fs_d4x = (-(fs[s][m][z][y_k][x-2][v] + fs[s][m][z][y_k][x+2][v]) + 4. * (fs[s][m][z][y_k][x-1][v] + fs[s][m][z][y_k][x+1][v]) - 6.*fs[s][m][z][y_k][x][v])/_16_dx4;
       
       
            ///////////////   The time derivative of the Vlasov equation      //////////////////////
       
            const CComplex dg_dt = 
                NL[y_k][x][v]                                                         // Non-linear ( array is zero for linear simulations) 
             +  ky* (-(w_n + w_T * (((V[v]*V[v])+ M[m])*kw_T  - sub)) * f0_ * phi_    // Driving term (Temperature/Density gradient)
             +  eta_kp2_phi * f0_ )                                                   // Contributions from gyro-1 (0 if not neq Gyro-1)
             -  alpha  * V[v]* kp  * ( g + sigma * phi_ * f0_)                        // Linear Landau damping
             +  collisionBeta  * (g  + alpha * V[v] * dfs_dv + v2_rms * ddfs_dv);     // Lennard-Bernstein Collision term
         
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        
        ft [s][m][z][y_k][x][v] = rk[0] * ft[s][m][z][y_k][x][v] + rk[1] * dg_dt                                ;
        fss[s][m][z][y_k][x][v] = f1[s][m][z][y_k][x][v]         + (rk[2] * ft[s][m][z][y_k][x][v] + dg_dt) * dt;
     
       }
         
         }} }
      
      }
   }
}

void VlasovAux::Vlasov_2D_Island(
                           CComplex fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           CComplex fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex f0 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           CComplex ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex phi[NsLD][NmLD][NzLB][NkyLD][NxLB+4]      ,
                           CComplex nonLinear                  [NkyLD][NxLD  ][NvLD],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           const double dt, const int rk_step, const double rk[3])
{ 

    const double w     = setup->get("Island.Width", 0.); 
    const double shear = setup->get("Geometry.Shear", 0.4); 
  
    const CComplex imag = ((CComplex) (0. + 1.j));


    for(int s = NsLlD; s <= NsLuD; s++) {
        
      // small abbrevations
      const double w_n   = plasma->species(s).w_n;
      const double w_T   = plasma->species(s).w_T;
      const double alpha = plasma->species(s).alpha;
      const double sigma = plasma->species(s).sigma;
      const double Temp  = plasma->species(s).T0;
      const double sub   = (plasma->species(s).doGyro) ? 3./2. : 1./2.;

      const double v2_rms = 1.;//pow2(alpha);


      for(int m=NmLlD; m<=NmLuD; m++) {

         // calculate non-linear term (rk_step == 0 for eigenvalue calculations)
         //if(nonLinear && (rk_step != 0)) calculatePoissonBracket(G, Xi, g, phi, z, m, s, ExB, Xi_max, true); 

      for(int z=NzLlD; z<= NzLuD;z++) {  omp_for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) {



            // Note : for negative modes we need to use complex conjugate value

            const CComplex ky     = imag * fft->ky(y_k);

            
            // We need to take care of boundaries. For polidal numbers y_k > N_k-1, we  use zero.
            // For y_k < 0, the corresponding complex conjugate value is used.
            const CComplex ky_p1  = (y_k == Nky-1) ? 0.                 :  imag * fft->ky(y_k+1);
            const CComplex ky_m1  = (y_k == 0    ) ? imag * -fft->ky(1) : imag * fft->ky(y_k-1); 
            const CComplex ky_1   = imag * fft->ky(1);

          simd_for(int x=NxLlD; x<= NxLuD;x++) {  

               // calculate for estimation of CFL condition
             const CComplex phi_ = phi[s][m][z][y_k][x];

             const CComplex dphi_dx  = (8.*(phi[s][m][z][y_k][x+1] - phi[s][m][z][y_k][x-1]) - (phi[s][m][z][y_k][x+2] - phi[s][m][z][y_k][x-2]))/(12.*dx)  ;  

             updateCFL(dphi_dx, ky*phi_, 0.);
             
        /////////////////////////////////////////////////// Magnetic Island Contribution    /////////////////////////////////////////
        
        // NOTE :  at the Nyquist frequency we have no coupling with higher frequencies (actually phi(m=Ny) = 0. anyway)

        const double zeta      = 2.*M_PI/Ly;
        
        //const double MagIs     = -  0.5 * w*w*shear/16. * cos(zeta * X[x]);
        //const double dMagIs_dx = +  0.5 * w*w*shear/16. * sin(zeta * X[x]) * zeta;
        const double p[] = { 0.13828847,  0.70216594, -0.01033686 };
        const double xx = pow2(X[x]);
        
        const double psi  = (1. + p[0]*pow(xx,p[1])) * exp( p[2] * xx);
        const double dpsi = (X[x] == 0.) ? 0. :  (p[0] * 2. * X[x] * p[1]*pow(xx, p[1]-1.) + (1. + p[0] * pow(xx,p[1])) * p[2] * 2. * X[x]) * exp(p[2]*xx);


        //const double MagIs     = - 0.5 * w*w*shear/16.  * psi ;//cos(zeta * X[x]);
        //const double dMagIs_dx = - 0.5 * w*w*shear/16. * dpsi;// - sin(zeta * X[x]) * zeta;
        // changed sign , check !!!
        const double MagIs     = - 0.5 * w*w*shear/16.  * psi ;//cos(zeta * X[x]);
        const double dMagIs_dx = - 0.5 * w*w*shear/16. * dpsi;// - sin(zeta * X[x]) * zeta;

        
        const CComplex     phi_p1 = ( y_k == Nky-1) ? 0.                       : phi[s][m][z][y_k+1][x] ;
        const CComplex     phi_m1 = ( y_k ==  0   ) ? conj(phi[s][m][z][1][x]) : phi[s][m][z][y_k-1][x] ;


        const CComplex dphi_dx_p1 = ( y_k == Nky-1) ? 0. 
                                                  : (8.*(phi[s][m][z][y_k+1][x+1] - phi[s][m][z][y_k+1][x-1]) - (phi[s][m][z][y_k+1][x+2] - phi[s][m][z][y_k+1][x-2]))/(12.*dx)  ;

        const CComplex dphi_dx_m1 = ( y_k ==    0 ) ?  conj(8.*(phi[s][m][z][    1][x+1] - phi[s][m][z][    1][x-1]) - (phi[s][m][z][    1][x+2] - phi[s][m][z][    1][x-2]))/(12.*dx) 
                                                  :      (8.*(phi[s][m][z][y_k-1][x+1] - phi[s][m][z][y_k-1][x-1]) - (phi[s][m][z][y_k-1][x+2] - phi[s][m][z][y_k-1][x-2]))/(12.*dx) ;
        
       
        // The magnetic island
    
        /*
                *  For the latter term, the intrigate derivative is the \partial_y * Island
                *  remember the island structure is 
                *  \partial_y (e^{imx} + e^{-imx}) = (i m) * ( e^{imx} - e^{-imx} )
                *
        */
        const CComplex Island_A_phi =   dMagIs_dx * ( ky_m1 * phi_m1 + ky_p1 * phi_p1) - MagIs *  ky_1 * ( dphi_dx_m1 -  dphi_dx_p1);
        
             ///////////////////////////////////////////////////////////////////////////////
            
        //const Complex kp = 0. ; //geo->get_kp(x, ky, z);
        const CComplex kp = geo->shear * ky * X[x];

               // velocity space magic
        simd_for(int v=NvLlD; v<= NvLuD;v++) {

            const CComplex g   = fs[s][m][z][y_k][x][v];
            const CComplex F0  = f0[s][m][z][y_k][x][v];


        /////////////////////////////////////////////////// Magnetic Island Contribution    /////////////////////////////////////////
      
        const CComplex dfs_dx_p1  =  (y_k == Nky-1) 
                            ? 0.
                            : (8. *(fs[s][m][z][y_k+1][x+1][v] - fs[s][m][z][y_k+1][x-1][v])  - (fs[s][m][z][y_k+1][x+2][v] - fs[s][m][z][y_k+1][x-2][v]))/(12.*dx) ;

        const CComplex dfs_dx_m1  =  ( y_k == 0   ) 
                            ?     conj(8. *(fs[s][m][z][    1][x+1][v] - fs[s][m][z][    1][x-1][v])  - (fs[s][m][z][    1][x+2][v] - fs[s][m][z][    1][x-2][v]))/(12.*dx) 
                            :       (8. *(fs[s][m][z][y_k-1][x+1][v] - fs[s][m][z][y_k-1][x-1][v])  - (fs[s][m][z][y_k-1][x+2][v] - fs[s][m][z][y_k-1][x-2][v]))/(12.*dx) ;

        // Note Nky-1 is the maximum mode number Nky = 6 i-> [ 0, 1, 2, 3, 4, 5] 
        const CComplex fs_p1      = (y_k == Nky-1) ? 0.                         : fs[s][m][z][y_k+1][x][v] ;
        const CComplex fs_m1      = (y_k ==  0   ) ? conj(fs[s][m][z][1][x][v]) : fs[s][m][z][y_k-1][x][v] ;
         
        // not at the Nyquist frequency we have no coupling with higher frequencies
   
        // mode-mode connections
        register const CComplex Island_A_F1 =  dMagIs_dx * (ky_m1 * fs_m1  + ky_p1 * fs_p1 )  -  MagIs  * ky_1 *  (dfs_dx_m1  - dfs_dx_p1 )  ;


   
        /////////// Collisions ////////////////////////////////////////////////////////////////////

        const CComplex dfs_dv    = (8.  *(fs[s][m][z][y_k][x][v+1] - fs[s][m][z][y_k][x][v-1]) - 1. *(fs[s][m][z][y_k][x][v+2] - fs[s][m][z][y_k][x][v-2]))/(12.*dv);
        const CComplex ddfs_dvv  = (16. *(fs[s][m][z][y_k][x][v+1] + fs[s][m][z][y_k][x][v-1]) - 1. *(fs[s][m][z][y_k][x][v+2] + fs[s][m][z][y_k][x][v-2]) - 30.*fs[s][m][z][y_k][x][v])/(12.*pow2(dv));


        /////////////// Finally the Vlasov equation calculate the time derivatve      //////////////////////
        const CComplex dg_dt = 
             // Island
             - alpha * V[v] * (Island_A_F1 + sigma * Island_A_phi * F0) +   
             // driving term
             ky* (-(w_n + w_T * (((V[v]*V[v])+ M[m])/Temp  - sub)) * F0 * phi_

             // add first order gyro-average term (zero when full-gyro)
                 )//    + 0.5 * w_T  * k2_phi[1][z][y_k][x] * F0)
              // Landau Damping term 
           - alpha  * V[v]* kp  * ( g + sigma * phi_ * F0);
            // collisional term
            + collisionBeta * (g  + V[v] * dfs_dv + v2_rms * ddfs_dvv)
           + nonLinear[y_k][x][v]
          ;

        // screen out Nyquiest frequeny
        //if(y_k == Ny/2) dg_dt = 0.;


        //////////////////////////// Vlasov End ////////////////////////////
        //  time-integrate the distribution function    
        ft [s][m][z][y_k][x][v] = rk[0] * ft[s][m][z][y_k][x][v] + rk[1] * dg_dt             ;
        fss[s][m][z][y_k][x][v] = f1[s][m][z][y_k][x][v]         + (rk[2] * ft[s][m][z][y_k][x][v] + dg_dt) * dt;
      
        }}} }}
   }

}

void VlasovAux::Vlasov_EM(
                           CComplex fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           CComplex fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex f0 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           CComplex ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex phi[NsLD][NmLD][NzLB][NkyLD][NxLB+4]      ,
                           const CComplex Ap [NsLD][NmLD][NzLB][NkyLD][NxLB+4]      ,
                           const CComplex Bp [NsLD][NmLD][NzLB][NkyLD][NxLB+4]      ,
                           CComplex    nonLinear               [NkyLD][NxLD  ][NvLD],
                           CComplex Xi       [NzLB][NkyLD][NxLB+4][NvLB],
                           CComplex G        [NzLB][NkyLD][NxLB][NvLB],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           const double dt, const int rk_step, const double rk[3])
{ 


   for(int s = NsLlD; s <= NsLuD; s++) {
        
      // abbrevations
      const double w_n   = plasma->species(s).w_n;
      const double w_T   = plasma->species(s).w_T;
      const double alpha = plasma->species(s).alpha;
      const double sigma = plasma->species(s).sigma;
      const double Temp  = plasma->species(s).T0;
    
      const double sub   = (plasma->species(s).doGyro) ? 3./2. : 1./2.;
      

      for(int m=NmLlD; m<= NmLuD;m++) { 
 
//         setupXiAndG(fs, f0 , phi, Ap, Bp, Xi, G, V, M, m , s);

//         if(calculate_nonLinear && (rk_step != 0)) calculatePoissonBracket(Vlasov::Xi, _fs,m, s);
       
         // calculate for estimation of CFL condition
         for(int z=NzLlD; z<= NzLuD;z++) { omp_for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int x=NxLlD; x<= NxLuD;x++) { 
       
               const CComplex phi_ = (__extension__ 1.0i) * phi[s][m][z][y_k][x];
               //const CComplex phi_ = _Imaginary * phi[s][m][z][y_k][x];
               const CComplex dphi_dx = (8.*(phi[s][m][z][y_k][x+1] - phi[s][m][z][y_k][x-1]) - (phi[s][m][z][y_k][x+2] - phi[s][m][z][y_k][x-2]))/(12.*dx)  ;  

               const CComplex ky = ((CComplex) (0. + 1.i)) *  fft->ky(y_k);
               const CComplex kp = geo->get_kp(x, ky, z);

               updateCFL(dphi_dx, ky*phi_, 0.);


   #pragma ivdep
   #pragma vector always 
   for(int v=NvLlD; v<= NvLuD;v++) {
        

      const CComplex g    = fs[s][m][z][y_k][x][v];
      const CComplex F0   = f0 [s][m][z][y_k][x][v];
      const CComplex G_   = G[z][y_k][x][v];
      const CComplex Xi_  = Xi[z][y_k][x][v];

      // Velocity derivaties for Lennard-Bernstein Collisional Model
      const CComplex dfs_dv   = (8. *(fs[s][m][z][y_k][x][v+1] - fs[s][m][z][y_k][x][v-1]) - (fs[s][m][z][y_k][x][v+2] - fs[s][m][z][y_k][x][v-2]))/(12.*dv);
      const CComplex ddfs_dvv = (16.*(fs[s][m][z][y_k][x][v+1] + fs[s][m][z][y_k][x][v-1]) - (fs[s][m][z][y_k][x][v+2] + fs[s][m][z][y_k][x][v-2]) - 30.*fs[s][m][z][y_k][x][v])/(12.*dv*dv);
      const double v2_rms = 1.;//pow2(alpha)
    
      //Complex k2_Xi = 0.;
      /* 
      if(plasma->species(s).gyroModel == "Gyro-1") { 
        const Complex ddXi_dx = (16.*(Xi[z][y_k][x+1][v] + Xi[z][y_k][x-1][v]) - (Xi[z][y_k][x+2][v] + Xi[z][y_k][x-2][v]) - 30.*Xi[z][y_k][x][v])/(12.*dx*dx);
        k2_Xi = (ky*ky + ddXi_dx);
     }
       * */
   
     /////////////// Finally the Vlasov equation calculate the time derivatve      //////////////////////
   
     const CComplex dg_dt = 
             
      // nonLinearTerm
      nonLinear[y_k][x][v] +
             
      // driving term (use dphi_dy instead of dXi_dy, because v * A does not vanish due to numerical errors)
      //ky* (-(w_n + w_T * (((V[v]*V[v])+ M[m])/Temp  - sub)) * F0 * Xi_
      ky* (-(w_n + w_T * (((V[v]*V[v])+ M[m])/Temp  - sub)) * F0 * phi_
          )//             + 0.5 * w_T  * k2_phi[1][z][y_k][x] * F0)
      // add first order gyro-average term (zero when full-gyro)
      // Landau Damping term and parallel ... ? 
      - alpha  * V[v]* kp * G_
      // Collisional terms 
      + collisionBeta * (g  + V[v] * dfs_dv + v2_rms * ddfs_dvv)    ;

        //////////////////////////// Vlasov End ////////////////////////////
  
      //  time-integrate the distribution function    
      ft [s][m][z][y_k][x][v] = rk[0] * ft[s][m][z][y_k][x][v] + rk[1] * dg_dt             ;
      fss[s][m][z][y_k][x][v] = f1[s][m][z][y_k][x][v]         + (rk[2] * ft[s][m][z][y_k][x][v] + dg_dt) * dt;



      }}} }}
   }
}

/*
 *
void VlasovAux::Vlasov_2D_Island(
                           Complex fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           Complex fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const Complex f0 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const Complex f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           Complex ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const Complex phi[NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           Complex k2_phi[plasma->nfields][NzLD][NkyLD][NxLD],
                           Complex nonLinear[NzLD][NkyLD][NxLD][NvLD],
                           Complex dphi_dx[NzLB][NkyLD][NxLB],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           const double dt, const int rk_step, Array6C _fs)
{ 


            const Complex ky     = Complex(0.,fft->ky(y_k));
            const Complex ky_p1  = (y_k+1 <= Nky-1) ? Complex(0.,fft->ky(y_k+1)) : 0.;
            const Complex ky_m1  = (y_k-1 >= 0    ) ? Complex(0.,fft->ky(y_k-1)) : 0.; 
            const Complex ky_1    = Complex(0., fft->ky(1));

          //#pragma simd
          for(int x=NxLlD; x<= NxLuD;x++) {  

        
               // calculate for estimation of CFL condition
             const Complex phi_ = phi[s][m][z][y_k][x];
             dphi_dx[z][y_k][x]  = (8.*(phi[s][m][z][y_k][x+1] - phi[s][m][z][y_k][x-1]) - (phi[s][m][z][y_k][x+2] - phi[s][m][z][y_k][x-2]))/(12.*dx)  ;  

             updateCFL(dphi_dx[z][y_k][x], ky*phi_, 0.);
             
        /////////////////////////////////////////////////// Magnetic Island Contribution    /////////////////////////////////////////

        const double zeta=2.*M_PI/Ly;
        const double MagIs     =  0.5 * w*w*shear/16. * cos(zeta * X[x]);
        const double dMagIs_dx = -  0.5 * w*w*shear/16. * sin(zeta * X[x]) * zeta;

        // Note : f, fftw ordering [ 0, 1, 2, 3, 4, -3, -2, -1 ]
        const Complex     phi_p1 = (( y_k+1) <= Nky-1) ? phi[s][m][z][y_k+1][x] : 0.;
        const Complex     phi_m1 = (( y_k-1) >=    0) ? phi[s][m][z][y_k-1][x] : 0.;
        const Complex dphi_dx_p1 = (( y_k+1) <= Nky-1) ?
                    (8.*(phi[s][m][z][y_k+1][x+1] - phi[s][m][z][y_k+1][x-1]) - (phi[s][m][z][y_k+1][x+2] - phi[s][m][z][y_k+1][x-2]))/(12.*dx)  : 0.;
        const Complex dphi_dx_m1 = (( y_k-1) >=    0) ?
                    (8.*(phi[s][m][z][y_k-1][x+1] - phi[s][m][z][y_k-1][x-1]) - (phi[s][m][z][y_k-1][x+2] - phi[s][m][z][y_k-1][x-2]))/(12.*dx)  : 0.;
        
        // NOTE :  at the Nyquist frequency we have no coupling with higher frequencies (actually phi(m=Ny) = 0. anyway)
       
        const Complex Island_A_phi =  ( - dMagIs_dx * ( ky_m1 * phi_m1 + ky_p1 * phi_p1) + MagIs * ky_1 * (dphi_dx_m1 - dphi_dx_p1));
        
             ///////////////////////////////////////////////////////////////////////////////
            
        const Complex kp = geo->get_kp(x, ky, z);


        /////////////////////////////////////////////////// Magnetic Island Contribution    /////////////////////////////////////////
      
        // how does they couple with Nyquist frequencty ??   // take care of boundaries !!
        const Complex dfs_dx_p1  =  ((y_k+1) <= Nky-1) ?
                                    (8. *(fs[s][m][z][y_k+1][x+1][v] - fs[s][m][z][y_k+1][x-1][v])  - (fs[s][m][z][y_k+1][x+2][v] - fs[s][m][z][y_k+1][x-2][v]))/(12.*dx) : 0.;

        const Complex dfs_dx_m1  =  ((y_k-1) >= 0   ) ? 
                                    (8. *(fs[s][m][z][y_k-1][x+1][v] - fs[s][m][z][y_k-1][x-1][v])  - (fs[s][m][z][y_k-1][x+2][v] - fs[s][m][z][y_k-1][x-2][v]))/(12.*dx) : 0.;

        const Complex fs_p1      = ((y_k+1) <= Nky-1) ? fs[s][m][z][y_k+1][x][v] : 0.;
        const Complex fs_m1      = ((y_k-1) >= 0   ) ? fs[s][m][z][y_k-1][x][v] :  0.;
         

        // not at the Nyquist frequency we have no coupling with higher frequencies
   
        // mode-mode connections
        register const Complex Island_A_F1 = ( - dMagIs_dx * (ky_m1 *fs_m1  + ky_p1 * fs_p1 )  +  MagIs  * ky_1 *  (dfs_dx_m1  - dfs_dx_p1 ))  ;
   
        /////////////// Finally the Vlasov equation calculate the time derivatve      //////////////////////
        Complex dg_dt = 
             // Island
             alpha * V[v] * (Island_A_F1 + Island_A_phi * F0) +   
*/

void VlasovAux::Landau_Damping(
                           CComplex fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           CComplex fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex f0 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           CComplex ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex phi[NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           const double dt, const int rk_step, const double rk[3])
{ 

   for(int s = NsLlD; s <= NsLuD; s++) { for(int m=NmLlD; m<= NmLuD;m++) { 
      
     const double alpha = plasma->species(s).alpha;
  
       for(int z=NzLlD; z<= NzLuD;z++) {       for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) {
   omp_for(int x=NxLlD; x<= NxLuD;x++) {  simd_for(int v=NvLlD; v<= NvLuD;v++)         {
       
           const CComplex f0_     = f0 [s][m][z][y_k][x][v];
           const CComplex dphi_dx = (8.*(phi[s][m][z][y_k][x+1] - phi[s][m][z][y_k][x-1]) - (phi[s][m][z][y_k][x+2] - phi[s][m][z][y_k][x-2]))/(12.*dx)  ;  
           const CComplex   df_dv = (8.*(fs[s][m][z][y_k][x][v+1] - fs[s][m][z][y_k][x][v-1]) - (fs[s][m][z][y_k][x][v+2] - fs[s][m][z][y_k][x][v-2]))/(12.*dv);
           const CComplex   df_dx = (8. *(fs[s][m][z][y_k][x+1][v] - fs[s][m][z][y_k][x-1][v])  - (fs[s][m][z][y_k][x+2][v] - fs[s][m][z][y_k][x-2][v]))/(12.*dx) ;
           
           /////////////// 2D (x,v) Landau Damping Test       //////////////////////
        
           const CComplex dg_dt = - alpha * V[v] * (df_dx +  f0_ * dphi_dx)   // linear damping term
                                  + dphi_dx  * df_dv;                        // non-linear term

           //////////////////////////// Landau End ////////////////////////////

        //  time-integrate the distribution function    
        ft [s][m][z][y_k][x][v] = rk[0] * ft[s][m][z][y_k][x][v] + rk[1] * dg_dt             ;
        fss[s][m][z][y_k][x][v] = f1[s][m][z][y_k][x][v]         + (rk[2] * ft[s][m][z][y_k][x][v] + dg_dt) * dt;

      }}} }}
   }
}




void    VlasovAux::Vlasov_2D_Fullf(
                           Complex fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           Complex fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const Complex f0 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const Complex f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           Complex ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const Complex phi[NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           const double dt, const int rk_step, const double rk[3])
{ 

   for(int s = NsLlD; s <= NsLuD; s++) {
        
      // small abbrevations
      const double w_n   = plasma->species(s).w_n;
      const double w_T   = plasma->species(s).w_T;
      const double alpha = plasma->species(s).alpha;
      const double sigma = plasma->species(s).sigma;
      const double Temp  = plasma->species(s).T0;
    
      const double sub = (plasma->species(s).doGyro) ? 3./2. : 1./2.;
      
      for(int m=NmLlD; m<= NmLuD;m++) { 
  
       // gyro-fluid model
 //      if(calculate_nonLinear && (rk_step != 0)) calculatePoissonBracket(fields->phi, _fs,m, s);
       //if(nonLinear)   calculatePoissonBracket(fields->phi(RxLD, RkyLD, RzLD, m, s), _fs,m, s);
       
       // calculate for estimation of CFL condition
       for(int z=NzLlD; z<= NzLuD;z++) { omp_for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int x=NxLlD; x<= NxLuD;x++) { 
       
             const Complex dphi_dx  = (8.*(phi[s][m][z][y_k][x+1] - phi[s][m][z][y_k][x-1]) - (phi[s][m][z][y_k][x+2] - phi[s][m][z][y_k][x-2]))/(12.*dx)  ;  

             const Complex ky = Complex(0., fft->ky(y_k));
             const Complex kp = 0.; //geo->get_kp(x, ky, z);

  //           Xi_max(DIR_X) = max(Xi_max(DIR_X), abs(dphi_dx(x,y_k,z)));
  //           Xi_max(DIR_Y) = max(Xi_max(DIR_Y), abs(ky*phi[s][m][z][y_k][x])); 
 //            Xi_max(DIR_Z) = max(Xi_max(DIR_Z), abs(kp));// * abs(phi[s][m][z][y_k][x])); 


            for(int v=NvLlD; v<= NvLuD;v++) {
        



        // abbreviations to not clutter the equations
        const Complex g_   = fs [s][m][z][y_k][x][v];
        const Complex f0_  = f0 [s][m][z][y_k][x][v];
        const Complex phi_ = phi[s][m][z][y_k][x];

        // Hyper diffusion terms
        //const Complex d4g_dv    =  0.;
//        const Complex d4g_dv    =  -1.e-3 * (-39. *(fs[s][m][z][y_k][x][v+1] - fs[s][m][z][y_k][x][v-1])  + 12. *(fs[s][m][z][y_k][x][v+2] - fs[s][m][z][y_k][x][v-2]) + 56. * fs[s][m][z][y_k][x][v]);///pow4(dv);

        const Complex dfs_dv    = (8.  *(fs[s][m][z][y_k][x][v+1] - fs[s][m][z][y_k][x][v-1]) - (fs[s][m][z][y_k][x][v+2] - fs[s][m][z][y_k][x][v-2]))/(12.*dv);
        const Complex ddfs_dvv  = (16. *(fs[s][m][z][y_k][x][v+1] + fs[s][m][z][y_k][x][v-1]) - (fs[s][m][z][y_k][x][v+2] + fs[s][m][z][y_k][x][v-2]) - 30.*fs[s][m][z][y_k][x][v])/(12.*pow2(dv));
        const double v2_rms = 1.;//pow2(alpha)
        ;
        /////////////// Finally the Vlasov equation calculate the time derivatve      //////////////////////
        Complex dg_dt = 


             nonLinearTerms(x,y_k,v) +
             
             // driving term (use dphi_dy instead of dXi_dy, because v * A does not vanish due to numerical errors)
             //ky * (-(w_n + w_T * ((pow2(V(v))+ M(m))/Temp  - sub)) * F0 * phi_
             ky * (-(w_T * pow2(V(v)) )* f0_ * phi_
                    -(w_n + w_T * (M(m)/Temp  - sub)) * f0_ * phi_)

             // add first order gyro-average term (zero when full-gyro)
//             + 0.5 * w_T  * k2_phi[z][y_k][x] * F0

              // Landau Damping term and parallel ... ? - alpha  * V(v)* geo->get_kp(x)  * ( g + sigma * phi * F0)) 
           - alpha  * V(v)* kp * ( g_ + sigma * phi_ * f0_ )
           // Collisional terms 
             + collisionBeta * (g_  + V(v) * dfs_dv + v2_rms * ddfs_dvv)
          ;

         

// non-linear Landau damping is not included, but when is it important ?
          // simple pitch angle scattering in n without
          //+ Coll(x,y,z,v,m,s);
          //
            // Energy evolution term
//          + rhoOverLn * 1./mass*(shear * X(x) + theta) * dXi_dy *df1_dv;


        //////////////////////////// Vlasov End ////////////////////////////
        //  time-integrate the distribution function    
        ft [s][m][z][y_k][x][v] = rk[0] * ft[s][m][z][y_k][x][v] + rk[1] * dg_dt             ;
        fss[s][m][z][y_k][x][v] = f1[s][m][z][y_k][x][v]         + (rk[2] * ft[s][m][z][y_k][x][v] + dg_dt) * dt;

      }}} }}
   }
}




