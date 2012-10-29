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

#include "Vlasov/Vlasov_Island.h"

  
VlasovIsland::VlasovIsland(Grid *_grid, Parallel *_parallel, Setup *_setup, FileIO *fileIO, Geometry *_geo, FFTSolver *fft, Benchmark *_bench, Collisions *_coll)    
: VlasovAux(_grid, _parallel, _setup, fileIO, _geo, fft, _bench, _coll) 
{

    /// Set Magnetic Island structure
    width = setup->get("Island.Width"  , 0.); 
    shear = setup->get("Geometry.Shear", 0.4); 

    

    //// Setup Magnetic Island structure
    nct::allocate(nct::Range(NxGlB, NxGB))(&MagIs, &dMagIs_dx);

    for(int x=NxLlD; x<= NxLuD;x++) {  

    const double zeta      = 2.*M_PI/Ly;
       
    // Island cos(k_T y) = 0.5 [ exp(-i m) + exp( i m ) ]
    // The poloidal island term is direclty included per 
    // mode-mode coupling
    
    // we used simple fitting model to get least square coefficient
    const double p[] = { 0.13828847,  0.70216594, -0.01033686 };

    const double xx   = pow2(X[x]);
    const double psi  = (1. + p[0]*pow(xx,p[1])) * exp( p[2] * xx);
    const double dpsi = (X[x] == 0.) ? 0. :  (p[0] * 2. * X[x] * p[1]*pow(xx, p[1]-1.) + (1. + p[0] * pow(xx,p[1])) * p[2] * 2. * X[x]) * exp(p[2]*xx);

    
    MagIs[x]     = 0.5 * width*width*shear/16.  * psi ;//cos(zeta * X[x]);
    dMagIs_dx[x] = 0.5 * width*width*shear/16. * dpsi;// - sin(zeta * X[x]) * zeta;
    //MagIs[x]     = 0.5 * width*width*shear/16.  * cos(zeta * X[x]);
    //dMagIs_dx[x] = 0.5 * width*width*shear/16. *  (-sin(zeta * X[x])) * zeta;
 
    }

}


void VlasovIsland::solve(std::string equation_type, Fields *fields, CComplex *f_in, CComplex *f_out, double dt, int rk_step, const double rk[3]) 
{

  // do I need both, we can stick to e-m ? Speed penality ?
  if(0) ;  
  
  else if(equation_type == "2D_Island") 
    
      Vlasov_2D_Island((A6zz) f_in, (A6zz) f_out, (A6zz) f0, (A6zz) f, 
                       (A6zz) ft  , (A6zz) Coll, (A6zz) fields->Field, (A3zz) nonLinearTerms,
                       MagIs, dMagIs_dx, X, V, M, dt, rk_step, rk);
  
  else   check(-1, DMESG("No Such Equation"));



  return;

}

void VlasovIsland::Vlasov_2D_Island(
                           CComplex fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           CComplex fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex f0 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           CComplex ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex Coll      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex Fields[Nq][NsLD][NmLD][NzLB][NkyLD][NxLB+4]      ,
                           CComplex nonLinear                  [NkyLD][NxLD  ][NvLD],
                           const double MagIs[NxGB], const double dMagIs[NxGB], 
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           const double dt, const int rk_step, const double rk[3])
{ 

  
    const CComplex imag = ((CComplex) (0. + 1.j));

    for(int s = NsLlD; s <= NsLuD; s++) {
        
      // small abbrevations
      const double w_n   = plasma->species[s].w_n;
      const double w_T   = plasma->species[s].w_T;
      const double alpha = plasma->species[s].alpha;
      const double sigma = plasma->species[s].sigma;
      const double Temp  = plasma->species[s].T0;
      const double sub   = (plasma->species[s].doGyro) ? 3./2. : 1./2.;

      const double v2_rms = 1.;//pow2(alpha);


      const double kw_T = 1./Temp;

      bool isGyro1 = (plasma->species[s].gyroModel == "Gyro-1");
      
      const double rho_t2 = plasma->species[s].T0 * plasma->species[s].m / (pow2(plasma->species[s].q) * plasma->B0); 


      for(int m=NmLlD; m<=NmLuD; m++) { for(int z=NzLlD; z<= NzLuD;z++) {  
        
        //calculate non-linear term (rk_step == 0 for eigenvalue calculations)
        if(nonLinear && (rk_step != 0)) calculatePoissonBracket(nullptr, nullptr, fs, Fields, z, m, s, nonLinear, Xi_max, false); 
        
        omp_for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) {

        
        // Note : for negative modes we need to use complex conjugate value
        const CComplex ky     = imag * fft->ky(y_k);
        
        // We need to take care of boundaries. For poloidal numbers y_k > N_k-1, we  use zero.
        // For y_k < 0, the corresponding complex conjugate value is used.
        const CComplex ky_p1  = (y_k == Nky-1) ? 0.                  : imag * fft->ky(y_k+1);
        const CComplex ky_m1  = (y_k == 0    ) ? - imag * fft->ky(1) : imag * fft->ky(y_k-1); 
        const CComplex ky_1   = imag * fft->ky(1);

        
        for(int x=NxLlD; x<= NxLuD;x++) {  

          
          const CComplex phi_ = Fields[Field::phi][s][m][z][y_k][x];
          
          const CComplex dphi_dx  = (8.*(Fields[Field::phi][s][m][z][y_k][x+1] - Fields[Field::phi][s][m][z][y_k][x-1])
                                      - (Fields[Field::phi][s][m][z][y_k][x+2] - Fields[Field::phi][s][m][z][y_k][x-2])) * _kw_12_dx  ;  

          /////////////////////////////////////////////////// Magnetic Island Contribution    /////////////////////////////////////////
        
          // NOTE :  at the Nyquist frequency we have no coupling with higher frequencies (actually phi(m=Ny) = 0. anyway)

          const CComplex     phi_p1 = ( y_k == Nky-1) ? 0.                                      : Fields[Field::phi][s][m][z][y_k+1][x] ;
          const CComplex     phi_m1 = ( y_k ==  0   ) ? conj(Fields[Field::phi][s][m][z][1][x]) : Fields[Field::phi][s][m][z][y_k-1][x] ;

          
          // X-derivative (First Derivative with Central Difference 4th) of phi for poloidal mode +1, take care of Nyquist frequency
          const CComplex dphi_dx_p1 = ( y_k == Nky-1) 
                                   ? 0.

                                   :      (8.*(Fields[Field::phi][s][m][z][y_k+1][x+1] - Fields[Field::phi][s][m][z][y_k+1][x-1]) 
                                            - (Fields[Field::phi][s][m][z][y_k+1][x+2] - Fields[Field::phi][s][m][z][y_k+1][x-2])) * _kw_12_dx  ;

          // X-derivative (1st deriv. CD-4 )of phi for poloidal mode -1, take care of complex conjugate relation for y_k=-1
          const CComplex dphi_dx_m1 = ( y_k ==    0 ) 
                                   ?  conj(8.*(Fields[Field::phi][s][m][z][    1][x+1] - Fields[Field::phi][s][m][z][    1][x-1]) 
                                            - (Fields[Field::phi][s][m][z][    1][x+2] - Fields[Field::phi][s][m][z][    1][x-2])) * _kw_12_dx

                                   :      (8.*(Fields[Field::phi][s][m][z][y_k-1][x+1] - Fields[Field::phi][s][m][z][y_k-1][x-1]) 
                                            - (Fields[Field::phi][s][m][z][y_k-1][x+2] - Fields[Field::phi][s][m][z][y_k-1][x-2])) * _kw_12_dx ;
        
       
        // The magnetic island
    
        //  For the latter term, the intrigate derivative is the \partial_y * Island
        //  remember the island structure is 
        //  \partial_y (e^{imx} + e^{-imx}) = (i m) * ( e^{imx} - e^{-imx} )
        //
        const CComplex Island_phi =   dMagIs_dx[x] * ( ky_m1 * phi_m1 + ky_p1 * phi_p1) - MagIs[x] *  ky_1 * ( dphi_dx_m1 -  dphi_dx_p1);
        
             ///////////////////////////////////////////////////////////////////////////////
            
        const CComplex kp = geo->get_kp(x, ky, z);
           
        CComplex half_eta_kperp2_phi = 0;

        if(isGyro1) { // first order approximation for gyro-kinetics
             const CComplex ddphi_dx_dx = (16. *(Fields[Field::phi][s][m][z][y_k][x+1] + Fields[Field::phi][s][m][z][y_k][x-1]) 
                                              - (Fields[Field::phi][s][m][z][y_k][x+2] + Fields[Field::phi][s][m][z][y_k][x-2]) 
                                         - 30.*phi_) * _kw_12_dx_dx;

             half_eta_kperp2_phi     = rho_t2 * 0.5 * w_T  * ( (ky*ky) * phi_ + ddphi_dx_dx ) ; 
         }
             
        
        // velocity space magic
        simd_for(int v=NvLlD; v<= NvLuD;v++) {

            const CComplex g    = fs[s][m][z][y_k][x][v];
            const CComplex f0_  = f0[s][m][z][y_k][x][v];


        /////////////////////////////////////////////////// Magnetic Island Contribution    /////////////////////////////////////////
      
          
        // X-derivative of f1 (1-CD4) for poloidal mode +1, take care of Nyquist frequency
        const CComplex dfs_dx_p1  =  (y_k == Nky-1)

                            ? 0.

                            : (8. *(fs[s][m][z][y_k+1][x+1][v] - fs[s][m][z][y_k+1][x-1][v])  
                                 - (fs[s][m][z][y_k+1][x+2][v] - fs[s][m][z][y_k+1][x-2][v])) * _kw_12_dx;


        // X-derivative of f1 (1-CD4) for poloidal mode -1, take care of complex conjugate relation for y_k=-1 
        const CComplex dfs_dx_m1  =  ( y_k == 0   )

                        ? conj(8. *(fs[s][m][z][    1][x+1][v] - fs[s][m][z][    1][x-1][v])  
                                 - (fs[s][m][z][    1][x+2][v] - fs[s][m][z][    1][x-2][v])) * _kw_12_dx 

                        :     (8. *(fs[s][m][z][y_k-1][x+1][v] - fs[s][m][z][y_k-1][x-1][v])  
                                 - (fs[s][m][z][y_k-1][x+2][v] - fs[s][m][z][y_k-1][x-2][v])) * _kw_12_dx;

        // Note Nky-1 is the maximum mode number Nky = 6 i-> [ 0, 1, 2, 3, 4, 5] 
        const CComplex fs_p1      = (y_k == Nky-1) ? 0.                         : fs[s][m][z][y_k+1][x][v] ;
        const CComplex fs_m1      = (y_k ==  0   ) ? conj(fs[s][m][z][1][x][v]) : fs[s][m][z][y_k-1][x][v] ;
         
        // Coupling of phase-space with Island mode-mode coupling
        const CComplex Island_g =  dMagIs_dx[x] * (ky_m1 * fs_m1  + ky_p1 * fs_p1 )  -  MagIs[x]  * ky_1 *  (dfs_dx_m1  - dfs_dx_p1 )  ;

        
        /////////// Collisions ////////////////////////////////////////////////////////////////////
/*
        const CComplex dfs_dv    = (8.  *(fs[s][m][z][y_k][x][v+1] - fs[s][m][z][y_k][x][v-1]) - 
                                         (fs[s][m][z][y_k][x][v+2] - fs[s][m][z][y_k][x][v-2])) * _kw_12_dv;

        const CComplex ddfs_dvv  = (16. *(fs[s][m][z][y_k][x][v+1] + fs[s][m][z][y_k][x][v-1]) -
                                         (fs[s][m][z][y_k][x][v+2] + fs[s][m][z][y_k][x][v-2]) 
                                    - 30.*fs[s][m][z][y_k][x][v  ]) * _kw_12_dv_dv;
*/
        //// Hypervisocisty to stabilize simulation
        const double hypvisc_phi_val = -1.e-5;
        
        const CComplex d4_dx_phi    = (-39. *(Fields[Field::phi][s][m][z][y_k][x+1] - Fields[Field::phi][s][m][z][y_k][x-1])  
                                      + 12. *(Fields[Field::phi][s][m][z][y_k][x+2] - Fields[Field::phi][s][m][z][y_k][x-2]) 
                                      + 56. * Fields[Field::phi][s][m][z][y_k][x  ])/pow4(dx);

        const CComplex hypvisc_phi    = hypvisc_phi_val * ( d4_dx_phi + pow4(ky) * phi_);


        /////////////// Finally the Vlasov equation calculate the time derivatve      //////////////////////
             
        const CComplex dg_dt =
             +  hypvisc_phi
             -  alpha * V[v] * (Island_g + sigma * Island_phi * f0_) +           // Island term
             +  nonLinear[y_k][x][v]                                                  // Non-linear ( array is zero for linear simulations) 
             +  ky* (-(w_n + w_T * (((V[v]*V[v])+ M[m])*kw_T  - sub)) * f0_ * phi_    // Driving term (Temperature/Density gradient)
             -  half_eta_kperp2_phi * f0_)                                            // Contributions from gyro-1 (0 if not neq Gyro-1)
             -  alpha  * V[v]* kp  * ( g + sigma * phi_ * f0_)                        // Linear Landau damping
             +  Coll[s][m][z][y_k][x][v]  ;                                           // Collisional operator
         
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        


        //////////////////////////// Vlasov End ////////////////////////////
        //  time-integrate the distribution function    
        ft [s][m][z][y_k][x][v] = rk[0] * ft[s][m][z][y_k][x][v] + rk[1] * dg_dt             ;
        fss[s][m][z][y_k][x][v] = f1[s][m][z][y_k][x][v]         + (rk[2] * ft[s][m][z][y_k][x][v] + dg_dt) * dt;
      
        }}} }}
   }

}

