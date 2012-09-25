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


  
VlasovIsland::VlasovIsland(Grid *_grid, Parallel *_parallel, Setup *_setup, FileIO *fileIO, Geometry *_geo, FFTSolver *fft, Benchmark *_bench)    
: VlasovAux(_grid, _parallel, _setup, fileIO, _geo, fft, _bench) 
{

    /// Set Magnetic Island structure
    width = setup->get("Island.Width"  , 0.); 
    shear = setup->get("Geometry.Shear", 0.4); 

    //p = { 0.13828847,  0.70216594, -0.01033686 };
    p[0] =  0.13828847;
    p[1] =  0.70216594;
    p[2] = -0.01033686;


}


int VlasovIsland::solve(std::string equation_type, Fields *fields, Array6C f_in, Array6C f_out, double dt, int rk_step, const double rk[3]) 
{

  // do I need both, we can stick to e-m ? Speed penality ?
  if(0) ;  
  
  else if(equation_type == "2D_Island") 
    
      Vlasov_2D_Island((A6zz) f_in.dataZero(), (A6zz) f_out.dataZero(), (A6zz) f0.dataZero(), (A6zz) f.dataZero(), 
                       (A6zz) ft.dataZero(), (A6zz) fields->Field.dataZero(), (A3zz) nonLinearTerms,
                       X, V, M, dt, rk_step, rk);
  
  else   check(-1, DMESG("No Such Equation"));



  return GKC_SUCCESS;

}

void VlasovIsland::Vlasov_2D_Island(
                           CComplex fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           CComplex fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex f0 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           CComplex ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex Fields[Nq][NsLD][NmLD][NzLB][NkyLD][NxLB+4]      ,
                           CComplex nonLinear                  [NkyLD][NxLD  ][NvLD],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           const double dt, const int rk_step, const double rk[3])
{ 

  
    const CComplex imag = ((CComplex) (0. + 1.j));

//    const double _kw_12_dx_dx = 1./(12.*dx*dx);

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


      for(int m=NmLlD; m<=NmLuD; m++) {

       //calculate non-linear term (rk_step == 0 for eigenvalue calculations)
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
             const CComplex phi_ = Fields[Field::phi][s][m][z][y_k][x];

             const CComplex dphi_dx  = (8.*(Fields[Field::phi][s][m][z][y_k][x+1] - Fields[Field::phi][s][m][z][y_k][x-1])
                                         - (Fields[Field::phi][s][m][z][y_k][x+2] - Fields[Field::phi][s][m][z][y_k][x-2]))/(12.*dx)  ;  

        /////////////////////////////////////////////////// Magnetic Island Contribution    /////////////////////////////////////////
        
        // NOTE :  at the Nyquist frequency we have no coupling with higher frequencies (actually phi(m=Ny) = 0. anyway)

        const double zeta      = 2.*M_PI/Ly;
        
        const double xx = pow2(X[x]);
        
        const double psi  = (1. + p[0]*pow(xx,p[1])) * exp( p[2] * xx);
        const double dpsi = (X[x] == 0.) ? 0. :  (p[0] * 2. * X[x] * p[1]*pow(xx, p[1]-1.) + (1. + p[0] * pow(xx,p[1])) * p[2] * 2. * X[x]) * exp(p[2]*xx);


        //const double MagIs     = - 0.5 * w*w*shear/16.  * psi ;//cos(zeta * X[x]);
        //const double dMagIs_dx = - 0.5 * w*w*shear/16. * dpsi;// - sin(zeta * X[x]) * zeta;
        // changed sign , check !!!
        const double MagIs     = - 0.5 * width*width*shear/16.  * psi ;//cos(zeta * X[x]);
        const double dMagIs_dx = - 0.5 * width*width*shear/16. * dpsi;// - sin(zeta * X[x]) * zeta;

        
        const CComplex     phi_p1 = ( y_k == Nky-1) ? 0.                       : Fields[Field::phi][s][m][z][y_k+1][x] ;
        const CComplex     phi_m1 = ( y_k ==  0   ) ? conj(Fields[Field::phi][s][m][z][1][x]) : Fields[Field::phi][s][m][z][y_k-1][x] ;


        const CComplex dphi_dx_p1 = ( y_k == Nky-1) ? 0. 
                                                  : (8.*(Fields[Field::phi][s][m][z][y_k+1][x+1] - Fields[Field::phi][s][m][z][y_k+1][x-1]) - (Fields[Field::phi][s][m][z][y_k+1][x+2] - Fields[Field::phi][s][m][z][y_k+1][x-2]))/(12.*dx)  ;

        const CComplex dphi_dx_m1 = ( y_k ==    0 ) ?  conj(8.*(Fields[Field::phi][s][m][z][    1][x+1] - Fields[Field::phi][s][m][z][    1][x-1]) 
                                                             - (Fields[Field::phi][s][m][z][    1][x+2] - Fields[Field::phi][s][m][z][    1][x-2]))/(12.*dx) 
                                                  :        (8.*(Fields[Field::phi][s][m][z][y_k-1][x+1] - Fields[Field::phi][s][m][z][y_k-1][x-1]) 
                                                             - (Fields[Field::phi][s][m][z][y_k-1][x+2] - Fields[Field::phi][s][m][z][y_k-1][x-2]))/(12.*dx) ;
        
       
        // The magnetic island
    
        /*
                *  For the latter term, the intrigate derivative is the \partial_y * Island
                *  remember the island structure is 
                *  \partial_y (e^{imx} + e^{-imx}) = (i m) * ( e^{imx} - e^{-imx} )
                *
        */
        const CComplex Island_A_phi =   dMagIs_dx * ( ky_m1 * phi_m1 + ky_p1 * phi_p1) - MagIs *  ky_1 * ( dphi_dx_m1 -  dphi_dx_p1);
        
             ///////////////////////////////////////////////////////////////////////////////
            
        const CComplex kp = geo->get_kp(x, ky, z);
           
           CComplex half_eta_kperp2_phi = 0;
           if(isGyro1) { // first order approximation for gyro-kinetics
             const CComplex ddphi_dx_dx = (16. *(Fields[Field::phi][s][m][z][y_k][x+1] + Fields[Field::phi][s][m][z][y_k][x-1]) - (Fields[Field::phi][s][m][z][y_k][x+2] + Fields[Field::phi][s][m][z][y_k][x-2]) - 30.*phi_) * _kw_12_dx_dx;
             half_eta_kperp2_phi     = 0.5 * w_T  * ( (ky*ky) * phi_ + ddphi_dx_dx ) ; 
           }
             
           // Sign has no influence on result ...

               // velocity space magic
        simd_for(int v=NvLlD; v<= NvLuD;v++) {

            const CComplex g    = fs[s][m][z][y_k][x][v];
            const CComplex f0_  = f0[s][m][z][y_k][x][v];


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
             -  alpha * V[v] * (Island_A_F1 + sigma * Island_A_phi * f0_) +           // Island term
             +  nonLinear[y_k][x][v]                                                  // Non-linear ( array is zero for linear simulations) 
             +  ky* (-(w_n + w_T * (((V[v]*V[v])+ M[m])*kw_T  - sub)) * f0_ * phi_    // Driving term (Temperature/Density gradient)
             -  half_eta_kperp2_phi * f0_)                                            // Contributions from gyro-1 (0 if not neq Gyro-1)
             -  alpha  * V[v]* kp  * ( g + sigma * phi_ * f0_)                        // Linear Landau damping
             +  collisionBeta  * (g  + alpha * V[v] * dfs_dv + 2. * ddfs_dvv);        // Lennard-Bernstein Collision term
         
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        


        //////////////////////////// Vlasov End ////////////////////////////
        //  time-integrate the distribution function    
        ft [s][m][z][y_k][x][v] = rk[0] * ft[s][m][z][y_k][x][v] + rk[1] * dg_dt             ;
        fss[s][m][z][y_k][x][v] = f1[s][m][z][y_k][x][v]         + (rk[2] * ft[s][m][z][y_k][x][v] + dg_dt) * dt;
      
        }}} }}
   }

}

