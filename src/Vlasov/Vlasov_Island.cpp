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
  i     = setup->get("Island.Mode"   , 1 ); 
  omega = setup->get("Island.Omega"  , 0.0); 
  shear = setup->get("Geometry.Shear", 0.2); 
    

  //// Setup Magnetic Island structure
  ArrayX = nct::allocate(nct::Range(NxGlB-2, Nx+8))(&MagIs, &dMagIs_dx);
  ArrayY = nct::allocate(nct::Range(NkyLlD, Nky))(&ky_filter);
    
  const double ky0   = setup->get("Island.Filter.ky0", 1.2); 
  const double filter_gradient = setup->get("Island.Filter.Gradient  ", 10.); 
  const double signf = setup->get("Island.Filter.Sign", 0.5); 
    
  for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) {  
      
    ky_filter[y_k] = 0.5 + signf*tanh(filter_gradient*(fft->ky(y_k)-ky0));

  }  

  for(int x = NxGlD-4; x <= NxGuD+4; x++) {  

    const double zeta      = 2.*M_PI/Ly;
       
    // Island cos(k_T y) = 0.5 [ exp(-i m) + exp( i m ) ]
    // The poloidal island term is direclty included per 
    // mode-mode coupling
    
    // we used simple fitting model to get least square coefficient
    const double p[] = { 0.13828847,  0.70216594, -0.01033686 };

    const double xx   = pow2(X[x]);
    const double psi  = (1. + p[0]*pow(xx,p[1])) * exp( p[2] * xx);
    const double dpsi = (X[x] == 0.) ? 0. :  (p[0] * 2. * X[x] * p[1]*pow(xx, p[1]-1.) + (1. + p[0] * pow(xx,p[1])) * p[2] * 2. * X[x]) * exp(p[2]*xx);

    // Corresponding fit function to translate from sparatrix island width to
    
   
    // using Horner's rule for polynomial evaluation
    auto evalPoly = [=](int N, const double *polyConst, double x) -> double 
    {
      double res = 0.;
      for(int n=0; n < N-1; n++) res += (res + polyConst[n])*x;
      return res + polyConst[N-1];

    };
    
    const double polyIsland[] = 
    { -6.28828806e-13, 1.60965096e-10, -1.77107115e-08, 1.09512432e-06,
      -4.17566040e-05, 1.01434574e-03, -1.56855609e-02, 1.50685821e-01,
      -8.66414053e-01, 3.18134960e+00,  6.27800358e-01 };

    const double width_scale = (width == 0.) ? 0. : evalPoly(11, polyIsland, abs(width));

    auto sign = [=](double a) { return a >= 0. ? 1. : -1.;} ;

    // we need to translate island width size to translation size. 
     MagIs   [x] = 0.5 * sign(width) * pow2(width_scale) * shear/16. * psi ;//cos(zeta * X[x]);
    dMagIs_dx[x] = 0.5 * sign(width) * pow2(width_scale) * shear/16. * dpsi;// - sin(zeta * X[x]) * zeta;
    
    //MagIs[x]     = 0.5 * width*width*shear/16.  * cos(zeta * X[x]);
    //dMagIs_dx[x] = 0.5 * width*width*shear/16. *  (-sin(zeta * X[x])) * zeta;
 
  }
        
   
  // if electro-magnetic version is used intialize fields
  if(Nq >= 2) {
  
    ArrayAp_mod = nct::allocate(grid->RzLB, grid->RkyLD, nct::Range(NxLlD-4, NxLD+8))(&Ap_mod);

    ArrayXi(&Xi_lin);
    ArrayG(&G_lin );

    [=](CComplex Ap_mod[NzLB][Nky][NxLB+4]) {

      for(int z = NzLlD  ; z <= NzLuD  ; z++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { 
      for(int x = NxLlB-2; x <= NxLuB+2; x++) {

       Ap_mod[z][y_k][x] = - ((y_k == i ) ? 1. : 0.) * MagIs[x];

      } } }
      }( (A3zz) Ap_mod);
  }

  initDataOutput(setup, fileIO);

}


void VlasovIsland::solve(std::string equation_type, Fields *fields, CComplex *f_in, CComplex *f_out, double dt, int rk_step, const double rk[3]) 
{
  // do I need both, we can stick to e-m ? Speed penality ?
  if(0) ;  
  else if(equation_type == "2D_Island") 
    
      Vlasov_2D_Island((A6zz) f_in, (A6zz) f_out, (A6zz) f0, (A6zz) f, 
                       (A6zz) ft  , (A6zz) Coll, (A6zz) fields->Field, (A3zz) nonLinearTerm,
                       MagIs, dMagIs_dx, X, V, M, (A3zz) Ap_mod, (A4zz) fields->Field0, dt, rk_step, rk);
  
  else if(equation_type == "2D_Island_EM") 
    
      Vlasov_2D_Island_EM   ((A6zz) f_in, (A6zz) f_out, (A6zz) f0, (A6zz) f,
                   (A6zz) ft, (A6zz) Coll, (A6zz) fields->Field, (A3zz) nonLinearTerm,
                   (A4zz) Xi, (A4zz) G, (A4zz) Xi_lin, (A4zz) G_lin, (A3zz) Ap_mod, (A4zz) fields->Field0, dt, rk_step, rk);

  else if(equation_type == "2D_Island_Filter") 
    
      Vlasov_2D_Island_filter((A6zz) f_in, (A6zz) f_out, (A6zz) f0, (A6zz) f, 
                       (A6zz) ft  , (A6zz) Coll, (A6zz) fields->Field, (A3zz) nonLinearTerm,
                       MagIs, dMagIs_dx, X, V, M, dt, rk_step, rk);

  else if(equation_type == "2D_Island_Equi")

      Vlasov_2D_Island_Equi((A6zz) f_in, (A6zz) f_out, (A6zz) f0, (A6zz) f, 
                       (A6zz) ft  , (A6zz) Coll, (A6zz) fields->Field, (A3zz) nonLinearTerm,
                       X, V, M, dt, rk_step, rk);

  else if(equation_type == "2D_Island_Rotation") 
    
      Vlasov_2D_Island_Rotation((A6zz) f_in, (A6zz) f_out, (A6zz) f0, (A6zz) f, 
                       (A6zz) ft  , (A6zz) Coll, (A6zz) fields->Field, (A3zz) nonLinearTerm,
                       MagIs, dMagIs_dx, X, V, M, dt, rk_step, rk);
    
  else   check(-1, DMESG("No Such Equation"));


  return;

}

void VlasovIsland::Vlasov_2D_Island(
                           CComplex fs        [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex fss       [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f0  [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f1  [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex ft        [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Coll[NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Fields[Nq][NsLD][NmLD][NzLB][Nky][NxLB+4]      ,
                           CComplex nonLinearTerm                  [Nky][NxLD  ][NvLD],
                           const double MagIs[NxGB], const double dMagIs_dx[NxGB], 
                           const double X[NxGB+4], const double V[NvGB], const double M[NmGB],
                           const CComplex Ap_mod                    [NzLB][Nky][NxLB+4]      ,
                           CComplex Field0[Nq][NzLD][Nky][NxLD]   ,
                           const double dt, const int rk_step, const double rk[3])
{ 

  if((Nq > 1)) {
   Field0[Field::Ap][NzLlD:NzLD][:][NxLlD:NxLD] = 0. ;//Ap_mod[NzLlD:NzLD][:][NxLlD:NxLD];
  }

  for(int s = NsLlD; s <= NsLuD; s++) {
        
      // small abbrevations
      const double w_n   = species[s].w_n;
      const double w_T   = species[s].w_T;
      const double alpha = species[s].alpha;
      const double sigma = species[s].sigma;
      const double Temp  = species[s].T0;
      const double sub   = (species[s].doGyro) ? 3./2. : 1./2.;

      const double kw_T = 1./Temp;

      bool isGyro1 = (species[s].gyroModel == "Gyro-1");
      
      const double rho_t2 = species[s].T0 * species[s].m / (pow2(species[s].q) * plasma->B0); 


      for(int m=NmLlD; m<=NmLuD; m++) { for(int z=NzLlD; z<= NzLuD;z++) {  
        
        //calculate non-linear term (rk_step == 0 for eigenvalue calculations)
        if(doNonLinear && (rk_step != 0)) calculatePoissonBracket(nullptr, nullptr, fs, Fields, z, m, s, nonLinearTerm, Xi_max, false); 
        
        #pragma omp for
        for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) {

        
        // Note : for negative modes we need to use complex conjugate value
        const CComplex ky     = _imag * fft->ky(y_k);
        
        // We need to take care of boundaries. For poloidal numbers y_k > N_k-1, we  use zero.
        // For y_k < 0, the corresponding complex conjugate value is used.
        // is ky_p1 = ky or not ?!
        const CComplex ky_1   = _imag * fft->ky(i);
        const CComplex ky_p1  = ((y_k+i) >= Nky-1) ? 0.    : _imag * fft->ky(y_k+i);
        const CComplex ky_m1  = _imag * fft->ky(y_k-i); 
        
        //const CComplex ky_1   = _imag * fft->ky(1);
        //const CComplex ky_m1  = _imag * (y_k == 0    ) ? -fft-ky_1 : _imag * fft->ky(y_k-i); 
        //const CComplex ky_p1  = ((y_k+i) >= Nky-1) ? 0.    : _imag * fft->ky(y_k+i);
        //const CComplex ky_m1  = _imag * fft->ky(y_k-i);

        
        for(int x = NxLlD; x <= NxLuD; x++) {  

          
          const CComplex phi_ = Fields[Field::phi][s][m][z][y_k][x];
          
          const CComplex dphi_dx  = (8.*(Fields[Field::phi][s][m][z][y_k][x+1] - Fields[Field::phi][s][m][z][y_k][x-1])
                                      - (Fields[Field::phi][s][m][z][y_k][x+2] - Fields[Field::phi][s][m][z][y_k][x-2])) * _kw_12_dx  ;  

          /////////////////////////////////////////////////// Magnetic Island Contribution    /////////////////////////////////////////
        
          // NOTE :  at the Nyquist frequency we have no coupling with higher frequencies (actually phi(m=Ny) = 0. anyway)

          const CComplex     phi_p1 = ((y_k+i) >= Nky-1) ? 0.                                          : Fields[Field::phi][s][m][z][y_k+i][x] ;
          const CComplex     phi_m1 = ((y_k-i) <  0    ) ? conj(Fields[Field::phi][s][m][z][i-y_k][x]) : Fields[Field::phi][s][m][z][y_k-i][x] ;
          //const CComplex     phi_m1 = ( y_k ==  0   ) ? (Fields[Field::phi][s][m][z][1][x]) : Fields[Field::phi][s][m][z][y_k-1][x] ;

          
          // X-derivative (First Derivative with Central Difference 4th) of phi for poloidal mode +1, take care of Nyquist frequency
          const CComplex dphi_dx_p1 = ( (y_k+i) >= Nky-1) 
                                   ? 0.

                                   :      (8.*(Fields[Field::phi][s][m][z][y_k+i][x+1] - Fields[Field::phi][s][m][z][y_k+i][x-1]) 
                                            - (Fields[Field::phi][s][m][z][y_k+i][x+2] - Fields[Field::phi][s][m][z][y_k+i][x-2])) * _kw_12_dx  ;

          // X-derivative (1st deriv. CD-4 )of phi for poloidal mode -1, take care of complex conjugate relation for y_k=-1
          const CComplex dphi_dx_m1 = ( (y_k-i) < 0 ) 
                                   ?  conj(8.*(Fields[Field::phi][s][m][z][i-y_k][x+1] - Fields[Field::phi][s][m][z][i-y_k][x-1]) 
                                            - (Fields[Field::phi][s][m][z][i-y_k][x+2] - Fields[Field::phi][s][m][z][i-y_k][x-2])) * _kw_12_dx
                                   :      (8.*(Fields[Field::phi][s][m][z][y_k-i][x+1] - Fields[Field::phi][s][m][z][y_k-i][x-1]) 
                                            - (Fields[Field::phi][s][m][z][y_k-i][x+2] - Fields[Field::phi][s][m][z][y_k-i][x-2])) * _kw_12_dx ;
        
       
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
        const CComplex dfs_dx_p1  =  ((y_k+i) >= Nky-i)

                            ? 0.

                            : (8. *(fs[s][m][z][y_k+i][x+1][v] - fs[s][m][z][y_k+i][x-1][v])  
                                 - (fs[s][m][z][y_k+i][x+2][v] - fs[s][m][z][y_k+i][x-2][v])) * _kw_12_dx;


        // X-derivative of f1 (1-CD4) for poloidal mode -1, take care of complex conjugate relation for y_k=-1 
        const CComplex dfs_dx_m1  =  ( (y_k-i) < 0  )

                        ? conj(8. *(fs[s][m][z][i-y_k][x+1][v] - fs[s][m][z][i-y_k][x-1][v])  
                                 - (fs[s][m][z][i-y_k][x+2][v] - fs[s][m][z][i-y_k][x-2][v])) * _kw_12_dx 

                        :     (8. *(fs[s][m][z][y_k-i][x+1][v] - fs[s][m][z][y_k-i][x-1][v])  
                                 - (fs[s][m][z][y_k-i][x+2][v] - fs[s][m][z][y_k-i][x-2][v])) * _kw_12_dx;

        // Note Nky-1 is the maximum mode number Nky = 6 i-> [ 0, 1, 2, 3, 4, 5] 
        const CComplex fs_p1      = ((y_k+i) >= Nky-1) ? 0.                             : fs[s][m][z][y_k+i][x][v] ;
        const CComplex fs_m1      = ((y_k-i) <  0    ) ? conj(fs[s][m][z][i-y_k][x][v]) : fs[s][m][z][y_k-i][x][v] ;
        //const CComplex fs_m1      = (y_k ==  0   ) ? (fs[s][m][z][1][x][v]) : fs[s][m][z][y_k-1][x][v] ;
         
        // Coupling of phase-space with Island mode-mode coupling
        const CComplex Island_g =  dMagIs_dx[x] * (ky_m1 * fs_m1  + ky_p1 * fs_p1 )  -  MagIs[x]  * ky_1 *  (dfs_dx_m1  - dfs_dx_p1 )  ;


        //// Hypervisocisty to stabilize simulation
        /*
        const double hypvisc_phi_val = -1.e-5;
        
        const CComplex d4_dx_phi    = (-39. *(Fields[Field::phi][s][m][z][y_k][x+1] - Fields[Field::phi][s][m][z][y_k][x-1])  
                                      + 12. *(Fields[Field::phi][s][m][z][y_k][x+2] - Fields[Field::phi][s][m][z][y_k][x-2]) 
                                      + 56. * Fields[Field::phi][s][m][z][y_k][x  ])/pow4(dx);

        const CComplex hypvisc_phi    = hypvisc_phi_val * ( d4_dx_phi + pow4(ky) * phi_);
*/

        /////////////// Finally the Vlasov equation calculate the time derivatve      //////////////////////
             
        const CComplex dg_dt =
  //           +  hypvisc_phi
             -  alpha * V[v] * (Island_g + sigma * Island_phi * f0_) +                // Island term
             +  nonLinearTerm[y_k][x][v]                                             // Non-linear ( array is zero for linear simulations) 
             +  ky* (-(w_n + w_T * (((V[v]*V[v])+ M[m])*kw_T  - sub)) * f0_ * 
//                 phi_    )                     // Source term (Temperature/Density gradient)
//                (phi_ - (y_k == 1 ? V[v] * MagIs[x] : 0.  )  ) )                       // Source term (Temperature/Density gradient)
             phi_)                        // Source term (Temperature/Density gradient)
     //        -  half_eta_kperp2_phi * f0_)                                            // Contributions from gyro-1 (0 if not neq Gyro-1)
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


void VlasovIsland::initDataOutput(Setup *setup, FileIO *fileIO) 
{
    
   hid_t islandGroup = check(H5Gcreate(fileIO->getFileID(), "/Islands",H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), DMESG("H5Gcreate"));
          
   // Length scale
   check(H5LTset_attribute_double(islandGroup, ".", "MagIs"    , &MagIs[NxGlD]    , Nx), DMESG("Attribute"));
   check(H5LTset_attribute_double(islandGroup, ".", "dMagIs_dx", &dMagIs_dx[NxGlD], Nx), DMESG("Attribute"));
   
   check(H5LTset_attribute_double(islandGroup, ".", "Width", &width, 1), DMESG("Attribute"));
   check(H5LTset_attribute_double(islandGroup, ".", "Shear", &shear, 1), DMESG("Attribute"));
   check(H5LTset_attribute_double(islandGroup, ".", "Omega", &omega, 1), DMESG("Attribute"));
   check(H5LTset_attribute_int   (islandGroup, ".", "Mode" , &i    , 1), DMESG("Attribute"));
         
    H5Gclose(islandGroup);

          

}   

void VlasovIsland::printOn(std::ostream &output) const 
{

         output   << " Island    |  Width : " << width << std::endl;
}


void VlasovIsland::Vlasov_2D_Island_Equi(
                           CComplex fs       [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex fss      [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f0 [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f1 [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex ft       [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Coll      [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Fields[Nq][NsLD][NmLD][NzLB][Nky][NxLB+4]      ,
                           CComplex nonLinearTerm                  [Nky][NxLD  ][NvLD],
                           const double X[NxGB+4], const double V[NvGB], const double M[NmGB],
                           const double dt, const int rk_step, const double rk[3])
{ 



    for(int s = NsLlD; s <= NsLuD; s++) {
        
      // small abbrevations
      const double w_n   = species[s].w_n;
      const double w_T   = species[s].w_T;
      const double alpha = species[s].alpha;
      const double sigma = species[s].sigma;
      const double Temp  = species[s].T0;
      const double sub   = (species[s].doGyro) ? 3./2. : 1./2.;

      const double kw_T = 1./Temp;

      bool isGyro1 = (species[s].gyroModel == "Gyro-1");
      
      const double rho_t2 = species[s].T0 * species[s].m / (pow2(species[s].q) * plasma->B0); 


      for(int m=NmLlD; m<=NmLuD; m++) { for(int z=NzLlD; z<= NzLuD;z++) {  
        
        //calculate non-linear term (rk_step == 0 for eigenvalue calculations)
        if(doNonLinear && (rk_step != 0)) calculatePoissonBracket(nullptr, nullptr, fs, Fields, z, m, s, nonLinearTerm, Xi_max, false); 
        
        for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) {

        
        // Note : for negative modes we need to use complex conjugate value
        const CComplex ky     = _imag * fft->ky(y_k) ;
        
        

        
        for(int x=NxLlD; x<= NxLuD;x++) {  

          
          const CComplex phi_ = Fields[Field::phi][s][m][z][y_k][x];
          
          const CComplex dphi_dx  = (8.*(Fields[Field::phi][s][m][z][y_k][x+1] - Fields[Field::phi][s][m][z][y_k][x-1])
                                      - (Fields[Field::phi][s][m][z][y_k][x+2] - Fields[Field::phi][s][m][z][y_k][x-2])) * _kw_12_dx  ;  

        //  std::cout << "------1----> " << dMagIs_dx[x] << std::endl;  
        const CComplex kp = geo->get_kp(x, ky, z) - dMagIs_dx[x] * ky;
        //  std::cout << "------2----> " << kp << std::endl;  
           
        CComplex half_eta_kperp2_phi = 0;

        
        // velocity space magic
        simd_for(int v=NvLlD; v<= NvLuD;v++) {

            const CComplex g    = fs[s][m][z][y_k][x][v];
            const CComplex f0_  = f0[s][m][z][y_k][x][v];


        /////////////////////////////////////////////////// Magnetic Island Contribution    /////////////////////////////////////////
      

        /////////// Collisions ////////////////////////////////////////////////////////////////////

        const CComplex dfs_dv    = (8.  *(fs[s][m][z][y_k][x][v+1] - fs[s][m][z][y_k][x][v-1]) - 
                                         (fs[s][m][z][y_k][x][v+2] - fs[s][m][z][y_k][x][v-2])) * _kw_12_dv;

        const CComplex ddfs_dvv  = (16. *(fs[s][m][z][y_k][x][v+1] + fs[s][m][z][y_k][x][v-1]) -
                                         (fs[s][m][z][y_k][x][v+2] + fs[s][m][z][y_k][x][v-2]) 
                                    - 30.*fs[s][m][z][y_k][x][v  ]) * _kw_12_dv_dv;

        //// Hypervisocisty to stabilize simulation
        /*
        const double hypvisc_phi_val = -1.e-5;
        
        const CComplex d4_dx_phi    = (-39. *(Fields[Field::phi][s][m][z][y_k][x+1] - Fields[Field::phi][s][m][z][y_k][x-1])  
                                      + 12. *(Fields[Field::phi][s][m][z][y_k][x+2] - Fields[Field::phi][s][m][z][y_k][x-2]) 
                                      + 56. * Fields[Field::phi][s][m][z][y_k][x  ])/pow4(dx);

        const CComplex hypvisc_phi    = hypvisc_phi_val * ( d4_dx_phi + pow4(ky) * phi_);
*/

        /////////////// Finally the Vlasov equation calculate the time derivatve      //////////////////////
             
        const CComplex dg_dt =
             +  nonLinearTerm[y_k][x][v]                                             // Non-linear ( array is zero for linear simulations) 
             +  ky* (-(w_n + w_T * (((V[v]*V[v])+ M[m])*kw_T  - sub)) * f0_ * phi_     // Driving term (Temperature/Density gradient)
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
//  std::cout << "OUT" << std::flush;  

}


void VlasovIsland::Vlasov_2D_Island_filter(
                           CComplex fs        [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex fss       [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f0  [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f1  [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex ft        [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Coll[NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Fields[Nq][NsLD][NmLD][NzLB][Nky][NxLB+4]      ,
                           CComplex nonLinearTerm                  [Nky][NxLD  ][NvLD],
                           const double MagIs[NxGB], const double dMagIs_dx[NxGB], 
                           const double X[NxGB+4], const double V[NvGB], const double M[NmGB],
                           const double dt, const int rk_step, const double rk[3])
{ 

  
    for(int s = NsLlD; s <= NsLuD; s++) {
        
      // small abbrevations
      const double w_n   = species[s].w_n;
      const double w_T   = species[s].w_T;
      const double alpha = species[s].alpha;
      const double sigma = species[s].sigma;
      const double Temp  = species[s].T0;
      const double sub   = (species[s].doGyro) ? 3./2. : 1./2.;

      const double kw_T = 1./Temp;

      bool isGyro1 = (species[s].gyroModel == "Gyro-1");
      
      const double rho_t2 = species[s].T0 * species[s].m / (pow2(species[s].q) * plasma->B0); 


      for(int m=NmLlD; m<=NmLuD; m++) { for(int z=NzLlD; z<= NzLuD;z++) {  
        
        //calculate non-linear term (rk_step == 0 for eigenvalue calculations)
        if(doNonLinear && (rk_step != 0)) calculatePoissonBracket(nullptr, nullptr, fs, Fields, z, m, s, nonLinearTerm, Xi_max, false); 
        
        for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) {

        
        // Note : for negative modes we need to use complex conjugate value
        const CComplex ky     = _imag * fft->ky(y_k);
        
        // We need to take care of boundaries. For poloidal numbers y_k > N_k-1, we  use zero.
        // For y_k < 0, the corresponding complex conjugate value is used.
        // is ky_p1 = ky or not ?!
        const CComplex ky_1   = _imag * fft->ky(i);
        const CComplex ky_p1  = ((y_k+i) >= Nky-1) ? 0.    : _imag * fft->ky(y_k+i);
        const CComplex ky_m1  = _imag * fft->ky(y_k-i); 
        
        //const CComplex ky_1   = _imag * fft->ky(1);
        //const CComplex ky_m1  = _imag * (y_k == 0    ) ? -fft-ky_1 : _imag * fft->ky(y_k-i); 
        //const CComplex ky_p1  = ((y_k+i) >= Nky-1) ? 0.    : _imag * fft->ky(y_k+i);
        //const CComplex ky_m1  = _imag * fft->ky(y_k-i);

        
        for(int x=NxLlD; x<= NxLuD;x++) {  

          
          const CComplex phi_ = Fields[Field::phi][s][m][z][y_k][x];
          
          const CComplex dphi_dx  = (8.*(Fields[Field::phi][s][m][z][y_k][x+1] - Fields[Field::phi][s][m][z][y_k][x-1])
                                      - (Fields[Field::phi][s][m][z][y_k][x+2] - Fields[Field::phi][s][m][z][y_k][x-2])) * _kw_12_dx  ;  

          /////////////////////////////////////////////////// Magnetic Island Contribution    /////////////////////////////////////////
        
          // NOTE :  at the Nyquist frequency we have no coupling with higher frequencies (actually phi(m=Ny) = 0. anyway)

          const CComplex     phi_p1 = ((y_k+i) >= Nky-1) ? 0.                                          : Fields[Field::phi][s][m][z][y_k+i][x] ;
          const CComplex     phi_m1 = ((y_k-i) <  0    ) ? conj(Fields[Field::phi][s][m][z][i-y_k][x]) : Fields[Field::phi][s][m][z][y_k-i][x] ;
          //const CComplex     phi_m1 = ( y_k ==  0   ) ? (Fields[Field::phi][s][m][z][1][x]) : Fields[Field::phi][s][m][z][y_k-1][x] ;

          
          // X-derivative (First Derivative with Central Difference 4th) of phi for poloidal mode +1, take care of Nyquist frequency
          const CComplex dphi_dx_p1 = ( (y_k+i) >= Nky-1) 
                                   ? 0.

                                   :      (8.*(Fields[Field::phi][s][m][z][y_k+i][x+1] - Fields[Field::phi][s][m][z][y_k+i][x-1]) 
                                            - (Fields[Field::phi][s][m][z][y_k+i][x+2] - Fields[Field::phi][s][m][z][y_k+i][x-2])) * _kw_12_dx  ;

          // X-derivative (1st deriv. CD-4 )of phi for poloidal mode -1, take care of complex conjugate relation for y_k=-1
          const CComplex dphi_dx_m1 = ( (y_k-i) < 0 ) 
                                   ?  conj(8.*(Fields[Field::phi][s][m][z][i-y_k][x+1] - Fields[Field::phi][s][m][z][i-y_k][x-1]) 
                                            - (Fields[Field::phi][s][m][z][i-y_k][x+2] - Fields[Field::phi][s][m][z][i-y_k][x-2])) * _kw_12_dx
                                   :      (8.*(Fields[Field::phi][s][m][z][y_k-i][x+1] - Fields[Field::phi][s][m][z][y_k-i][x-1]) 
                                            - (Fields[Field::phi][s][m][z][y_k-i][x+2] - Fields[Field::phi][s][m][z][y_k-i][x-2])) * _kw_12_dx ;
        
       
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
        const CComplex dfs_dx_p1  =  ((y_k+i) >= Nky-i)

                            ? 0.

                            : (8. *(fs[s][m][z][y_k+i][x+1][v] - fs[s][m][z][y_k+i][x-1][v])  
                                 - (fs[s][m][z][y_k+i][x+2][v] - fs[s][m][z][y_k+i][x-2][v])) * _kw_12_dx;


        // X-derivative of f1 (1-CD4) for poloidal mode -1, take care of complex conjugate relation for y_k=-1 
        const CComplex dfs_dx_m1  =  ( (y_k-i) < 0  )

                        ? conj(8. *(fs[s][m][z][i-y_k][x+1][v] - fs[s][m][z][i-y_k][x-1][v])  
                                 - (fs[s][m][z][i-y_k][x+2][v] - fs[s][m][z][i-y_k][x-2][v])) * _kw_12_dx 

                        :     (8. *(fs[s][m][z][y_k-i][x+1][v] - fs[s][m][z][y_k-i][x-1][v])  
                                 - (fs[s][m][z][y_k-i][x+2][v] - fs[s][m][z][y_k-i][x-2][v])) * _kw_12_dx;

        // Note Nky-1 is the maximum mode number Nky = 6 i-> [ 0, 1, 2, 3, 4, 5] 
        const CComplex fs_p1      = ((y_k+i) >= Nky-1) ? 0.                             : fs[s][m][z][y_k+i][x][v] ;
        const CComplex fs_m1      = ((y_k-i) <  0    ) ? conj(fs[s][m][z][i-y_k][x][v]) : fs[s][m][z][y_k-i][x][v] ;
        //const CComplex fs_m1      = (y_k ==  0   ) ? (fs[s][m][z][1][x][v]) : fs[s][m][z][y_k-1][x][v] ;
         
        // Coupling of phase-space with Island mode-mode coupling
        const CComplex Island_g =  dMagIs_dx[x] * (ky_m1 * fs_m1  + ky_p1 * fs_p1 )  -  MagIs[x]  * ky_1 *  (dfs_dx_m1  - dfs_dx_p1 )  ;

        
        /////////// Collisions ////////////////////////////////////////////////////////////////////

        const CComplex dfs_dv    = (8.  *(fs[s][m][z][y_k][x][v+1] - fs[s][m][z][y_k][x][v-1]) - 
                                         (fs[s][m][z][y_k][x][v+2] - fs[s][m][z][y_k][x][v-2])) * _kw_12_dv;

        const CComplex ddfs_dvv  = (16. *(fs[s][m][z][y_k][x][v+1] + fs[s][m][z][y_k][x][v-1]) -
                                         (fs[s][m][z][y_k][x][v+2] + fs[s][m][z][y_k][x][v-2]) 
                                    - 30.*fs[s][m][z][y_k][x][v  ]) * _kw_12_dv_dv;

        //// Hypervisocisty to stabilize simulation
        /*
        const double hypvisc_phi_val = -1.e-5;
        
        const CComplex d4_dx_phi    = (-39. *(Fields[Field::phi][s][m][z][y_k][x+1] - Fields[Field::phi][s][m][z][y_k][x-1])  
                                      + 12. *(Fields[Field::phi][s][m][z][y_k][x+2] - Fields[Field::phi][s][m][z][y_k][x-2]) 
                                      + 56. * Fields[Field::phi][s][m][z][y_k][x  ])/pow4(dx);

        const CComplex hypvisc_phi    = hypvisc_phi_val * ( d4_dx_phi + pow4(ky) * phi_);
*/

        /////////////// Finally the Vlasov equation calculate the time derivatve      //////////////////////
             
        const CComplex dg_dt =
          ky_filter[y_k] * (                                                          // filter artifically remove parts
             -  alpha * V[v] * (Island_g + sigma * Island_phi * f0_) +                // Island term
             +  nonLinearTerm[y_k][x][v]                                             // Non-linear ( array is zero for linear simulations) 
             +  ky* (-(w_n + w_T * (((V[v]*V[v])+ M[m])*kw_T  - sub)) * f0_ * phi_    // Driving term (Temperature/Density gradient)
             -  half_eta_kperp2_phi * f0_)                                            // Contributions from gyro-1 (0 if not neq Gyro-1)
             -  alpha  * V[v]* kp  * ( g + sigma * phi_ * f0_)                        // Linear Landau damping
             +  Coll[s][m][z][y_k][x][v]                                             // Collisional operator
             );
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        


        //////////////////////////// Vlasov End ////////////////////////////
        //  time-integrate the distribution function    
        ft [s][m][z][y_k][x][v] = rk[0] * ft[s][m][z][y_k][x][v] + rk[1] * dg_dt             ;
        fss[s][m][z][y_k][x][v] = f1[s][m][z][y_k][x][v]         + (rk[2] * ft[s][m][z][y_k][x][v] + dg_dt) * dt;
      
        }}} }}
   }

}


void VlasovIsland::Vlasov_2D_Island_Rotation(
                           CComplex fs        [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex fss       [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f0  [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f1  [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex ft        [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Coll[NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Fields[Nq][NsLD][NmLD][NzLB][Nky][NxLB+4]      ,
                           CComplex nonLinearTerm                  [Nky][NxLD  ][NvLD],
                           const double MagIs[NxGB], const double dMagIs_dx[NxGB], 
                           const double X[NxGB+4], const double V[NvGB], const double M[NmGB],
                           const double dt, const int rk_step, const double rk[3])
{ 

    static double Time = 0.;  if(rk_step == 1) Time += dt;

    const CComplex IslandPhase = cexp(_imag * Time * omega);


    for(int s = NsLlD; s <= NsLuD; s++) {
        
      // small abbrevations
      const double w_n   = species[s].w_n;
      const double w_T   = species[s].w_T;
      const double alpha = species[s].alpha;
      const double sigma = species[s].sigma;
      const double Temp  = species[s].T0;
      const double sub   = (species[s].doGyro) ? 3./2. : 1./2.;

      const double kw_T = 1./Temp;

      bool isGyro1 = (species[s].gyroModel == "Gyro-1");
      
      const double rho_t2 = species[s].T0 * species[s].m / (pow2(species[s].q) * plasma->B0); 


      for(int m=NmLlD; m<=NmLuD; m++) { for(int z=NzLlD; z<= NzLuD;z++) {  
        
        //calculate non-linear term (rk_step == 0 for eigenvalue calculations)
        if(doNonLinear && (rk_step != 0)) calculatePoissonBracket(nullptr, nullptr, fs, Fields, z, m, s, nonLinearTerm, Xi_max, false); 
        
        #pragma omp for
        for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) {

        
        // Note : for negative modes we need to use complex conjugate value
        const CComplex ky     = _imag * fft->ky(y_k);
        
        // We need to take care of boundaries. For poloidal numbers y_k > N_k-1, we  use zero.
        // For y_k < 0, the corresponding complex conjugate value is used.
        // is ky_p1 = ky or not ?!
        const CComplex ky_1   = _imag * fft->ky(i);
        const CComplex ky_p1  = ((y_k+i) >= Nky-1) ? 0.    : _imag * fft->ky(y_k+i);
        const CComplex ky_m1  = _imag * fft->ky(y_k-i); 
        
        //const CComplex ky_1   = _imag * fft->ky(1);
        //const CComplex ky_m1  = _imag * (y_k == 0    ) ? -fft-ky_1 : _imag * fft->ky(y_k-i); 
        //const CComplex ky_p1  = ((y_k+i) >= Nky-1) ? 0.    : _imag * fft->ky(y_k+i);
        //const CComplex ky_m1  = _imag * fft->ky(y_k-i);

        
        for(int x=NxLlD; x<= NxLuD;x++) {  

          
          const CComplex phi_ = Fields[Field::phi][s][m][z][y_k][x];
          
          const CComplex dphi_dx  = (8.*(Fields[Field::phi][s][m][z][y_k][x+1] - Fields[Field::phi][s][m][z][y_k][x-1])
                                      - (Fields[Field::phi][s][m][z][y_k][x+2] - Fields[Field::phi][s][m][z][y_k][x-2])) * _kw_12_dx  ;  

          /////////////////////////////////////////////////// Magnetic Island Contribution    /////////////////////////////////////////
        
          // NOTE :  at the Nyquist frequency we have no coupling with higher frequencies (actually phi(m=Ny) = 0. anyway)

          const CComplex     phi_p1 = ((y_k+i) >= Nky-1) ? 0.                                          : Fields[Field::phi][s][m][z][y_k+i][x] ;
          const CComplex     phi_m1 = ((y_k-i) <  0    ) ? conj(Fields[Field::phi][s][m][z][i-y_k][x]) : Fields[Field::phi][s][m][z][y_k-i][x] ;
          //const CComplex     phi_m1 = ( y_k ==  0   ) ? (Fields[Field::phi][s][m][z][1][x]) : Fields[Field::phi][s][m][z][y_k-1][x] ;

          
          // X-derivative (First Derivative with Central Difference 4th) of phi for poloidal mode +1, take care of Nyquist frequency
          const CComplex dphi_dx_p1 = ( (y_k+i) >= Nky-1) 
                                   ? 0.

                                   :      (8.*(Fields[Field::phi][s][m][z][y_k+i][x+1] - Fields[Field::phi][s][m][z][y_k+i][x-1]) 
                                            - (Fields[Field::phi][s][m][z][y_k+i][x+2] - Fields[Field::phi][s][m][z][y_k+i][x-2])) * _kw_12_dx  ;

          // X-derivative (1st deriv. CD-4 )of phi for poloidal mode -1, take care of complex conjugate relation for y_k=-1
          const CComplex dphi_dx_m1 = ( (y_k-i) < 0 ) 
                                   ?  conj(8.*(Fields[Field::phi][s][m][z][i-y_k][x+1] - Fields[Field::phi][s][m][z][i-y_k][x-1]) 
                                            - (Fields[Field::phi][s][m][z][i-y_k][x+2] - Fields[Field::phi][s][m][z][i-y_k][x-2])) * _kw_12_dx
                                   :      (8.*(Fields[Field::phi][s][m][z][y_k-i][x+1] - Fields[Field::phi][s][m][z][y_k-i][x-1]) 
                                            - (Fields[Field::phi][s][m][z][y_k-i][x+2] - Fields[Field::phi][s][m][z][y_k-i][x-2])) * _kw_12_dx ;
        
       
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
        const CComplex dfs_dx_p1  =  ((y_k+i) >= Nky-i)

                            ? 0.

                            : (8. *(fs[s][m][z][y_k+i][x+1][v] - fs[s][m][z][y_k+i][x-1][v])  
                                 - (fs[s][m][z][y_k+i][x+2][v] - fs[s][m][z][y_k+i][x-2][v])) * _kw_12_dx;


        // X-derivative of f1 (1-CD4) for poloidal mode -1, take care of complex conjugate relation for y_k=-1 
        const CComplex dfs_dx_m1  =  ( (y_k-i) < 0  )

                        ? conj(8. *(fs[s][m][z][i-y_k][x+1][v] - fs[s][m][z][i-y_k][x-1][v])  
                                 - (fs[s][m][z][i-y_k][x+2][v] - fs[s][m][z][i-y_k][x-2][v])) * _kw_12_dx 

                        :     (8. *(fs[s][m][z][y_k-i][x+1][v] - fs[s][m][z][y_k-i][x-1][v])  
                                 - (fs[s][m][z][y_k-i][x+2][v] - fs[s][m][z][y_k-i][x-2][v])) * _kw_12_dx;

        // Note Nky-1 is the maximum mode number Nky = 6 i-> [ 0, 1, 2, 3, 4, 5] 
        const CComplex fs_p1      = ((y_k+i) >= Nky-1) ? 0.                             : fs[s][m][z][y_k+i][x][v] ;
        const CComplex fs_m1      = ((y_k-i) <  0    ) ? conj(fs[s][m][z][i-y_k][x][v]) : fs[s][m][z][y_k-i][x][v] ;
        //const CComplex fs_m1      = (y_k ==  0   ) ? (fs[s][m][z][1][x][v]) : fs[s][m][z][y_k-1][x][v] ;
         
        // Coupling of phase-space with Island mode-mode coupling
        const CComplex Island_g =  dMagIs_dx[x] * (ky_m1 * fs_m1  + ky_p1 * fs_p1 )  -  MagIs[x]  * ky_1 *  (dfs_dx_m1  - dfs_dx_p1 )  ;


        //// Hypervisocisty to stabilize simulation
        /*
        const double hypvisc_phi_val = -1.e-5;
        
        const CComplex d4_dx_phi    = (-39. *(Fields[Field::phi][s][m][z][y_k][x+1] - Fields[Field::phi][s][m][z][y_k][x-1])  
                                      + 12. *(Fields[Field::phi][s][m][z][y_k][x+2] - Fields[Field::phi][s][m][z][y_k][x-2]) 
                                      + 56. * Fields[Field::phi][s][m][z][y_k][x  ])/pow4(dx);

        const CComplex hypvisc_phi    = hypvisc_phi_val * ( d4_dx_phi + pow4(ky) * phi_);
*/

        /////////////// Finally the Vlasov equation calculate the time derivatve      //////////////////////
             
        const CComplex dg_dt = (y_k == 0) ? 0. : 
  //           +  hypvisc_phi
             -  IslandPhase * alpha * V[v] * (Island_g + sigma * Island_phi * f0_) +                // Island term
             +  nonLinearTerm[y_k][x][v]                                             // Non-linear ( array is zero for linear simulations) 
             +  ky* (-(w_n + w_T * (((V[v]*V[v])+ M[m])*kw_T  - sub)) * f0_ * 
//                 phi_    )                     // Source term (Temperature/Density gradient)
           phi_)   //  (phi_ - (y_k == 1 ? V[v] * MagIs[x] : 0.  )  ) )                       // Source term (Temperature/Density gradient)
     //        -  half_eta_kperp2_phi * f0_)                                            // Contributions from gyro-1 (0 if not neq Gyro-1)
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




void VlasovIsland::Vlasov_2D_Island_EM(
                           CComplex fs       [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex fss      [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f0 [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f1 [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex ft       [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Coll      [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           CComplex Fields[Nq][NsLD][NmLD][NzLB][Nky][NxLB+4]      ,
                           CComplex    nonLinearTerm               [Nky][NxLD  ][NvLD],
                           CComplex Xi       [NzLB][Nky][NxLB+4][NvLB],
                           CComplex G        [NzLB][Nky][NxLB][NvLB],
                           CComplex Xi_lin       [NzLB][Nky][NxLB+4][NvLB],
                           CComplex G_lin        [NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Ap_mod  [NzLB][Nky][NxLB+4]      ,
                           CComplex Field0[Nq][NzLD][Nky][NxLD]   ,
                           const double dt, const int rk_step, const double rk[3])
{ 

  // Add modified Ap from the Island to the field equations (gyro-averaging neglected)
   for(int s = NsLlD; s <= NsLuD; s++) { for(int m = NmLlD; m <= NmLuD; m++) { 
  
   Fields[Field::Ap][s][m][NzLlB:NzLB][0:Nky][NxLlB-2:NxLB+4] = Ap_mod[NzLlB:NzLB][0:Nky][NxLlB-2:NxLB+4];
   
   } }
   
   Field0[Field::Ap][NzLlD:NzLD][0:Nky][NxLlD:NxLD] = Ap_mod[NzLlD:NzLD][0:Nky][NxLlD:NxLD];


   ////
   for(int s = NsLlD; s <= NsLuD; s++) {
        
      // abbrevations
      const double w_n   = species[s].w_n;
      const double w_T   = species[s].w_T;
      const double alpha = species[s].alpha;
      const double sigma = species[s].sigma;
      const double Temp  = species[s].T0;
    
      const double sub   = (species[s].doGyro) ? 3./2. : 1./2.;
      
      bool isGyro1 = (species[s].gyroModel == "Gyro-1");
      

    for(int m = NmLlD; m <= NmLuD; m++) { 
 
          if(doNonLinear) setupXiAndG    (fs, f0 , Fields, Xi    , G    , m , s);
          else            setupXiAndG_lin(fs, f0 , Fields, Xi_lin, G_lin, m , s);
       
    for(int z = NzLlD; z <= NzLuD; z++) { 
      
          if(doNonLinear && (rk_step != 0)) calculatePoissonBracket(G    , Xi    , nullptr, nullptr, z, m, s, nonLinearTerm, Xi_max, true);
          else                              calculatePoissonBracket(G_lin, Xi_lin, nullptr, nullptr, z, m, s, nonLinearTerm, Xi_max, true);
      
    for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { for(int x = NxLlD; x <= NxLuD; x++) { 

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
      const CComplex f0_  = f0[s][m][z][y_k][x][v];

      const CComplex G_   =  G[z][y_k][x][v];
      const CComplex Xi_  = Xi[z][y_k][x][v];

    
    /////////////// Finally the Vlasov equation calculate the time derivatve      //////////////////////
   
            
    const CComplex dg_dt = 
    
    +  nonLinearTerm[y_k][x][v]                                             // Non-linear ( array is zero for linear simulations) 
    -  ky* (w_n + w_T * ((V[v]*V[v]+ M[m])/Temp  - sub)) * f0_ * phi_    // Driving term (Temperature/Density gradient)
    -  alpha  * V[v]* kp  * ( g + sigma * phi_ * f0_)                        // Linear Landau damping
    +  Coll[s][m][z][y_k][x][v]  ;                                          // Collisional operator
         
        
    //////////////////////////// Vlasov End ////////////////////////////
  
    //  time-integrate the distribution function    
    ft [s][m][z][y_k][x][v] = rk[0] * ft[s][m][z][y_k][x][v] +  rk[1] * dg_dt             ;
    fss[s][m][z][y_k][x][v] =         f1[s][m][z][y_k][x][v] + (rk[2] * ft[s][m][z][y_k][x][v] + dg_dt) * dt;



   } } } 
   } }
   }
}

void VlasovIsland::setupXiAndG_lin(
                           const CComplex g          [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f0         [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Fields [Nq][NsLD][NmLD][NzLB][Nky][NxLB+4]      ,
                           CComplex Xi                           [NzLB][Nky][NxLB+4][NvLB],
                           CComplex G                            [NzLB][Nky][NxLB  ][NvLB],
                           const int m, const int s) 
{

  // small abbrevations
  const double alpha = species[s].alpha;
  const double sigma = species[s].sigma;
  
  const double aeb   =  alpha * geo->eps_hat * plasma->beta; 
  const double saeb  =  sigma * aeb;

  const bool useAp   = (Nq >= 2);
  const bool useBp   = (Nq >= 3);


  // ICC vectorizes useAp/useBp into separate lopps, check for any speed penelity ? 
  #pragma omp for collapse(2)
  for(int z = NzLlB; z <= NzLuB; z++) {      for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { 
  for(int x = NxLlB; x <= NxLuB; x++) { simd_for(int v   = NvLlB ;   v <= NvLuB ;   v++) { 

    Xi[z][y_k][x][v] = - aeb*V[v]*Fields[Field::Ap][s][m][z][y_k][x]; 

    G [z][y_k][x][v] = g[s][m][z][y_k][x][v] + sigma * Fields[Field::phi][s][m][z][y_k][x] * f0[s][m][z][y_k][x][v];

    // substract canonical momentum to get "real" f1 (not used "yet")
    // f1[z][y_k][x][v] = g[s][m][z][y_k][x][v] - (useAp ? saeb * V[v] * f0[s][n][z][y_k][x][v] * Ap[s][m][z][y_k][x] : 0.);
      
  } } // v, x
     
  // Note we have extended boundaries in X (NxLlB-2 -- NxLuB+2) for fields
  // Intel Inspector complains about useAp ? ... memory violation, is it true? 
  for(int v = NvLlB; v <= NvLuB; v++) {

    Xi[z][y_k][NxLlB-2:2][v] =  - aeb*V[v]*Fields[Field::Ap][s][m][z][y_k][NxLlB-2:2];

    Xi[z][y_k][NxLuB+1:2][v] =  - aeb*V[v]*Fields[Field::Ap][s][m][z][y_k][NxLuB+1:2];

  }
  
  } } // y_k, z

}

