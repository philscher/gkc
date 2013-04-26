/*
 * =====================================================================================
 *
 *       Filename: Vlasov_Cilk.cpp
 *
 *    Description: Implementation of GK Vlasov equation using 
 *                 Intel Cilk (Array Notation)
 *
 *         Author: Paul P. Hilscher (2011-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#include "Vlasov_Cilk.h"


VlasovCilk::VlasovCilk(Grid *_grid, Parallel *_parallel, Setup *_setup, FileIO *fileIO, Geometry *_geo, FFTSolver *fft, Benchmark *_bench, Collisions *_coll)    
: Vlasov(_grid, _parallel, _setup, fileIO, _geo, fft, _bench, _coll)
{

}


void VlasovCilk::solve(std::string equation_type, Fields *fields, CComplex *f_in, CComplex *f_out, 
                       double dt, int rk_step, const double rk[3]) 
{

  if(0);
  else if(equation_type == "EM") Vlasov_EM((A6zz) f_in, (A6zz) f_out, (A6zz) f0, (A6zz) f, (A6zz) ft, (A6zz) Coll, 
                                           (A6zz) fields->Field, (A4zz) Xi, (A4zz) G, (A3zz) nonLinearTerm,
                                           (A2rr) geo->Kx, (A2rr) geo->Ky, (A2rr) geo->dB_dz,
                                           dt, rk_step, rk);
  else   check(-1, DMESG("No Such Equation"));

  return;

}


// Nice stuff http://stackoverflow.com/questions/841433/gcc-attribute-alignedx-explanation
//            http://stackoverflow.com/questions/3876758/working-around-the-char-array-on-stack-align-problem
// Align arrays, allocate on stack, (TAKE care of stackoverflow, if happens allocated
// dynamically with allocation in Constructor)
// Speed should not be a concern, as allocation happens instantaneously
// however, check this, what about alignment ?
// how to deal with indexes ?
// from stackoverflow guru : http://stackoverflow.com/questions/161053/c-which-is-faster-stack-allocation-or-heap-allocation
// Stack is hot, as it probably resides directly in cache. Good. 
// Hope stackoverflow is right and we don't get a stackoverflow...
//
// Stack variables are aligned per default on 8 bytes (at least for SSE2+ CPUs)
// Can be controlled on gcc with  -mpreferred-stack-boundary=n
// See : http://stackoverflow.com/questions/1061818/stack-allocation-padding-and-alignment
// short array[3] __attribute__ ((aligned (__BIGGEST_ALIGNMENT__)));
// __BIGGEST_ALIGNMENT automatically uses max alignments sizes supported by vector instructions
// take care, for electro-static simulations, G & Xi are null pointers (for em &phi respectively)
//
// Take care that OMP_STACKSIZE is large enough in case of thread parallelization
//
void VlasovCilk::calculateExBNonLinearity(const CComplex  G              [NzLB][Nky][NxLB  ][NvLB],  // in case of em
                                         const CComplex Xi              [NzLB][Nky][NxLB+4][NvLB],  // in case of em
                                         const CComplex  f [NsLD][NmLD ][NzLB][Nky][NxLB  ][NvLB],  // in case of es
                                         const CComplex Fields[Nq][NsLD][NmLD ][NzLB][Nky][NxLB+4], // in case of es
                                         const int z, const int m, const int s,
                                         CComplex ExB[Nky][NxLD][NvLD], double Xi_max[3], const bool electroMagnetic)
{
  // phase space function & Poisson bracket
  const double _kw_fft_Norm = 1./(fft->Norm_Y_Backward * fft->Norm_Y_Backward * fft->Norm_Y_Forward);
   
  const doubleAA _kw_12_dx  = 1./(12.*dx), _kw_12_dy=1./(12.*dy);
  const doubleAA _kw_24_dx  = 1./(24.*dx), _kw_24_dy=1./(24.*dy);
      
  CComplexAA  xky_Xi [Nky][NxLB+4];
  CComplexAA  xky_f1 [Nky][NxLB  ];
  CComplexAA  xky_ExB[Nky][NxLD  ];

  doubleAA    xy_Xi    [NyLD+8][NxLB+4]; // extended BC 
  doubleAA    xy_dXi_dy[NyLD+4][NxLB  ]; // normal BC
  doubleAA    xy_dXi_dx[NyLD+4][NxLB  ];
  doubleAA    xy_f1    [NyLD+4][NxLB  ];
  doubleAA    xy_ExB   [NyLD  ][NxLD  ];

  bool have_xy_dXi = false; // need for OpenMP parallelization over v
                            //  as all variables are on stack and thread local

  #pragma omp for
  for(int v = NvLlD; v <= NvLuD; v++) { 

    // Transform Xi to real space  
        
    // In case of electro-static simulations Xi[x][k_y] and no v-dependence
    // for electro-static field this has to be calculated only once
    if(electroMagnetic || !have_xy_dXi) {

      if(electroMagnetic) xky_Xi[:][:] =     Xi                  [z][:][NxLlB-2:NxLB+4][v];
      else                xky_Xi[:][:] = Fields[Field::phi][s][m][z][:][NxLlB-2:NxLB+4]   ;
       
      // xy_Xi[shift by +4][], as we also will use extended BC in Y
      // xy_Xi[Nky][:] = 0.; // we do not include Nyquist frequency
      fft->solve(FFT_Type::Y_FIELDS, FFT_Sign::Backward, xky_Xi, &xy_Xi[4][0]);
       
      // Set Periodic-Boundary in Y (in X is not necessary as we transform it too)
      //wonder if this is faster : xy_Xi[0     :4][4:NxLB] =  xy_Xi[NyLD  :4][4:NxLB];
      xy_Xi[0     :4][:] =  xy_Xi[NyLD  :4][:];
      xy_Xi[NyLD+4:4][:] =  xy_Xi[4     :4][:];

      // perform CD-4 derivative for dphi_dx , and dphi_dy (Note, we have extended GC in X&Y)
      for(int y=2; y < NyLB+2; y++) { simd_for(int x = 2; x < NxLB+2; x++)  {

         xy_dXi_dx[y-2][x-2] = (8.*(xy_Xi[y][x+1] - xy_Xi[y][x-1]) - (xy_Xi[y][x+2] - xy_Xi[y][x-2])) * _kw_12_dx;
         xy_dXi_dy[y-2][x-2] = (8.*(xy_Xi[y+1][x] - xy_Xi[y-1][x]) - (xy_Xi[y+2][x] - xy_Xi[y-2][x])) * _kw_12_dy;
      } } 
      
      // get maximum partial_nu chi value to calculate CFL condition
      {
        const double max_dXi_dx =  __sec_reduce_max(cabs(xy_dXi_dx[:][:]));
        const double max_dXi_dy =  __sec_reduce_max(cabs(xy_dXi_dy[:][:]));
       
        #pragma omp critical
        Xi_max[DIR_X] = std::max(Xi_max[DIR_X], max_dXi_dx);
        #pragma omp critical
        Xi_max[DIR_Y] = std::max(Xi_max[DIR_Y], max_dXi_dy);
      
        //Xi_max[DIR_X] = std::max(Xi_max[DIR_X], __sec_reduce_max(cabs(xy_dXi_dx[:][:])));
        //Xi_max[DIR_Y] = std::max(Xi_max[DIR_Y], __sec_reduce_max(cabs(xy_dXi_dy[:][:])));
      }
      have_xy_dXi = true; 
    }

    if(electroMagnetic) xky_f1[:][:] = G      [z][0:Nky][NxLlB:NxLB][v];
    else                xky_f1[:][:] = f[s][m][z][0:Nky][NxLlB:NxLB][v];

    fft->solve(FFT_Type::Y_PSF, FFT_Sign::Backward, (CComplex *) xky_f1, &xy_f1[2][0]);

    // Boundary in Y (In X is not necessary as we transformed it too), take care of FFT-boundary conditions
    xy_f1[0     :2][:] =  xy_f1[NyLD  :2][:];
    xy_f1[NyLD+2:2][:] =  xy_f1[2     :2][:];

    /////////////////   calculate cross terms using Morinishi scheme (Arakawa type) [Xi,G] (or [phi,F1])  /////////////////////
    for(int y = 2; y < NyLD+2; y++) { simd_for(int x = 2; x < NxLD+2; x++) {

      const double dXi_dy__dG_dx =  ( 8. * ( (xy_dXi_dy[y][x] + xy_dXi_dy[y][x+1]) * xy_f1[y][x+1]
                                           - (xy_dXi_dy[y][x] + xy_dXi_dy[y][x-1]) * xy_f1[y][x-1])
                                      -    ( (xy_dXi_dy[y][x] + xy_dXi_dy[y][x+2]) * xy_f1[y][x+2]
                                           - (xy_dXi_dy[y][x] + xy_dXi_dy[y][x-2]) * xy_f1[y][x-2]) ) * _kw_24_dx;
            
      const double dXi_dx__dG_dy =  ( 8. * ( (xy_dXi_dx[y][x] + xy_dXi_dx[y+1][x]) * xy_f1[y+1][x]
                                           - (xy_dXi_dx[y][x] + xy_dXi_dx[y-1][x]) * xy_f1[y-1][x])
                                      -    ( (xy_dXi_dx[y][x] + xy_dXi_dx[y+2][x]) * xy_f1[y+2][x] 
                                           - (xy_dXi_dx[y][x] + xy_dXi_dx[y-2][x]) * xy_f1[y-2][x]) ) * _kw_24_dy;
            
      // Take care of Fourier normalization : A*sqrt(N) * B*sqrt(N) 
      xy_ExB[y-2][x-2]  = (dXi_dy__dG_dx - dXi_dx__dG_dy) * _kw_fft_Norm;
    
    } } // x,y
   
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    fft->solve(FFT_Type::Y_NL, FFT_Sign::Forward, xy_ExB, (CComplex *) xky_ExB);

    // Done - store the non-linear term in ExB
    ExB[0:Nky][NxLlD:NxLD][v] = xky_ExB[:][:];
  }
    
  return;
}
                           
void VlasovCilk::setupXiAndG(
                           const CComplex g          [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex f0         [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                           const CComplex Fields [Nq][NsLD][NmLD][NzLB][Nky][NxLB+4]      ,
                           CComplex Xi                           [NzLB][Nky][NxLB+4][NvLB],
                           CComplex G                            [NzLB][Nky][NxLB  ][NvLB],
                           const int m, const int s) 
{

  // small abbreviations
  const double alpha = species[s].alpha;
  const double sigma = species[s].sigma;
  
  const double aeb   =  alpha * geo->eps_hat * plasma->beta; 
  const double saeb  =  sigma * aeb;

  const bool useAp   = (Nq >= 2);
  const bool useBp   = (Nq >= 3);

  // do we need boundaries in velocity space ?!

  // ICC vectorizes useAp/useBp into separate loops, check for any speed penalty ? 
  #pragma omp for collapse(2)
  for(int z = NzLlB; z <= NzLuB; z++) {      for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { 
  for(int x = NxLlB; x <= NxLuB; x++) { simd_for(int v   = NvLlB ;   v <= NvLuB ;   v++) { 

     Xi[z][y_k][x][v] = Fields[Field::phi][s][m][z][y_k][x] - (useAp ? aeb*V[v]*Fields[Field::Ap][s][m][z][y_k][x] : 0.) 
                                                            - (useBp ? aeb*M[m]*Fields[Field::Bp][s][m][z][y_k][x] : 0.);

     G [z][y_k][x][v] = g[s][m][z][y_k][x][v]  + sigma * Xi[z][y_k][x][v] * f0[s][m][z][y_k][x][v];

     // subtract canonical momentum to get "real" f1 (not used "yet")
     // f1[z][y_k][x][v] = g[s][m][z][y_k][x][v] - (useAp ? saeb * V[v] * f0[s][n][z][y_k][x][v] * Ap[s][m][z][y_k][x] : 0.);
      
  } } // v, x
     
  // Note we have extended boundaries in X (NxLlB-2 -- NxLuB+2) for fields
  // Intel Inspector complains about useAp ? ... memory violation, is it true? 
  simd_for(int v = NvLlB; v <= NvLuB; v++) {

    Xi[z][y_k][NxLlB-2:2][v] = Fields[Field::phi][s][m][z][y_k][NxLlB-2:2]  - (useAp ? aeb*V[v]*Fields[Field::Ap][s][m][z][y_k][NxLlB-2:2] : 0.)
                                                                            - (useBp ? aeb*M[m]*Fields[Field::Bp][s][m][z][y_k][NxLlB-2:2] : 0.);

    Xi[z][y_k][NxLuB+1:2][v] = Fields[Field::phi][s][m][z][y_k][NxLuB+1:2]  - (useAp ? aeb*V[v]*Fields[Field::Ap][s][m][z][y_k][NxLuB+1:2] : 0.) 
                                                                            - (useBp ? aeb*M[m]*Fields[Field::Bp][s][m][z][y_k][NxLuB+1:2] : 0.);
  }
  
  } } // y_k, z

}



void VlasovCilk::Vlasov_EM(
    const CComplex g   [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],  // Current step phase-space function
    CComplex       h   [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],  // Phase-space function for next step
    const CComplex f0  [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],  // Background Maxwellian
    const CComplex f1  [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],  // previous RK-Step
    CComplex       ft  [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],  // previous RK-Step
    CComplex       Coll[NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],  // Collisional corrections
    const CComplex Fields[Nq][NsLD][NmLD][NzLB][Nky][NxLB+4],
    CComplex Xi             [NzLB][Nky][NxLB+4][NvLB],
    CComplex G              [NzLB][Nky][NxLB  ][NvLB],
    CComplex NonLinearTerm        [Nky][NxLD  ][NvLD],        // Non-Linear Term (ExB)
    const double Kx[NzLD][NxLD], const double Ky[NzLD][NxLD], const double dB_dz[NzLD][NxLD], // Geometry stuff
    const double dt, const int rk_step, const double rk[3])
{ 

  const double B0 = plasma->B0;

  for(int s = NsLlD; s <= NsLuD; s++) {
        
    // small abbreviations
    const double w_n   = species[s].w_n;
    const double w_T   = species[s].w_T;
    const double alpha = species[s].alpha;
    const double sigma = species[s].sigma;
    const double Temp  = species[s].T0;
    

  for(int m = NmLlD; m <= NmLuD; m++) { 

    // Calculate before z loop as we use dg_dz and dXi_dz derivative
    setupXiAndG(g, f0, Fields, Xi, G, m, s);
      
    // Nested Parallelism (Poisson-bracket variables are allocated on stack)
    // Cannot collapse as we have no perfectly nested loops
          
  for(int z = NzLlD; z <= NzLuD; z++) { 
           
    // calculate non-linear term (rk_step == 0 for eigenvalue)
    if(doNonLinear         && (rk_step != 0)) calculateExBNonLinearity(G, Xi, nullptr, nullptr, z, m, s, NonLinearTerm, Xi_max, true); 
    if(doNonLinearParallel && (rk_step != 0)) calculateParallelNonLinearity(g, Fields, z, m, s, NonLinearTerm);
    
  // Note : we do not evolve highest mode (Nyquist)
  #pragma omp for collapse(2) 
  for(int y_k = 0; y_k < Nky-1; y_k++) { for(int x = NxLlD; x <= NxLuD; x++) { 
           
       
    const CComplex phi_    = Fields[Field::phi][s][m][z][y_k][x];

    const CComplex dphi_dx = (8.*(Fields[Field::phi][s][m][z][y_k][x+1] - Fields[Field::phi][s][m][z][y_k][x-1]) 
                               - (Fields[Field::phi][s][m][z][y_k][x+2] - Fields[Field::phi][s][m][z][y_k][x-2])) * _kw_12_dx  ;  
            
    const CComplex ky      = _imag * fft->ky(y_k);

    const double CoJB      = 1.;///geo->get_J(x,z);
    
  simd_for(int v = NvLlD; v <= NvLuD; v++) {

    const CComplex g_   =  g[s][m][z][y_k][x][v];
    const CComplex f0_  = f0[s][m][z][y_k][x][v];

    const CComplex G_   =  G[z][y_k][x][v];
    const CComplex Xi_  = Xi[z][y_k][x][v];

        
    const CComplex dg_dv = (8. *(g[s][m][z][y_k][x][v+1] - g[s][m][z][y_k][x][v-1]) 
                              - (g[s][m][z][y_k][x][v+2] - g[s][m][z][y_k][x][v-2])) * _kw_12_dv;


    const CComplex dG_dx = (8.*(G[z][y_k][x+1][v] - G[z][y_k][x-1][v]) 
                           -   (G[z][y_k][x+2][v] - G[z][y_k][x-2][v])) * _kw_12_dx;
    
    const CComplex dG_dz = (8.*(G[z+1][y_k][x][v] - G[z-1][y_k][x][v])    
                           -   (G[z+2][y_k][x][v] - G[z-2][y_k][x][v])) * _kw_12_dz;


    // magnetic pre-factor defined as  $ \hat{B}_0 / \hat{B}_{0\parallel}^\star = \left[ 1 + \beta_{ref} \sqrt{\frac{\hat{m_\sigma T_{0\sigma}{2}}}}
    // note j0 is calculated and needs to be replaced, or ? no we calculate j1 ne ?!
    const double j0   = 0.;
    const double Bpre = 1.; //1./(1. + plasma->beta * sqrt(m * T/2.) * j0 / (q * pow2(geo->B(x,y,z))) * V[v]);

        
    ///////////////   The time derivative of the Vlasov equation      //////////////////////
        
    const CComplex dg_dt = 
            
      NonLinearTerm[y_k][x][v]                                                              // Non-linear ( array is zero for linear simulations) 
    - Bpre * (w_n + w_T * ((pow2(V[v])+ M[m] * B0)/Temp - 3./2.)) * f0_ * Xi_ * ky          // Source Term
    - Bpre * sigma * ((M[m] * B0 + 2.*pow2(V[v]))/B0) *                                   
      (Kx[z][x] * dG_dx - Ky[z][x] * ky * G_)                                               // Magnetic curvature term
    + alpha * pow2(V[v]) * plasma->beta * plasma->w_p * G_ * ky                             // Plasma pressure gradient
    - CoJB * alpha * V[v]* dG_dz                                                            // Linear Landau damping term
    + alpha  / 2. * M[m] * dB_dz[z][x] * dg_dv                                              // Magnetic mirror term    
//    + Bpre *  sigma * (M[m] * B0 + 2. * pow2(V[v]))/B0 * Kx[z][x] * 
//    + ((w_n + w_T * (pow2(V[v]) + M[m] * B0)/Temp - 3./2.) * dG_dx + sigma * dphi_dx * f0_) // ??
    + Coll[s][m][z][y_k][x][v];                                                             // Collision term
          
    //////////////////////////// Vlasov End ////////////////////////////

    //  time-integrate the distribution function    
    ft[s][m][z][y_k][x][v] = rk[0] * ft[s][m][z][y_k][x][v] +  rk[1] * dg_dt             ;
    h [s][m][z][y_k][x][v] =         f1[s][m][z][y_k][x][v] + (rk[2] * ft[s][m][z][y_k][x][v] + dg_dt) * dt;
        
    } } // v  y_k

    } } }
  }
}



void VlasovCilk::printOn(std::ostream &output) const
{
  Vlasov::printOn(output);

  return;
}


void VlasovCilk::initData(Setup *setup, FileIO *fileIO) 
{
  Vlasov::initData(setup, fileIO);

  return;
}

   
void VlasovCilk::calculateParallelNonLinearity(
                                const CComplex f          [NsLD][NmLD][NzLB][Nky][NxLB   ][NvLB],
                                const CComplex Fields [Nq][NsLD][NmLD][NzLB][Nky][NxLB+4], 
                                const int z, const int m, const int s                     ,
                                CComplex NonLinearTerm[Nky][NxLD][NvLD])
{

  const double _kw_fft_Norm = 1./(fft->Norm_Y_Backward * fft->Norm_Y_Backward * fft->Norm_Y_Forward);
   
  CComplexAA  xky_dg_dv  [Nky][NxLB];
  CComplexAA  xky_dphi_dz[Nky][NxLB];
  CComplexAA  xky_v_NL   [Nky][NxLD];
  
  doubleAA    xy_dg_dv  [NyLD+4][NxLB  ];
  doubleAA    xy_dphi_dz[NyLD+4][NxLB  ];
  doubleAA    xy_v_NL   [NyLD  ][NxLD  ];
 
  // phi
  xky_dphi_dz[:][0:NxLB] = (8.*(Fields[Field::phi][s][m][z+1][:][NxLlB:NxLB] - Fields[Field::phi][s][m][z-1][:][NxLlB:NxLB]) 
                            -  (Fields[Field::phi][s][m][z+2][:][NxLlB:NxLB] - Fields[Field::phi][s][m][z-2][:][NxLlB:NxLB])) * _kw_12_dz  ; 

  fft->solve(FFT_Type::Y_PSF, FFT_Sign::Backward, xky_dphi_dz, &xy_dphi_dz[2][0]);

  #pragma omp for
  for(int v = NvLlD; v <= NvLuD; v++) { 
  
   xky_dg_dv[:][0:NxLB]  = (8. *(f[s][m][z][:][NxLlB:NxLB][v+1] - f[s][m][z][:][NxLlB:NxLB][v-1]) 
                              - (f[s][m][z][:][NxLlB:NxLB][v+2] - f[s][m][z][:][NxLlB:NxLB][v-2])) * _kw_12_dv;
    
   fft->solve(FFT_Type::Y_PSF, FFT_Sign::Backward, xky_dg_dv, &xy_dg_dv[2][0]);
  
   // Multiply in real space
   xy_v_NL[:][:] = xy_dphi_dz[2:NyLD][2:NxLD] * xy_dg_dv[2:NyLD][2:NxLD] * _kw_fft_Norm;

   fft->solve(FFT_Type::Y_NL, FFT_Sign::Forward, xy_v_NL, (CComplex *) xky_v_NL);

   if(doNonLinear) NonLinearTerm[0:Nky][NxLlD:NxLD][v] += xky_v_NL[:][:];
   else            NonLinearTerm[0:Nky][NxLlD:NxLD][v]  = xky_v_NL[:][:];
  }
}

