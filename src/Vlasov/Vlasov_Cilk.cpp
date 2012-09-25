/*
 * =====================================================================================
 *
 *       Filename: Vlasov_Cilk.cpp
 *
 *    Description: Implementation of GK Vlasov's equation using 
 *                 Intel Cilk (Array Notation)
 *
 *         Author: Paul P. Hilscher (2011-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#include "Vlasov_Cilk.h"


VlasovCilk::VlasovCilk(Grid *_grid, Parallel *_parallel, Setup *_setup, FileIO *fileIO, Geometry *_geo, FFTSolver *fft, Benchmark *_bench)    
: Vlasov(_grid, _parallel, _setup, fileIO, _geo, fft, _bench)
{

    //nonLinearTerms.resize(RxLD , RkyLD , RvLD);  nonLinearTerms = 0.e0;
        
    collisionBeta = setup->get("Vlasov.CollisionBeta", 0.);
    
    Vlasov::initDataOutput(fileIO);    
}


int VlasovCilk::solve(std::string equation_type, Fields *fields, Array6C _fs, Array6C _fss, double dt, int rk_step, const double rk[3]) 
{
  if(0);
  else if(equation_type == "Vlasov_EM") Vlasov_EM((A6zz) _fs.dataZero(), (A6zz) _fss.dataZero(), (A6zz) f0.dataZero(), (A6zz) f.dataZero(), (A6zz) ft.dataZero(), 
                                                  (A6zz) fields->Field.dataZero(), (A4zz) Xi, (A4zz) G, (A3zz) nonLinearTerms,
                                                  X, V, M,  dt, rk_step, rk);
  else   check(-1, DMESG("No Such Equation"));

  return GKC_SUCCESS;
}


/*
struct my_float {
          float number;
}  __attribute__((aligned(0x1000)));
}
 my_float a[4] = { ...} 
 aligned them all

Nice stuff http://stackoverflow.com/questions/841433/gcc-attribute-alignedx-explanation
           http://stackoverflow.com/questions/3876758/working-around-the-char-array-on-stack-align-problem
           godd bless sof

*/




// Align arrays, allocate on stack, (TAKE care of stackoverflow, if happens allocated
// dynamically with alloc in Constructor)
// Speed should not be a concern, as allocation happends instantenously 
// however, check this, what about alignement ?
// how to deal with indexes ?
// from stackoverflow guru : http://stackoverflow.com/questions/161053/c-which-is-faster-stack-allocation-or-heap-allocation
// Stack is hot, as it probably resides direclty in cache. Good. 
// Hope stackoverflow is right and we dont get a stackoverflow...
//
// Stack variables are alligned per default on 8 bytes (at least for SSE2+ CPUs)
// Can be controlled on gcc with  -mpreferred-stack-boundary=n
// See : http://stackoverflow.com/questions/1061818/stack-allocation-padding-and-alignment
// short array[3] __attribute__ ((aligned (__BIGGEST_ALIGNMENT__)));
// __BIGGEST_ALIGNMENT autmatically uses max alignments sizes supported by vector instructions

// take care, for electro-static simulations, G & Xi are null pointers (for e-m f&phi respectively)
void VlasovCilk::calculatePoissonBracket(const CComplex  G              [NzLB][NkyLD][NxLB  ][NvLB],  // in case of e-m
                                         const CComplex Xi              [NzLB][NkyLD][NxLB+4][NvLB],  // in case of e-m
                                         const CComplex  f [NsLD][NmLD ][NzLB][NkyLD][NxLB  ][NvLB],  // in case of e-s
                                         const CComplex Fields[Nq][NsLD][NmLD ][NzLB][NkyLD][NxLB+4], // in case of e-s
                                         const int z, const int m, const int s,
                                         CComplex ExB[NkyLD][NxLD][NvLD], double Xi_max[3], const bool electroMagnetic)
{
   // phase space function & Poisson bracket
   const double _kw_fft_Norm = 1./(fft->Norm_Y_Backward * fft->Norm_Y_Backward * fft->Norm_Y_Forward);
  
   // all direclty at cache-lines
   typedef __declspec(align(64)) double     doubleAA;
   typedef __declspec(align(64)) CComplex CComplexAA;

   CComplexAA  xky_Xi [NkyLD][NxLB+4];
   CComplexAA  xky_f1 [NkyLD][NxLB  ];
   CComplexAA  xky_ExB[NkyLD][NxLD  ];

   doubleAA    xy_Xi    [NyLD+8][NxLB+4]; // extended BC 
   doubleAA    xy_dXi_dy[NyLD+4][NxLB  ]; // normal BC
   doubleAA    xy_dXi_dx[NyLD+4][NxLB  ];
   doubleAA    xy_f1    [NyLD+4][NxLB  ];
   doubleAA    xy_ExB   [NyLD  ][NxLD  ];

   const doubleAA _kw_12_dx = 1./(12.*dx), _kw_12_dy=1./(12.*dy);
   const doubleAA _kw_24_dx = 1./(24.*dx), _kw_24_dy=1./(24.*dy);

   // stride is not good
   for(int v=NvLlD; v<=NvLuD;v++) { 

        // Transform Xi to real space  
        
        // In case of electro-static simulations Xi[x][k_y] and no v-dependence
        // for electro-static field this has to be calculated only once
        if(electroMagnetic || (v == NvLlD)) {

        if(electroMagnetic) xky_Xi[:][:] =     Xi                  [z][NkyLlD:NkyLD][NxLlB-2:NxLB+4][v];
        else                xky_Xi[:][:] = Fields[Field::phi][s][m][z][NkyLlD:NkyLD][NxLlB-2:NxLB+4]   ;
       
        // xy_Xi[shift by +4][], as we also will use extended BC in Y
        fft->solve(FFT_Y_FIELDS, FFT_BACKWARD, xky_Xi, &xy_Xi[4][0]);
       
        // Set Periodic-Boundary in Y (in X is not necessary as we transform it too)
        //wonder if this is faster : xy_Xi[0     :4][4:NxLB] =  xy_Xi[NyLD  :4][4:NxLB];
        xy_Xi[0     :4][:] =  xy_Xi[NyLD  :4][:];
        xy_Xi[NyLD+4:4][:] =  xy_Xi[4     :4][:];

        // perform CD-4 derivative for dphi_dx , and dphi_dy (Note, we have extendend GC in X&Y)
        omp_for(int y=2; y < NyLB+2; y++) { simd_for(int x=2; x < NxLB+2; x++)  {

         xy_dXi_dx[y-2][x-2] = (8.*(xy_Xi[y][x+1] - xy_Xi[y][x-1]) - (xy_Xi[y][x+2] - xy_Xi[y][x-2])) * _kw_12_dx;
         xy_dXi_dy[y-2][x-2] = (8.*(xy_Xi[y+1][x] - xy_Xi[y-1][x]) - (xy_Xi[y+2][x] - xy_Xi[y-2][x])) * _kw_12_dy;
            

        } } 
        
        // get maximum value to calculate CFL condition 
        // OPTIM : (is there sec_reduce_max_abs ? if not create)
        Xi_max[DIR_Y] = max(Xi_max[DIR_Y], max(__sec_reduce_max(xy_dXi_dy[:][:]), -__sec_reduce_min(xy_dXi_dy[:][:])));
        Xi_max[DIR_X] = max(Xi_max[DIR_X], max(__sec_reduce_max(xy_dXi_dx[:][:]), -__sec_reduce_min(xy_dXi_dx[:][:])));
        

        }

        if(electroMagnetic) xky_f1[:][:] = G      [z][NkyLlD:NkyLD][NxLlB:NxLB][v];
        else                xky_f1[:][:] = f[s][m][z][NkyLlD:NkyLD][NxLlB:NxLB][v];

        fft->solve(FFT_Y_PSF, FFT_BACKWARD, (CComplex *) xky_f1, &xy_f1[2][0]);

        // Boundary in Y (In X is not necessary as we transformed it too), take care of FFT-boundary conditions
        xy_f1[0     :2][:] =  xy_f1[NyLD  :2][:];
        xy_f1[NyLD+2:2][:] =  xy_f1[2     :2][:];

     /////////////////   calculate cross terms using Morinishi scheme (Arakawa type) [Xi,G] (or [phi,F1])  /////////////////////
      
     omp_for(int y=2; y < NyLD+2; y++) { simd_for(int x=2; x < NxLD+2; x++) {

            const double dXi_dy__dG_dx =  ( 8. * ( (xy_dXi_dy[y][x] + xy_dXi_dy[y][x+1]) * xy_f1[y][x+1]
                                                 - (xy_dXi_dy[y][x] + xy_dXi_dy[y][x-1]) * xy_f1[y][x-1])
                                            -    ( (xy_dXi_dy[y][x] + xy_dXi_dy[y][x+2]) * xy_f1[y][x+2]
                                                 - (xy_dXi_dy[y][x] + xy_dXi_dy[y][x-2]) * xy_f1[y][x-2]) ) * _kw_24_dx;
            
            const double dXi_dx__dG_dy =  ( 8. * ( (xy_dXi_dx[y][x] + xy_dXi_dx[y+1][x]) * xy_f1[y+1][x]
                                                 - (xy_dXi_dx[y][x] + xy_dXi_dx[y-1][x]) * xy_f1[y-1][x])
                                            -    ( (xy_dXi_dx[y][x] + xy_dXi_dx[y+2][x]) * xy_f1[y+2][x] 
                                                 - (xy_dXi_dx[y][x] + xy_dXi_dx[y-2][x]) * xy_f1[y-2][x]) ) * _kw_24_dy;
            
           // Take care of Fourier normalization : A*sqrt(N) * B*sqrt(N) 
           xy_ExB[y-2][x-2]  = -(dXi_dy__dG_dx - dXi_dx__dG_dy) * _kw_fft_Norm;
      } } // x,y
   
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      fft->solve(FFT_Y_NL, FFT_FORWARD, xy_ExB, (CComplex *) xky_ExB);

      // Done - store the non-linear term in ExB
      ExB[NkyLlD:NkyLD][NxLlD:NxLD][v] = xky_ExB[:][:];

   }
    
   return;

}
                           

void VlasovCilk::setupXiAndG(
                           const CComplex g   [NsLD][NmLD ][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex f0  [NsLD][NmLD ][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex Fields [Nq][NsLD][NmLD ][NzLB][NkyLD][NxLB+4],
                           CComplex Xi        [NzLB][NkyLD][NxLB+4][NvLB ],
                           CComplex G         [NzLB][NkyLD][NxLB][NvLB ],
                           const double V[NvGB], const double M[NmGB],
                           const int m, const int s) 
{

  // small abbrevations
  const double alpha = plasma->species[s].alpha;
  const double sigma = plasma->species[s].sigma;
  
  const double aeb   =  alpha* geo->eps_hat * plasma->beta; 
  const double saeb  =  sigma * alpha * geo->eps_hat * plasma->beta;

  const bool useAp = (plasma->nfields >= 2);
  const bool useBp = (plasma->nfields >= 3);


  // ICC vectorizes useAp/useBp into seperate lopps, check for any speed penelity ? 
  for(int z = NzLlB; z <= NzLuB; z++) {  omp_for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { 
  for(int x = NxLlB; x <= NxLuB; x++) { simd_for(int v   = NvLlB ;   v <= NvLuB ;   v++) { 

     Xi[z][y_k][x][v] = Fields[Field::phi][s][m][z][y_k][x] - (useAp ? aeb*V[v]*Fields[Field::Ap][s][m][z][y_k][x] : 0.) 
                                                            - (useBp ? aeb*M[m]*Fields[Field::Bp][s][m][z][y_k][x] : 0.);

     G [z][y_k][x][v] = g[s][m][z][y_k][x][v]  + sigma * Xi[z][y_k][x][v] * f0[s][m][z][y_k][x][v];

     // substract canonical momentum to get "real" f1 (not used "yet")
     // f1[z][y_k][x][v] = g[s][m][z][y_k][x][v] - (useAp ? saeb * V[v] * f0[s][n][z][y_k][x][v] * Ap[s][m][z][y_k][x] : 0.);
      
  } } // v, x
     
  // Note we have extended boundaries in X (NxLlB-2 -- NxLuB+2) for fields
  simd_for(int v   = NvLlB ;   v <= NvLuB ;   v++) {

     Xi[z][y_k][NxLlB-2:2][v] = Fields[Field::phi][s][m][z][y_k][NxLlB-2:2] - (useAp ? aeb*V[v]*Fields[Field::Ap][s][m][z][y_k][NxLlB-2:2] : 0.)
                                                                           - (useBp ? aeb*M[m]*Fields[Field::Bp][s][m][z][y_k][NxLlB-2:2] : 0.);

     Xi[z][y_k][NxLuB+1:2][v] = Fields[Field::phi][s][m][z][y_k][NxLuB+1:2] - (useAp ? aeb*V[v]*Fields[Field::Ap][s][m][z][y_k][NxLuB+1:2] : 0.) 
                                                                           - (useBp ? aeb*M[m]*Fields[Field::Bp][s][m][z][y_k][NxLuB+1:2] : 0.);
  }
  
  
  }} // y_k, z

};



void VlasovCilk::Vlasov_EM(
                                 CComplex g  [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB], // Current step phase-space function
                                 CComplex h  [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],  // Phase-space function for next step
                           const CComplex f0 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                                 CComplex ft [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex Fields[Nq][NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                                 CComplex Xi             [NzLB][NkyLD][NxLB+4][NvLB],
                                 CComplex G              [NzLB][NkyLD][NxLB  ][NvLB],
                                 CComplex ExB                  [NkyLD][NxLD  ][NvLB],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           const double dt, const int rk_step, const double rk[3])
{ 

   
   const double B0 = plasma->B0;

   const bool nonLinear = false;

   for(int s = NsLlD; s <= NsLuD; s++) {
        
      // small abbrevations
      const double w_n   = plasma->species[s].w_n;
      const double w_T   = plasma->species[s].w_T;
      const double alpha = plasma->species[s].alpha;
      const double sigma = plasma->species[s].sigma;
      const double Temp  = plasma->species[s].T0;
    
      const double sub = (plasma->species[s].doGyro) ? 3./2. : 1./2.;
      

      for(int m=NmLlD; m<= NmLuD;m++) { 

          // Calculate before z loop as we use dg_dz and dXi_dz derivative
          setupXiAndG(g, f0 , Fields, Xi, G, V, M, m , s);
      
          // Nested Parallelism (PoissonBracket variables are allocated on stack)
          // Cannot collapse as we have no perfetly nested loops
         omp_for(int z=NzLlD; z<= NzLuD;z++) { 
           
           
         // calculate non-linear term (rk_step == 0 for eigenvalue calculatinullptr)
         // CFL condition is calculated inside calculatePoissonBracket
         if(nonLinear && (rk_step != 0)) calculatePoissonBracket(G, Xi, nullptr, nullptr, z, m, s, ExB, Xi_max, true); 
           
         omp_for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int x=NxLlD; x<= NxLuD;x++) { 
           
       
            const CComplex phi_ = Fields[Field::phi][s][m][z][y_k][x];

            const CComplex dphi_dx = (8.*(Fields[Field::phi][s][m][z][y_k][x+1] - Fields[Field::phi][s][m][z][y_k][x-1]) - (Fields[Field::phi][s][m][z][y_k][x+2] - Fields[Field::phi][s][m][z][y_k][x-2]))/(12.*dx)  ;  

            const CComplex ky = (CComplex (0. + 1.j)) * fft->ky(y_k);

     
      simd_for(int v=NvLlD; v<= NvLuD;v++) {
        

          
        const CComplex g_   = g[s][m][z][y_k][x][v];
        const CComplex f0_  = f0 [s][m][z][y_k][x][v];
        const CComplex G_   = G[z][y_k][x][v];
        const CComplex Xi_  = Xi[z][y_k][x][v];

        
        // Velocity derivaties for Lennard-Bernstein Collisional Model
        const CComplex dg_dv   = (8. *(g[s][m][z][y_k][x][v+1] - g[s][m][z][y_k][x][v-1]) - (g[s][m][z][y_k][x][v+2] - g[s][m][z][y_k][x][v-2]))/(12.*dv);
        const CComplex ddg_dvv = (16.*(g[s][m][z][y_k][x][v+1] + g[s][m][z][y_k][x][v-1]) - (g[s][m][z][y_k][x][v+2] + g[s][m][z][y_k][x][v-2]) - 30.*g[s][m][z][y_k][x][v])/(12.*dv*dv);
        const double v2_rms = 1.;//pow2(alpha)
    
        
        /////////////// Finally the Vlasov equation calculate the time derivatve      //////////////////////

        // We use CD-4 (central difference fourth order for every variable)

        const CComplex dG_dz   = (8.*(G[z+1][y_k][x][v] - G[z-1][y_k][x][v])    -1.*(G[z+2][y_k][x][v] - G[z-2][y_k][x][v]))/(12.*dz);
        const CComplex dG_dx   = (8.*(G[z][y_k][x+1][v] - G[z][y_k][x-1][v])    -1.*(G[z][y_k][x+2][v] - G[z][y_k][x-2][v]))/(12.*dx);

        
        // magnetic prefactor defined as  $ \hat{B}_0 / \hat{B}_{0\parallel}^\star = \left[ 1 + \beta_{ref} \sqrt{\frac{\hat{m_\sigma T_{0\sigma}{2}}}}
        // note j0 is calculated and needs to be replaced, or ? no we calculate j1 ne ?!
        const double j0 = 0.;
        const double Bpre  = 1.; //1./(1. + plasma->beta * sqrt(m * T/2.) * j0 / (q * pow2(geo->B(x,y,z))) * V[v]);
        const double CoJB = 1./geo->J(x,z);

        
        ///////////////   The time derivative of the Vlasov equation      //////////////////////
        
        const CComplex dg_dt = 
            
                ExB[y_k][x][v]                                                                 // Non-linear ( array is zero for linear simulations) 
          + Bpre * (w_n + w_T * ((pow2(V[v])+ M[m] * B0)/Temp - sub)) * f0_ * Xi_ * ky         // Driving Term
          - Bpre * sigma * ((M[m] * B0 + 2.*pow2(V[v]))/B0) *                                   
            (geo->Kx(x,z) * dG_dx - geo->Ky(x,z) * ky * G_)                                    // Magnetic curvature term
          //- alpha * pow2(V[v]) * plasma->beta * plasma->w_p * G_ * ky                        // Plasma pressure gradient
          -  CoJB *  alpha * V[v]* dG_dz                                                       // Landau damping term
          + alpha  / 2. * M[m] * geo->dB_dz(x,z) * dg_dv                                       // Magnetic mirror term    
          + Bpre *  sigma * (M[m] * B0 + 2. * pow2(V[v]))/B0 * geo->Kx(x,z) * 
          ((w_n + w_T * (pow2(V[v]) + M[m] * B0)/Temp - sub) * dG_dx + sigma * dphi_dx * f0_); // ??
           + collisionBeta  * (g_  + alpha * V[v] * dg_dv + v2_rms * ddg_dvv);                 // Lennard-Bernstein Collision term

          
        //////////////////////////// Vlasov End ////////////////////////////

        //  time-integrate the distribution function    
        ft[s][m][z][y_k][x][v] = rk[0] * ft[s][m][z][y_k][x][v] + rk[1] * dg_dt             ;
        h [s][m][z][y_k][x][v] = f1[s][m][z][y_k][x][v]         + (rk[2] * ft[s][m][z][y_k][x][v] + dg_dt) * dt;
        
      }}} }}
   }
}


//const double Krook_nu = krook_nu * ( (X[x] > 0.8 * Lx/2.) || (X[x] < -0.8 * Lx/2.))  ?  0.1 * pow2(abs(X[x]) - 0.8 * Lx/2.): 0.;



void VlasovCilk::printOn(std::ostream &output) const
{
   Vlasov::printOn(output);

};


void VlasovCilk::initDataOutput(FileIO *fileIO) 
{
                Vlasov::initDataOutput(fileIO); 

};

