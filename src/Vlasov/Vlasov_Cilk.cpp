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
#include "Benchmark.h"


VlasovCilk::VlasovCilk(Grid *_grid, Parallel *_parallel, Setup *_setup, FileIO *fileIO, Geometry *_geo, FFTSolver *fft, Benchmark *_bench)    
: Vlasov(_grid, _parallel, _setup, fileIO, _geo, fft, _bench)
{

    nonLinearTerms.resize(RxLD , RkyLD , RvLD);  nonLinearTerms = 0.e0;
        
    collisionBeta = setup->get("Vlasov.CollisionBeta", 0.);
    
    Vlasov::initDataOutput(fileIO);    
}


int VlasovCilk::solve(std::string equation_type, Fields *fields, Array6C _fs, Array6C _fss, double dt, int rk_step, const double rk[3], int user_boundary_type) 
{
  if(0);
  else if(equation_type == "Vlasov_EM") Vlasov_EM((A6zz) _fs.dataZero(), (A6zz) _fss.dataZero(), (A6zz) f0.dataZero(), (A6zz) f.dataZero(), (A6zz) ft.dataZero(), 
                                                  (A5zz) fields->phi.dataZero(), (A5zz) fields->Ap.dataZero(), (A5zz) fields->Bp.dataZero(),
                                                  (A4zz) Xi.dataZero(), (A4zz) G.dataZero(),
                                                  (A3zz) nonLinearTerms.dataZero(),
                                                  X.dataZero(), V.dataZero(), M.dataZero(),  dt, rk_step, rk);
  else   check(-1, DMESG("No Such Equation"));

  return GKC_SUCCESS;
}


void VlasovCilk::calculatePoissonBracket(CComplex Xi        [NzLB][NkyLD][NxLB][NvLB ],
                                         CComplex  f [NsLD][NmLD ][NzLB][NkyLD][NxLB][NvLB],
                                         const int z, const int m, const int s,
                                         CComplex ExB[NxLD][NkyLD][NvLD], double Xi_max[3])
{
   // phase space function & Poisson bracket
   const double fft_Norm = fft->Norm_Y_Backward * fft->Norm_Y_Backward * fft->Norm_Y_Forward;

   // align arrays, allocate on stack, (TAKE care of stackoverflow, if happens allocated
   // dynamically with alloc in Constructor)
   // Speed should not be a concern, as allocation happends instantenously 
   // however, check this, what about alignement ?
   CComplex  xky_Xi [NkyLD][NxLD+8];
   CComplex  xky_f1 [NkyLD][NxLD+4];
   CComplex  xky_ExB[NkyLD][NxLD  ];

   double    xy_Xi    [NkyLD][NxLD+8];
   double    xy_dXi_dy[NyLD] [NxLD+8];
   double    xy_dXi_dx[NyLD] [NxLD+8];
   double    xy_f1    [NyLD] [NxLD+4];
   double    xy_ExB   [NyLD] [NxLD  ];

   const double _kw_12_dx = 1./(12.*dx), _kw_12_dy=1./(12.*dy);
   const double _kw_24_dx = 1./(24.*dx), _kw_24_dy=1./(24.*dy);

   // stride is not good
   for(int v=NvLlD; v<=NvLuD;v++) { 

        // Transform Xi to real space   
        xky_Xi[:][:] = Xi[z][NkyLlD:NkyLD][NxLlD-4:NxLD+8][v];
        
        fft->solve(FFT_BACKWARD, FFT_Y_FIELDS, (CComplex *) xky_Xi, xy_Xi);
       
        // Boundary in Y (In X is not necessary as we transformed it too)
        xy_Xi[NkyLuD+1:2][NxLB+1] =  xy_Xi[NyLlD  :2][NxLB+1];
        xy_Xi[NkyLlD-1:2][NxLB+1] =  xy_Xi[NyLuD-1:2][NxLB+1];

        // perform CD-4 derivative for dphi_dx , and dphi_dy
        omp_for(int y=NyLlD-2; y<= NyLuD+2;y++) { simd_for(int x=NxLlD-2; x<= NxLuD+2;x++)  {

         xy_dXi_dx[y][x] = (8.*(xy_Xi[y][x+1] - xy_Xi[y][x-1]) - (xy_Xi[y][x+2] - xy_Xi[y][x-2])) * _kw_12_dx;
         xy_dXi_dy[y][x] = (8.*(xy_Xi[y+1][x] - xy_Xi[y-1][x]) - (xy_Xi[y+2][x] - xy_Xi[y-2][x])) * _kw_12_dy;
            

        } } 
        

        // get maximum value to calculate CFL condition 
        // OPTIM : (is there sec_reduce_max_abs ? if not create)
        Xi_max[DIR_Y] = max(Xi_max[DIR_Y], max(__sec_reduce_max(xy_dXi_dx[:][:]), -__sec_reduce_min(xy_dXi_dx[:][:])));
        Xi_max[DIR_X] = max(Xi_max[DIR_X], max(__sec_reduce_max(xy_dXi_dy[:][:]), -__sec_reduce_min(xy_dXi_dy[:][:])));
        
        // for electro-static field this has to be calculated only once

        xky_f1[:][:] = f[s][m][z][NkyLlD:NkyLD][NxLlB:NxLB][v];
        fft->solve(FFT_BACKWARD, FFT_Y_PSF, (CComplex *) xky_f1, xy_f1);

        // Boundary in Y (In X is not necessary as we transformed it too), take care of FFT-boundary conditions
        xy_f1[NkyLuD+1:2][NxLB+1] =  xy_f1[NyLlD  :2][NxLB+1];
        xy_f1[NkyLlD-1:2][NxLB+1] =  xy_f1[NyLuD-1:2][NxLB+1];
         
       
     /////////////////   calculate cross terms using Morinishi scheme (Arakawa type)  [ phi, F1]  /////////////////////
     omp_for(int y=NyLlD; y<=NyLuD;y++) { simd_for(int x= NxLlD; x <= NxLuD; x++) {

            const double dXi_dy__dG_dx =  ( 8.* ( (  xy_dXi_dy[y][x] + xy_dXi_dy[y][x+1]) * xy_f1[y][x+1]) 
                                                 - ( xy_dXi_dy[y][x] + xy_dXi_dy[y][x-1]) * xy_f1[y][x-1])
                                            - (  ( ( xy_dXi_dy[y][x] + xy_dXi_dy[y][x+2]) * xy_f1[y][x+2]) 
                                            -      ( xy_dXi_dy[y][x] + xy_dXi_dy[y][x-2]) * xy_f1[y][x-2]) * _kw_24_dx     ;
            
            const double dXi_dx__dG_dy =  ( 8.* ( (  xy_dXi_dx[y][x] + xy_dXi_dx[y+1][x]) * xy_f1[y+1][x]) 
                                                 - ( xy_dXi_dx[y][x] + xy_dXi_dx[y-1][x]) * xy_f1[y-1][x])
                                            - (  ( ( xy_dXi_dx[y][x] + xy_dXi_dx[y+2][x]) * xy_f1[y+2][x]) 
                                            -      ( xy_dXi_dx[y][x] + xy_dXi_dx[y-2][x]) * xy_f1[y-2][x]) * _kw_24_dy     ;
            
            // Take care of normalization : A*sqrt(N) * B*sqrt(N) 
            xy_ExB[y][x]  = (dXi_dy__dG_dx - dXi_dx__dG_dy)/fft_Norm;
      } } // x,y
   
      fft->solve(FFT_FORWARD, FFT_Y_NL, xy_ExB, (CComplex *) xky_ExB);

      // Done - stores the non-linear terms ExB
      ExB[:][:][v] = xky_ExB[:][:];
   }
   
    
   return;

}
                           

void VlasovCilk::setupXiAndG(
                           const CComplex g   [NsLD][NmLD ][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex f0  [NsLD][NmLD ][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex phi [NsLD][NmLD ][NzLB][NkyLD][NxLB+4],
                           const CComplex Ap  [NsLD][NmLD ][NzLB][NkyLD][NxLB+4],
                           const CComplex Bp  [NsLD][NmLD ][NzLB][NkyLD][NxLB+4],
                           CComplex Xi        [NzLB][NkyLD][NxLB][NvLB ],
                           CComplex G         [NzLB][NkyLD][NxLB][NvLB ],
                           const double V[NvGB], const double M[NmGB],
                           const int m, const int s) 
{

  // small abbrevations
  const double alpha = plasma->species(s).alpha;
  const double sigma = plasma->species(s).sigma;
  
  const double aeb   =  alpha* geo->eps_hat * plasma->beta; 
  const double saeb  =  sigma * alpha * geo->eps_hat * plasma->beta;

  const bool useAp = (plasma->nfields >= 2);
  const bool useBp = (plasma->nfields >= 3);


  // ICC vectorizes useAp/useBp into seperate lopps, check for any speed penelity ? 
  for(int z = NzLlB; z <= NzLuB; z++) {  omp_for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { 
  for(int x = NxLlB; x <= NxLuB; x++) { simd_for(int v   = NvLlB ;   v <= NvLuB ;   v++) { 

     Xi[z][y_k][x][v] = phi[s][m][z][y_k][x] - (useAp ? aeb*V[v]*Ap[s][m][z][y_k][x] : 0.) - (useBp ? aeb*M[m]*Bp[s][m][z][y_k][x] : 0.);
     G [z][y_k][x][v] = g[s][m][z][y_k][x][v]  + sigma * Xi[z][y_k][x][v] * f0[s][m][z][y_k][x][v];
 // f1[z][y_k][x][v] = g[s][m][z][y_k][x][v] - (useAp ? saeb * V(v) * f0[s][n][z][y_k][x][v] * Ap[s][m][z][y_k][x] : 0.);
      
  }} // v, x
     
  // Note we have extended boundaries in X (NxLlB-2 -- NxLuB+2) for fields
  simd_for(int v   = NvLlB ;   v <= NvLuB ;   v++) { 
     Xi[z][y_k][NxLlB-2:2][v] = phi[s][m][z][y_k][NxLlB-2:2] - (useAp ? aeb*V[v]*Ap[s][m][z][y_k][NxLlB-2:2] : 0.) - (useBp ? aeb*M[m]*Bp[s][m][z][y_k][NxLlB-2:2] : 0.);
     Xi[z][y_k][NxLuB+1:2][v] = phi[s][m][z][y_k][NxLuB+1:2] - (useAp ? aeb*V[v]*Ap[s][m][z][y_k][NxLuB+1:2] : 0.) - (useBp ? aeb*M[m]*Bp[s][m][z][y_k][NxLuB+1:2] : 0.);
  }
  
  
  }} // y_k, z

};



void VlasovCilk::Vlasov_EM(
                           CComplex   g       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB], // Current step phase-space function
                           CComplex   h      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],  // Phase-space function for next step
                           const CComplex f0 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           CComplex ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex phi[NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           const CComplex Ap [NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           const CComplex Bp [NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           CComplex Xi                    [NzLB][NkyLD][NxLB][NvLB],
                           CComplex G                     [NzLB][NkyLD][NxLB][NvLB],
                           CComplex ExB                         [NkyLD][NxLB][NvLB],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           const double dt, const int rk_step, const double rk[3])
{ 

   
   const double B0 = plasma->B0;

   const bool nonLinear = false;

   for(int s = NsLlD; s <= NsLuD; s++) {
        
      // small abbrevations
      const double w_n   = plasma->species(s).w_n;
      const double w_T   = plasma->species(s).w_T;
      const double alpha = plasma->species(s).alpha;
      const double sigma = plasma->species(s).sigma;
      const double Temp  = plasma->species(s).T0;
    
      const double sub = (plasma->species(s).doGyro) ? 3./2. : 1./2.;
      

      for(int m=NmLlD; m<= NmLuD;m++) { 
 
          setupXiAndG(g, f0 , phi, Ap, Bp, Xi, G, V, M, m , s);
       
         // calculate for estimation of CFL condition
         for(int z=NzLlD; z<= NzLuD;z++) { 
           
           
         // rk_step == 0 for eigenvalue calculations
         if(nonLinear && (rk_step == 0)) calculatePoissonBracket(Xi, g, s, m, z, ExB, Xi_max); 
         // what about parallel non-linearity ?     
           
         omp_for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int x=NxLlD; x<= NxLuD;x++) { 
           
       
            const CComplex phi_ = phi[s][m][z][y_k][x];

            const CComplex dphi_dx = (8.*(phi[s][m][z][y_k][x+1] - phi[s][m][z][y_k][x-1]) - (phi[s][m][z][y_k][x+2] - phi[s][m][z][y_k][x-2]))/(12.*dx)  ;  

            const CComplex ky = (CComplex (0. + 1.j)) * fft->ky(y_k);

     
      for(int v=NvLlD; v<= NvLuD;v++) {
        

          
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
        const double Bpre  = 1.; //1./(1. + plasma->beta * sqrt(m * T/2.) * j0 / (q * pow2(geo->B(x,y,z))) * V(v));
        const double CoJB = 1./geo->J(x,z);

        
        ///////////////   The time derivative of the Vlasov equation      //////////////////////
        
        const CComplex dg_dt = 
            
          Bpre * (w_n + w_T * ((pow2(V[v])+ M[m] * B0)/Temp - sub)) * f0_ * Xi_ * ky           // Driving Term
          - Bpre * sigma * ((M[m] * B0 + 2.*pow2(V[v]))/B0) * 
            (geo->Kx(x,z) * dG_dx - geo->Ky(x,z) * ky * G_)                                   // Magnetic curvature term
          //- alpha * pow2(V[v]) * plasma->beta * plasma->w_p * G_ * ky
          -  CoJB *  alpha * V[v]* dG_dz                                                      // Landau damping term
          + alpha  / 2. * M[m] * geo->dB_dz(x,z) * dg_dv                                     // Magnetic mirror term    
          + Bpre *  sigma * (M[m] * B0 + 2. * pow2(V[v]))/B0 * geo->Kx(x,z) * 
          ((w_n + w_T * (pow2(V[v]) + M[m] * B0)/Temp - sub) * dG_dx + sigma * dphi_dx * f0_); // ??
           + collisionBeta  * (g_  + alpha * V[v] * dg_dv + v2_rms * ddg_dvv);               // Lennard-Bernstein Collision term

          
        //////////////////////////// Vlasov End ////////////////////////////

        //  time-integrate the distribution function    
        ft[s][m][z][y_k][x][v] = rk[0] * ft[s][m][z][y_k][x][v] + rk[1] * dg_dt             ;
        h [s][m][z][y_k][x][v] = f1[s][m][z][y_k][x][v]         + (rk[2] * ft[s][m][z][y_k][x][v] + dg_dt) * dt;
        
      }}} }}
   }
}


//const double Krook_nu = krook_nu * ( (X[x] > 0.8 * Lx/2.) || (X[x] < -0.8 * Lx/2.))  ?  0.1 * pow2(abs(X[x]) - 0.8 * Lx/2.): 0.;



void VlasovCilk::printOn(ostream &output) const
{
   Vlasov::printOn(output);

};


void VlasovCilk::initDataOutput(FileIO *fileIO) 
{
                Vlasov::initDataOutput(fileIO); 

};

