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


#include "config.h"
#include "Vlasov_Cilk.h"

#include "Geometry/Geometry2D.h"

#ifndef __cilk
#include <cilk/cilk_stub.h>
#endif
#include <cilk/cilk.h>

#include "Benchmark.h"

typedef blitz::cmplxd(*A6z)[][][][][];
typedef blitz::cmplxd(*A5z)[][][][];
typedef blitz::cmplxd(*A4z)[][][];
typedef blitz::cmplxd(*A3z)[][];

void VlasovCilk::setBoundaryXY(Array3d A, int dir) {

  if((dir == DIR_X) || (dir == DIR_XY)) {
    if(parallel->decomposition(DIR_X) > 1) {
          SendXlR(RB, RyLD, RzLD) = A(Range(NxLlD, NxLlD+1),   RyLD, RzLD);
          SendXuR(RB, RyLD, RzLD) = A(Range(NxLuD-1, NxLuD),   RyLD, RzLD);
   	  
          parallel->updateNeighbours(SendXuR, SendXlR, RecvXuR, RecvXlR, DIR_X, false);
          
          A(Range(NxLlB, NxLlB+1), RyLD, RzLD) = RecvXlR(RB, RyLD, RzLD);
          A(Range(NxLuD+1, NxLuB), RyLD, RzLD) = RecvXuR(RB, RyLD, RzLD);
    } else {
         A(Range(NxLlB, NxLlB+1), RyLD, RzLD) = A(Range(NxLuD-1, NxLuD),   RyLD, RzLD);
         A(Range(NxLuD+1, NxLuB), RyLD, RzLD) = A(Range(NxLlD, NxLlD+1),   RyLD, RzLD);
    }
  }
  if((dir == DIR_Y) || (dir == DIR_XY)) {
        A(RxLD, Range(NyLlB, NyLlB+1), RzLD) = A(RxLD, Range(NyLuD-1, NyLuD), RzLD);
        A(RxLD, Range(NyLuD+1, NyLuB), RzLD) = A(RxLD, Range(NyLlD, NyLlD+1), RzLD);
  }

};

VlasovCilk::VlasovCilk(Grid *_grid, Parallel *_parallel, Setup *_setup, FileIO *fileIO, Geometry *_geo, FFTSolver *fft)    
    : Vlasov(_grid, _parallel, _setup, fileIO, _geo, fft),
      dphi_dx(FortranArray<3>()),   dphi_dy(FortranArray<3>()),   
      dAp_dx(FortranArray<3>()),    dAp_dy(FortranArray<3>()),    k2p_phi(FortranArray<4>()), nonLinearTerms(GKCStorage4)
{

    f1.resize(RxLB , RkyLD , RzLB, RvLB, RmLD, RsLD);  f1 = 0.e0;
    allocate(RxLB, RyLB, RzLD, xy_phi, xy_dphi_dx, xy_dphi_dy, xy_f1);
        
    nonLinearTerms.resize(RxLD , RkyLD , RzLD, RvLD);  nonLinearTerms = 0.e0;
    allocate(RxLB, RkyLD, RzLB, dphi_dy, dphi_dx, dAp_dy, dAp_dx);
    
    k2p_phi.resize(RxLD, RkyLD, RzLD, RFields) ; k2p_phi = 0.;
        
    collisionBeta = setup->get("Vlasov.CollisionBeta", 0.);

    
    // allocate boundary (mpi) buffers (for real values)
    allocate(RB, RyLD, RzLD, SendXlR, SendXuR, RecvXlR, RecvXuR);


    Vlasov::initDataOutput(fileIO);    
}


int VlasovCilk::solve(std::string equation_type, Fields *fields, Array6z _fs, Array6z _fss, double dt, int rk_step, const double rk[3], int user_boundary_type) 
{
   /// Enter Cilk Part here
  Xi_max = 0.;
 
  if((equation_type == "2D_ES")) {
  Benchmark bench(setup, parallel);
  bench.start("Vlasov");
  Vlasov_2D((A6z) _fs.dataZero(), (A6z) _fss.dataZero(), (A6z) f0.dataZero(), (A6z) f.dataZero(), (A6z) ft.dataZero(), (A5z) fields->phi.dataZero(), (A4z) k2p_phi.dataZero(), (A4z) nonLinearTerms.dataZero(), X.dataZero(), V.dataZero(), M.dataZero(), fields, dt, rk_step, rk, _fs);
  bench.stop("Vlasov");
  }
  else if((equation_type == "2D_EM")) Vlasov_EM_2D((A6z) _fs.dataZero(), (A6z) _fss.dataZero(), (A6z) f0.dataZero(), (A6z) f.dataZero(), (A6z) ft.dataZero(), (A5z) fields->phi.dataZero(), (A5z) fields->Ap.dataZero(), (A5z) fields->Bp.dataZero(), (A4z) k2p_phi.dataZero(), (A3z) dphi_dx.dataZero(), (A4z) Xi.dataZero(), (A4z) G.dataZero(), X.dataZero(), V.dataZero(), M.dataZero(), fields, dt, rk_step, rk);
  else if((equation_type == "Vlasov_EM")) Vlasov_EM((A6z) _fs.dataZero(), (A6z) _fss.dataZero(), (A6z) f0.dataZero(), (A6z) f.dataZero(), (A6z) ft.dataZero(), (A5z) fields->phi.dataZero(), (A5z) fields->Ap.dataZero(), (A5z) fields->Bp.dataZero(), (A4z) k2p_phi.dataZero(), (A3z) dphi_dx.dataZero(), (A4z) Xi.dataZero(), (A4z) G.dataZero(), X.dataZero(), V.dataZero(), M.dataZero(), fields, dt, rk_step, rk);
  else if((equation_type == "2DIsland")) Vlasov_2D_Island((A6z) _fs.dataZero(), (A6z) _fss.dataZero(), (A6z) f0.dataZero(), (A6z) f.dataZero(), (A6z) ft.dataZero(), (A5z) fields->phi.dataZero(), (A4z) k2p_phi.dataZero(), (A4z) nonLinearTerms.dataZero(), (A3z) dphi_dx.dataZero(), 
      X.dataZero(), V.dataZero(), M.dataZero(), fields, dt, rk_step, rk,  _fs);
  else if((equation_type == "2DLandauDamping")) Landau_Damping((A6z) _fs.dataZero(), (A6z) _fss.dataZero(), (A6z) f0.dataZero(), (A6z) f.dataZero(), (A6z) ft.dataZero(), (A5z) fields->phi.dataZero(), (A4z) k2p_phi.dataZero(), (A3z) dphi_dx.dataZero(), 
      X.dataZero(), V.dataZero(), M.dataZero(), fields, dt, rk_step, rk);
  else   check(-1, DMESG("No Such Equation"));

  return GKC_SUCCESS;
}


Array3d VlasovCilk::transform2Real(const Array5z A, const int m, const int s) {

   fft->kYIn(RxLD, RkyLD, RzLD, Field::phi)  = A(RxLD,RkyLD,RzLD, m, s);
   fft->solve(FFT_Y, FFT_BACKWARD, NxLD * NzLD);
   
   return  fft->rYOut(RxLD, RyLD, RzLD, Field::phi);

}

Array3d VlasovCilk::transform2Real(const Array6z A, const int v, const int m, const int s) {

   fft->kYIn(RxLD, RkyLD, RzLD, Field::phi)  = A(RxLD,RkyLD,RzLD, v, m, s);
   fft->solve(FFT_Y, FFT_BACKWARD, NxLD * NzLD);
   
   return  fft->rYOut(RxLD, RyLD, RzLD, Field::phi);

}



// We copy arrays many time, very inefficient. Improve
Array4z VlasovCilk::calculatePoissonBracket(Array5z phi, Array6z f1, const int m, const int s) 
{

   // Transform fields to real space
   transform2Real(phi, m , s); xy_phi(RxLD, RyLD, RzLD) = fft->rYOut(RxLD, RyLD, RzLD, Field::phi);
   setBoundaryXY(xy_phi,DIR_XY);

   // perform CD-4 derivative for dphi_dx , and dphi_dy
   for(int z = NzLlD; z <= NzLuD; z++) { omp_for(int y=NyLlD; y<= NyLuD;y++) { for(int x=NxLlD; x<= NxLuD;x++)  {

     xy_dphi_dx(x, y, z) = (8.*(xy_phi(x+1, y,z) - xy_phi(x-1, y, z)) - (xy_phi(x+2,y,z) - xy_phi(x-2,y,z)))/(12.*dx); 
     xy_dphi_dy(x, y, z) = (8.*(xy_phi(x, y+1,z) - xy_phi(x, y-1, z)) - (xy_phi(x,y+2,z) - xy_phi(x,y-2,z)))/(12.*dy); 

   } } }

   setBoundaryXY(xy_dphi_dx,DIR_Y);
   setBoundaryXY(xy_dphi_dy,DIR_X);
   
   
   // phase space function & Poisson bracket
   const double fft_Norm = fft->Norm_Y_Backward * fft->Norm_Y_Backward * fft->Norm_Y_Forward;
   for(int v=NvLlD; v<=NvLuD;v++) { 
   
   // Transform phase-space to real space   
   transform2Real(f1, v, m ,s); xy_f1(RxLD, RyLD, RzLD) = fft->rYOut(RxLD, RyLD, RzLD, Field::phi);
   setBoundaryXY(xy_f1,DIR_XY);
     
      
      ///////////////////////   calculate cross terms using Morinishi scheme    [ phi, F1]  //////////////////////////////////
      for(int z=NzLlD; z<=NzLuD;z++) { omp_for(int y=NyLlD; y<=NyLuD;y++) { for(int x= NxLlD; x <= NxLuD; x++) {

                const double dXi_dy__dG_dx =  ( 8.*((xy_dphi_dy(x,y,z)+xy_dphi_dy(x+1,y,z))*(xy_f1(x+1, y, z)) - (xy_dphi_dy(x,y,z)+xy_dphi_dy(x-1,y,z))*(xy_f1(x-1, y, z)))      
                                            - (     (xy_dphi_dy(x,y,z)+xy_dphi_dy(x+2,y,z))*(xy_f1(x+2, y, z)) - (xy_dphi_dy(x,y,z)+xy_dphi_dy(x-2,y,z))*(xy_f1(x-2, y, z))))/(24.*dx)     ;
                const double dXi_dx__dG_dy =  ( 8.*((xy_dphi_dx(x,y,z)+xy_dphi_dx(x,y+1,z))*(xy_f1(x, y+1, z)) - (xy_dphi_dx(x,y,z)+xy_dphi_dx(x,y-1,z))*(xy_f1(x, y-1, z)))      
                                            - (     (xy_dphi_dx(x,y,z)+xy_dphi_dx(x,y+2,z))*(xy_f1(x, y+2, z)) - (xy_dphi_dx(x,y,z)+xy_dphi_dx(x,y-2,z))*(xy_f1(x, y-2, z))))/(24.*dy)     ;
        
            // Take care of normalization : A*sqrt(N) * B*sqrt(N) 
	        fft->rYIn(x,y,z, Field::phi) =  (dXi_dy__dG_dx - dXi_dx__dG_dy)/fft_Norm;
	    
      } } }
    
      // back transform to Fourier space 
      fft->solve(FFT_Y, FFT_FORWARD, NxLD * NzLD);
      nonLinearTerms(RxLD,RkyLD,RzLD,v) = fft->kYOut(RxLD, RkyLD, RzLD, Field::phi);
    
      
    
   }
          
   return nonLinearTerms;

}
                           


void VlasovCilk::Vlasov_2D(
                           const cmplxd fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd vf0[NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd phi[NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           cmplxd    k2_phi[plasma->nfields][NzLD][NkyLD][NxLD],
                           cmplxd    nonLinear[NzLD][NkyLD][NxLD][NvLD],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           Fields *fields,
                           const double dt, const int rk_step, const double rk[3], Array6z _fs)
{ 

//  const Geometry2D *geo = static_cast<Geometry2D*>(Vlasov::geo);
       
  // increases speed by 5% 
  const double _12_dv    = 1./(12.*dv);
  const double _12_dv_dv = 1./(12.*dv*dv);
  const double _16_dx4 = 16.*pow4(dx);

/*
  __asm{
         nop
  }
  */ 
   #pragma ivdep
   #pragma vector always 
   for(int s = NsLlD; s <= NsLuD; s++) {
        
      // small abbrevations
      const double w_n   = plasma->species(s).w_n;
      const double w_T   = plasma->species(s).w_T;
      const double alpha = plasma->species(s).alpha;
      const double sigma = plasma->species(s).sigma;
      const double Temp  = plasma->species(s).T0;
    
      const double sub = (plasma->species(s).doGyro) ? 3./2. : 1./2.;
        
      const double v2_rms   = 1.;//pow2(alpha);
      
      for(int m=NmLlD; m<= NmLuD;m++) { 
  
       // gyro-fluid model
  //     if(plasma->species(s).gyroModel == "Gyro-1") k2p_phi(RxLD, RkyLD, RzLD, RFields) = fields->gyroAverage(fields->Field0(RxLD, RkyLD, RzLD, RFields), 2, s,  Field::phi, true);
  //     if(calculate_nonLinear && (rk_step != 0)) calculatePoissonBracket(fields->phi, _fs,m, s);
       
       

       // calculate for estimation of CFL condition
       for(int z=NzLlD; z<= NzLuD;z++) { omp_for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int x=NxLlD; x<= NxLuD;x++) { 
      


     //        const double Krook_nu = krook_nu * ( (X[x] > 0.8 * Lx/2.) || (X[x] < -0.8 * Lx/2.))  ?  0.1 * pow2(abs(X[x]) - 0.8 * Lx/2.): 0.;

             const cmplxd dphi_dx  = (8.*(phi[s][m][z][y_k][x+1] - phi[s][m][z][y_k][x-1]) - (phi[s][m][z][y_k][x+2] - phi[s][m][z][y_k][x-2]))/(12.*dx)  ;  
             const cmplxd phi_     = phi[s][m][z][y_k][x];
             
             const cmplxd ky = cmplxd(0.,fft->ky(y_k));
  //           const cmplxd kp = geo->get_kp(x, ky, z);
             const cmplxd kp = 0.4 * ky * X[x]; //geo->get_kp(x, ky, z);



             // BUG : does it needs a mutex ?
             updateCFL(dphi_dx, ky*phi_, 0.);

            #pragma vector always 
            #pragma ivdep
            for(int v=NvLlD; v<= NvLuD;v++) {
        
        const cmplxd g    = fs [s][m][z][y_k][x][v];
        const cmplxd F0   = vf0[s][m][z][y_k][x][v];

	     //const cmplxd dfs_dv   = (8.  *(fs[s][m][z][y_k][x][v+1] - fs[s][m][z][y_k][x][v-1]) - (fs[s][m][z][y_k][x][v+2] - fs[s][m][z][y_k][x][v-2]))/(12.*dv);
        //const cmplxd ddfs_dvv = (16. *(fs[s][m][z][y_k][x][v+1] + fs[s][m][z][y_k][x][v-1]) - (fs[s][m][z][y_k][x][v+2] + fs[s][m][z][y_k][x][v-2]) - 30.*fs[s][m][z][y_k][x][v])/(12.*pow2(dv));
	     const cmplxd dfs_dv   = (8.  *(fs[s][m][z][y_k][x][v+1] - fs[s][m][z][y_k][x][v-1]) - (fs[s][m][z][y_k][x][v+2] - fs[s][m][z][y_k][x][v-2]))*_12_dv;
        const cmplxd ddfs_dvv = (16. *(fs[s][m][z][y_k][x][v+1] + fs[s][m][z][y_k][x][v-1]) - (fs[s][m][z][y_k][x][v+2] + fs[s][m][z][y_k][x][v-2]) - 30.*fs[s][m][z][y_k][x][v]) * _12_dv_dv;
        
       
        // hyperdiffusion terms
        const cmplxd d4fs_d4x = (-(fs[s][m][z][y_k][x-2][v] + fs[s][m][z][y_k][x+2][v]) + 4. * (fs[s][m][z][y_k][x-1][v] + fs[s][m][z][y_k][x+1][v]) - 6.*fs[s][m][z][y_k][x][v])/_16_dx4;
       
	    /////////////// Finally the Vlasov equation calculate the time derivatve      //////////////////////
       const register cmplxd dg_dt = 
             
	         // driving term (use dphi_dy instead of dXi_dy, because v * A does not vanish due to numerical errors)
             ky* (-(w_n + w_T * (((V[v]*V[v])+ M[m])/Temp  - sub)) * F0 * phi_
             // add first order gyro-average term (zero when full-gyro)
                 )//              + 0.5 * w_T  * k2_phi[1][z][y_k][x] * F0)
      	     // Landau Damping term 
             - alpha  * V[v]* kp  * ( g + sigma * phi_ * F0)
             // collisional term (Lennard-Bernstein)
             + collisionBeta  * (g  + V[v] * dfs_dv + v2_rms * ddfs_dvv)
             // nonLinearTerm
	         + nonLinear[z][y_k][x][v]

             // hyperdiffusive terms
             + hyper_visc[DIR_X] * d4fs_d4x ;//- krook_nu * g
             //+ hyper_visc[DIR_X] * d4fs_d4x ;//- krook_nu * g


        //////////////////////////// Vlasov End ////////////////////////////
        
        //  time-integrate the distribution function    
        ft [s][m][z][y_k][x][v] = rk[0] * ft[s][m][z][y_k][x][v] + rk[1] * dg_dt             ;
        fss[s][m][z][y_k][x][v] = f1[s][m][z][y_k][x][v]         + (rk[2] * ft[s][m][z][y_k][x][v] + dg_dt) * dt;
            
       }}} }
      
//       if(plasma->species(s).gyroModel != "Gyro") break;

      }
   }
}

void VlasovCilk::Vlasov_2D_Island(
                           cmplxd fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd vf0[NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd phi[NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           cmplxd k2_phi[plasma->nfields][NzLD][NkyLD][NxLD],
                           cmplxd nonLinear[NzLD][NkyLD][NxLD][NvLD],
                           cmplxd dphi_dx[NzLB][NkyLD][NxLB],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           Fields *fields,
                           const double dt, const int rk_step, const double rk[3], Array6z _fs)
{ 

    const Geometry2D *geo = static_cast<Geometry2D*>(Vlasov::geo);

    const double w     = setup->get("Island.Width", 0.); 
    const double shear = setup->get("Geometry.Shear", 0.4); 
  

    Xi_max = 0.;

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

       // gyro-fluid model
       if(plasma->species(s).gyroModel == "Gyro-1") k2p_phi(RxLD, RkyLD, RzLD, RFields) = fields->gyroAverage(fields->Field0(RxLD, RkyLD, RzLD, RFields), 2, s,  Field::phi, true);
       
       if(calculate_nonLinear && (rk_step != 0)) calculatePoissonBracket(fields->phi, _fs,m, s);
       //if(calculate_nonLinear && (rk_step != 0)) calculatePoissonBracket(fields->Ap, _fs,m, s);

      for(int z=NzLlD; z<= NzLuD;z++) {  omp_for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) {



            // Note : for negative modes we need to use complex conjugate value

            const cmplxd ky     = cmplxd(0., fft->ky(y_k));

            
            // We need to take care of boundaries. For polidal numbers y_k > N_k-1, we  use zero.
            // For y_k < 0, the corresponding complex conjugate value is used.
            const cmplxd ky_p1  = (y_k == Nky-1) ? 0.                       : cmplxd(0.,fft->ky(y_k+1)) ;
            const cmplxd ky_m1  = (y_k == 0    ) ? cmplxd(0., -fft->ky(1))  : cmplxd(0.,fft->ky(y_k-1)); 
            const cmplxd ky_1   = cmplxd(0., fft->ky(1));

          //#pragma simd
          for(int x=NxLlD; x<= NxLuD;x++) {  

       	     // calculate for estimation of CFL condition
             const cmplxd phi_ = phi[s][m][z][y_k][x];

             dphi_dx[z][y_k][x]  = (8.*(phi[s][m][z][y_k][x+1] - phi[s][m][z][y_k][x-1]) - (phi[s][m][z][y_k][x+2] - phi[s][m][z][y_k][x-2]))/(12.*dx)  ;  

             updateCFL(dphi_dx[z][y_k][x], ky*phi_, 0.);
             
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

        
        const cmplxd     phi_p1 = ( y_k == Nky-1) ? 0.                       : phi[s][m][z][y_k+1][x] ;
        const cmplxd     phi_m1 = ( y_k ==  0   ) ? conj(phi[s][m][z][1][x]) : phi[s][m][z][y_k-1][x] ;


        const cmplxd dphi_dx_p1 = ( y_k == Nky-1) ? 0. 
                                                  : (8.*(phi[s][m][z][y_k+1][x+1] - phi[s][m][z][y_k+1][x-1]) - (phi[s][m][z][y_k+1][x+2] - phi[s][m][z][y_k+1][x-2]))/(12.*dx)  ;

        const cmplxd dphi_dx_m1 = ( y_k ==    0 ) ?  conj(8.*(phi[s][m][z][    1][x+1] - phi[s][m][z][    1][x-1]) - (phi[s][m][z][    1][x+2] - phi[s][m][z][    1][x-2]))/(12.*dx) 
                                                  :      (8.*(phi[s][m][z][y_k-1][x+1] - phi[s][m][z][y_k-1][x-1]) - (phi[s][m][z][y_k-1][x+2] - phi[s][m][z][y_k-1][x-2]))/(12.*dx) ;
        
	    
        // The magnetic island
    
        /*
                *  For the latter term, the intrigate derivative is the \partial_y * Island
                *  remember the island structure is 
                *  \partial_y (e^{imx} + e^{-imx}) = (i m) * ( e^{imx} - e^{-imx} )
                *
        */
        const cmplxd Island_A_phi =   dMagIs_dx * ( ky_m1 * phi_m1 + ky_p1 * phi_p1) - MagIs *  ky_1 * ( dphi_dx_m1 -  dphi_dx_p1);
        
             ///////////////////////////////////////////////////////////////////////////////
            
        const cmplxd kp = geo->get_kp(x, ky, z);
        //const cmplxd kp = geo->shear * ky * X[x];

               // velocity space magic
        #pragma simd
        for(int v=NvLlD; v<= NvLuD;v++) {

            const cmplxd g   =  fs[s][m][z][y_k][x][v];
            const cmplxd F0  = vf0[s][m][z][y_k][x][v];


        /////////////////////////////////////////////////// Magnetic Island Contribution    /////////////////////////////////////////
      
        const cmplxd dfs_dx_p1  =  (y_k == Nky-1) 
                            ? 0.
                            : (8. *(fs[s][m][z][y_k+1][x+1][v] - fs[s][m][z][y_k+1][x-1][v])  - (fs[s][m][z][y_k+1][x+2][v] - fs[s][m][z][y_k+1][x-2][v]))/(12.*dx) ;

        const cmplxd dfs_dx_m1  =  ( y_k == 0   ) 
                            ?  	conj(8. *(fs[s][m][z][    1][x+1][v] - fs[s][m][z][    1][x-1][v])  - (fs[s][m][z][    1][x+2][v] - fs[s][m][z][    1][x-2][v]))/(12.*dx) 
                            :       (8. *(fs[s][m][z][y_k-1][x+1][v] - fs[s][m][z][y_k-1][x-1][v])  - (fs[s][m][z][y_k-1][x+2][v] - fs[s][m][z][y_k-1][x-2][v]))/(12.*dx) ;

        // Note Nky-1 is the maximum mode number Nky = 6 i-> [ 0, 1, 2, 3, 4, 5] 
        const cmplxd fs_p1      = (y_k == Nky-1) ? 0.                         : fs[s][m][z][y_k+1][x][v] ;
        const cmplxd fs_m1      = (y_k ==  0   ) ? conj(fs[s][m][z][1][x][v]) : fs[s][m][z][y_k-1][x][v] ;
         
        // not at the Nyquist frequency we have no coupling with higher frequencies
	
        // mode-mode connections
        register const cmplxd Island_A_F1 =  dMagIs_dx * (ky_m1 * fs_m1  + ky_p1 * fs_p1 )  -  MagIs  * ky_1 *  (dfs_dx_m1  - dfs_dx_p1 )  ;


	
        /////////// Collisions ////////////////////////////////////////////////////////////////////

        const cmplxd dfs_dv    = (8.  *(fs[s][m][z][y_k][x][v+1] - fs[s][m][z][y_k][x][v-1]) - 1. *(fs[s][m][z][y_k][x][v+2] - fs[s][m][z][y_k][x][v-2]))/(12.*dv);
        const cmplxd ddfs_dvv  = (16. *(fs[s][m][z][y_k][x][v+1] + fs[s][m][z][y_k][x][v-1]) - 1. *(fs[s][m][z][y_k][x][v+2] + fs[s][m][z][y_k][x][v-2]) - 30.*fs[s][m][z][y_k][x][v])/(12.*pow2(dv));


        /////////////// Finally the Vlasov equation calculate the time derivatve      //////////////////////
        cmplxd dg_dt = 
             // Island
             - alpha * V[v] * (Island_A_F1 + sigma * Island_A_phi * F0) +   
             // driving term
             ky* (-(w_n + w_T * (((V[v]*V[v])+ M[m])/Temp  - sub)) * F0 * phi_

             // add first order gyro-average term (zero when full-gyro)
             + 0.5 * w_T  * k2_phi[1][z][y_k][x] * F0)
      	     // Landau Damping term 
           - alpha  * V[v]* kp  * ( g + sigma * phi_ * F0);
            // collisional term
            + collisionBeta * (g  + V[v] * dfs_dv + v2_rms * ddfs_dvv)
	        + nonLinear[z][y_k][x][v]
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




void VlasovCilk::Landau_Damping(
                           cmplxd fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd vf0[NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd phi[NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           cmplxd k2_phi[plasma->nfields][NzLD][NkyLD][NxLD],
                           cmplxd dphi_dx[NzLB][NkyLD][NxLB],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           Fields *fields,
                           const double dt, const int rk_step, const double rk[3])
{ 

   for(int s = NsLlD; s <= NsLuD; s++) { for(int m=NmLlD; m<= NmLuD;m++) { 
      
     const double alpha = plasma->species(s).alpha;
  
   for(int z=NzLlD; z<= NzLuD;z++) { for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) {
     
   omp_for(int x=NxLlD; x<= NxLuD;x++) {  for(int v=NvLlD; v<= NvLuD;v++) {
       
           const cmplxd F0          = vf0[s][m][z][y_k][x][v];
           const cmplxd  dphi_dx    = (8.*(phi[s][m][z][y_k][x+1] - phi[s][m][z][y_k][x-1]) - (phi[s][m][z][y_k][x+2] - phi[s][m][z][y_k][x-2]))/(12.*dx)  ;  
           const cmplxd   df_dv     = (8.*(fs[s][m][z][y_k][x][v+1] - fs[s][m][z][y_k][x][v-1]) - (fs[s][m][z][y_k][x][v+2] - fs[s][m][z][y_k][x][v-2]))/(12.*dv);
           const cmplxd   df_dx     = (8. *(fs[s][m][z][y_k][x+1][v] - fs[s][m][z][y_k][x-1][v])  - (fs[s][m][z][y_k][x+2][v] - fs[s][m][z][y_k][x-2][v]))/(12.*dx) ;
           /////////////// 2D (x,v) Landau Damping Test       //////////////////////
        
        cmplxd dg_dt = - alpha * V[v] * (df_dx +  F0 * dphi_dx) + dphi_dx  * df_dv;

        //////////////////////////// Landau End ////////////////////////////

        //  time-integrate the distribution function    
        ft [s][m][z][y_k][x][v] = rk[0] * ft[s][m][z][y_k][x][v] + rk[1] * dg_dt             ;
        fss[s][m][z][y_k][x][v] = f1[s][m][z][y_k][x][v]         + (rk[2] * ft[s][m][z][y_k][x][v] + dg_dt) * dt;

      }}} }}
   }
}

void    VlasovCilk::Vlasov_2D_Fullf(
                           cmplxd fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd vf0[NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd phi[NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           cmplxd k2_phi[plasma->nfields][NzLD][NkyLD][NxLD],
                           Fields *fields,
                           const double dt, const int rk_step, const double rk[3], Array6z _fs)
{ 
   /// Enter Cilk Part here

  Xi_max = 0.;

  Geometry2D *geo = static_cast<Geometry2D*>(Vlasov::geo);

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
       if(plasma->species(s).gyroModel == "Gyro-1") k2p_phi(RxLD, RkyLD, RzLD, RFields) = fields->gyroAverage(fields->Field0(RxLD, RkyLD, RzLD, RFields), 2, s,  Field::phi, true);
       if(calculate_nonLinear && (rk_step != 0)) calculatePoissonBracket(fields->phi, _fs,m, s);
       //if(nonLinear)   calculatePoissonBracket(fields->phi(RxLD, RkyLD, RzLD, m, s), _fs,m, s);
       
       // calculate for estimation of CFL condition
       for(int z=NzLlD; z<= NzLuD;z++) { omp_for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int x=NxLlD; x<= NxLuD;x++) { 
       
             dphi_dx(x,y_k,z)  = (8.*(phi[s][m][z][y_k][x+1] - phi[s][m][z][y_k][x-1]) - (phi[s][m][z][y_k][x+2] - phi[s][m][z][y_k][x-2]))/(12.*dx)  ;  

             const cmplxd ky = cmplxd(0., fft->ky(y_k));
             const cmplxd kp = geo->get_kp(x, ky, z);

             Xi_max(DIR_X) = max(Xi_max(DIR_X), abs(dphi_dx(x,y_k,z)));
             Xi_max(DIR_Y) = max(Xi_max(DIR_Y), abs(ky*phi[s][m][z][y_k][x])); 
             Xi_max(DIR_Z) = max(Xi_max(DIR_Z), abs(kp));// * abs(phi[s][m][z][y_k][x])); 


            for(int v=NvLlD; v<= NvLuD;v++) {
        



        const cmplxd g    = fs[s][m][z][y_k][x][v];
        const cmplxd F0   = vf0[s][m][z][y_k][x][v];
        const cmplxd phi_ = phi[s][m][z][y_k][x];

        // Hyper diffusion terms
        //const cmplxd d4g_dv    =  0.;
//        const cmplxd d4g_dv    =  -1.e-3 * (-39. *(fs[s][m][z][y_k][x][v+1] - fs[s][m][z][y_k][x][v-1])  + 12. *(fs[s][m][z][y_k][x][v+2] - fs[s][m][z][y_k][x][v-2]) + 56. * fs[s][m][z][y_k][x][v]);///pow4(dv);

        const cmplxd dfs_dv    = (8.  *(fs[s][m][z][y_k][x][v+1] - fs[s][m][z][y_k][x][v-1]) - (fs[s][m][z][y_k][x][v+2] - fs[s][m][z][y_k][x][v-2]))/(12.*dv);
        const cmplxd ddfs_dvv  = (16. *(fs[s][m][z][y_k][x][v+1] + fs[s][m][z][y_k][x][v-1]) - (fs[s][m][z][y_k][x][v+2] + fs[s][m][z][y_k][x][v-2]) - 30.*fs[s][m][z][y_k][x][v])/(12.*pow2(dv));
        const double v2_rms = 1.;//pow2(alpha)
        ;
        /////////////// Finally the Vlasov equation calculate the time derivatve      //////////////////////
        cmplxd dg_dt = 


             nonLinearTerms(x,y_k,z,v) +
             
             // driving term (use dphi_dy instead of dXi_dy, because v * A does not vanish due to numerical errors)
             //ky * (-(w_n + w_T * ((pow2(V(v))+ M(m))/Temp  - sub)) * F0 * phi_
             ky * (-(w_T * pow2(V(v)) )* F0 * phi_
                    -(w_n + w_T * (M(m)/Temp  - sub)) * F0 * phi_)

             // add first order gyro-average term (zero when full-gyro)
//             + 0.5 * w_T  * k2_phi[z][y_k][x] * F0

      	     // Landau Damping term and parallel ... ? - alpha  * V(v)* geo->get_kp(x)  * ( g + sigma * phi * F0)) 
           - alpha  * V(v)* kp * ( g + sigma * phi_ * F0 )
           // Collisional terms 
             + collisionBeta * (g  + V(v) * dfs_dv + v2_rms * ddfs_dvv)
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


void VlasovCilk::Vlasov_EM_2D(
                           cmplxd fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd vf0[NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd phi[NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           const cmplxd Ap [NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           const cmplxd Bp [NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           cmplxd    k2_phi[plasma->nfields][NzLD][NkyLD][NxLD],
                           cmplxd   dphi_dx[NzLB][NkyLD][NxLB],
                           cmplxd Xi       [NzLB][NkyLD][NxLB][NvLB],
                           cmplxd G        [NzLB][NkyLD][NxLB][NvLB],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           Fields *fields,
                           const double dt, const int rk_step, const double rk[3])
{ 


  Geometry2D *geo = static_cast<Geometry2D*>(Vlasov::geo);

  Xi_max = 0.;

   for(int s = NsLlD; s <= NsLuD; s++) {
        
      // abbrevations
      const double w_n   = plasma->species(s).w_n;
      const double w_T   = plasma->species(s).w_T;
      const double alpha = plasma->species(s).alpha;
      const double sigma = plasma->species(s).sigma;
      const double Temp  = plasma->species(s).T0;
    
      const double sub = (plasma->species(s).doGyro) ? 3./2. : 1./2.;
      

         for(int m=NmLlD; m<= NmLuD;m++) { 
 
            setupXiAndG(fs, vf0, phi, Ap, Bp, Xi, G, V, M, m , s);

            // gyro-fluid model
           // if(nonLinear)   calculatePoissonBracket(Xi, _fs,m, s);
       
            // calculate for estimation of CFL condition
            for(int z=NzLlD; z<= NzLuD;z++) { omp_for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int x=NxLlD; x<= NxLuD;x++) { 
       
               const cmplxd phi_ = phi[s][m][z][y_k][x];
               dphi_dx[z][y_k][x] = (8.*(phi[s][m][z][y_k][x+1] - phi[s][m][z][y_k][x-1]) - (phi[s][m][z][y_k][x+2] - phi[s][m][z][y_k][x-2]))/(12.*dx)  ;  

               const cmplxd ky = cmplxd(0., fft->ky(y_k));
               const cmplxd kp = geo->get_kp(x, ky, z);

               updateCFL(dphi_dx[z][y_k][x], ky*phi_, 0.);


   #pragma ivdep
   #pragma vector always 
   for(int v=NvLlD; v<= NvLuD;v++) {
        

      const cmplxd g    = fs[s][m][z][y_k][x][v];
      const cmplxd F0   = vf0[s][m][z][y_k][x][v];
      const cmplxd G_   = G[z][y_k][x][v];
      const cmplxd Xi_  = Xi[z][y_k][x][v];

      // Velocity derivaties for Lennard-Bernstein Collisional Model
      const cmplxd dfs_dv   = (8. *(fs[s][m][z][y_k][x][v+1] - fs[s][m][z][y_k][x][v-1]) - (fs[s][m][z][y_k][x][v+2] - fs[s][m][z][y_k][x][v-2]))/(12.*dv);
      const cmplxd ddfs_dvv = (16.*(fs[s][m][z][y_k][x][v+1] + fs[s][m][z][y_k][x][v-1]) - (fs[s][m][z][y_k][x][v+2] + fs[s][m][z][y_k][x][v-2]) - 30.*fs[s][m][z][y_k][x][v])/(12.*dv*dv);
      const double v2_rms = 1.;//pow2(alpha)
    
      cmplxd k2_Xi = 0.;
      if(plasma->species(s).gyroModel == "Gyro-1") { 
        const cmplxd ddXi_dx = (16.*(Xi[z][y_k][x+1][v] + Xi[z][y_k][x-1][v]) - (Xi[z][y_k][x+2][v] + Xi[z][y_k][x-2][v]) - 30.*Xi[z][y_k][x][v])/(12.*dx*dx);
        k2_Xi = (ky*ky + ddXi_dx);
     }
   
     /////////////// Finally the Vlasov equation calculate the time derivatve      //////////////////////
   
     cmplxd dg_dt = 
             
      //nonLinearTerms(x,y_k,z,v) +
             
      // driving term (use dphi_dy instead of dXi_dy, because v * A does not vanish due to numerical errors)
      ky* (-(w_n + w_T * (((V[v]*V[v])+ M[m])/Temp  - sub)) * F0 * Xi_
             
      // add first order gyro-average term (zero when full-gyro)
            + 0.5 * w_T  * k2_Xi * F0)
      // Landau Damping term and parallel ... ? 
            - alpha  * V[v]* kp * G_
      // Collisional terms 
      + collisionBeta * (g  + V[v] * dfs_dv + v2_rms * ddfs_dvv)    ;

        //////////////////////////// Vlasov End ////////////////////////////
  
      //  time-integrate the distribution function    
      ft [s][m][z][y_k][x][v] = rk[0] * ft[s][m][z][y_k][x][v] + rk[1] * dg_dt             ;
      fss[s][m][z][y_k][x][v] = f1[s][m][z][y_k][x][v]         + (rk[2] * ft[s][m][z][y_k][x][v] + dg_dt) * dt;



        }
      }}} }}
   }
}


void VlasovCilk::setupXiAndG(
                           const cmplxd g   [NsLD][NmLD ][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd f0  [NsLD][NmLD ][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd phi [NsLD][NmLD ][NzLB][NkyLD][NxLB+4],
                           const cmplxd Ap  [NsLD][NmLD ][NzLB][NkyLD][NxLB+4],
                           const cmplxd Bp  [NsLD][NmLD ][NzLB][NkyLD][NxLB+4],
                           cmplxd Xi        [NzLB][NkyLD][NxLB][NvLB ],
                           cmplxd G         [NzLB][NkyLD][NxLB][NvLB ],
                           const double V[NvGB], const double M[NmGB],
                           const int m, const int s) {

  // small abbrevations
  const double alpha = plasma->species(s).alpha;
  const double sigma = plasma->species(s).sigma;
  
  const double aeb   =  alpha* geo->eps_hat * plasma->beta; 

  const bool useAp = (plasma->nfields >= 2);
  const bool useBp = (plasma->nfields >= 3);


  // setup values
  for(int z = NzLlB; z <= NzLuB; z++) { omp_for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { 
  for(int x = NxLlB; x <= NxLuB; x++) {     for(int v   = NvLlB ;   v <= NvLuB ;   v++) { 

     Xi[z][y_k][x][v] = phi[s][m][z][y_k][x] - (useAp ? aeb*V[v]*Ap[s][m][z][y_k][x] : 0.) - (useBp ? aeb*M[m]*Bp[s][m][z][y_k][x] : 0.);

     G [z][y_k][x][v] = g[s][m][z][y_k][x][v]  + sigma * Xi[z][y_k][x][v] * f0[s][m][z][y_k][x][v];
             
          /**   f1(x,y,z,v) = g(x,y,z,v,m,s) - 
                          ((plasma->nfields >= 2) ? plasma->species(s).sigma * plasma->species(s).alpha * V(v)*geo->eps_hat 
                             * f0(x,y,z,v,m,s) * plasma->beta * fields->Ap(x,y,z,m,s) : 0.);
                             * */
      }} }}

};


void VlasovCilk::printOn(ostream &output) const
{
	Vlasov::printOn(output);
            output << "Vlasov     |   Spatial : Morinishi        Time : Runge-Kutta 4 " << std::endl;
            output << "Vlasov     |   Hyper Viscosity : " << hyper_visc << std::endl;

};




/*
 *
void VlasovCilk::Vlasov_2D_Island(
                           cmplxd fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd vf0[NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd phi[NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           cmplxd k2_phi[plasma->nfields][NzLD][NkyLD][NxLD],
                           cmplxd nonLinear[NzLD][NkyLD][NxLD][NvLD],
                           cmplxd dphi_dx[NzLB][NkyLD][NxLB],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           Fields *fields,
                           const double dt, const int rk_step, Array6z _fs)
{ 


            const cmplxd ky     = cmplxd(0.,fft->ky(y_k));
            const cmplxd ky_p1  = (y_k+1 <= Nky-1) ? cmplxd(0.,fft->ky(y_k+1)) : 0.;
            const cmplxd ky_m1  = (y_k-1 >= 0    ) ? cmplxd(0.,fft->ky(y_k-1)) : 0.; 
            const cmplxd ky_1    = cmplxd(0., fft->ky(1));

          //#pragma simd
          for(int x=NxLlD; x<= NxLuD;x++) {  

        
       	     // calculate for estimation of CFL condition
             const cmplxd phi_ = phi[s][m][z][y_k][x];
             dphi_dx[z][y_k][x]  = (8.*(phi[s][m][z][y_k][x+1] - phi[s][m][z][y_k][x-1]) - (phi[s][m][z][y_k][x+2] - phi[s][m][z][y_k][x-2]))/(12.*dx)  ;  

             updateCFL(dphi_dx[z][y_k][x], ky*phi_, 0.);
             
        /////////////////////////////////////////////////// Magnetic Island Contribution    /////////////////////////////////////////

        const double zeta=2.*M_PI/Ly;
        const double MagIs     =  0.5 * w*w*shear/16. * cos(zeta * X[x]);
        const double dMagIs_dx = -  0.5 * w*w*shear/16. * sin(zeta * X[x]) * zeta;

        // Note : f, fftw ordering [ 0, 1, 2, 3, 4, -3, -2, -1 ]
        const cmplxd     phi_p1 = (( y_k+1) <= Nky-1) ? phi[s][m][z][y_k+1][x] : 0.;
        const cmplxd     phi_m1 = (( y_k-1) >=    0) ? phi[s][m][z][y_k-1][x] : 0.;
        const cmplxd dphi_dx_p1 = (( y_k+1) <= Nky-1) ?
                    (8.*(phi[s][m][z][y_k+1][x+1] - phi[s][m][z][y_k+1][x-1]) - (phi[s][m][z][y_k+1][x+2] - phi[s][m][z][y_k+1][x-2]))/(12.*dx)  : 0.;
        const cmplxd dphi_dx_m1 = (( y_k-1) >=    0) ?
                    (8.*(phi[s][m][z][y_k-1][x+1] - phi[s][m][z][y_k-1][x-1]) - (phi[s][m][z][y_k-1][x+2] - phi[s][m][z][y_k-1][x-2]))/(12.*dx)  : 0.;
        
        // NOTE :  at the Nyquist frequency we have no coupling with higher frequencies (actually phi(m=Ny) = 0. anyway)
	    
        const cmplxd Island_A_phi =  ( - dMagIs_dx * ( ky_m1 * phi_m1 + ky_p1 * phi_p1) + MagIs * ky_1 * (dphi_dx_m1 - dphi_dx_p1));
        
             ///////////////////////////////////////////////////////////////////////////////
            
        const cmplxd kp = geo->get_kp(x, ky, z);


        /////////////////////////////////////////////////// Magnetic Island Contribution    /////////////////////////////////////////
      
        // how does they couple with Nyquist frequencty ??   // take care of boundaries !!
        const cmplxd dfs_dx_p1  =  ((y_k+1) <= Nky-1) ?
                              	   (8. *(fs[s][m][z][y_k+1][x+1][v] - fs[s][m][z][y_k+1][x-1][v])  - (fs[s][m][z][y_k+1][x+2][v] - fs[s][m][z][y_k+1][x-2][v]))/(12.*dx) : 0.;

        const cmplxd dfs_dx_m1  =  ((y_k-1) >= 0   ) ? 
                              	   (8. *(fs[s][m][z][y_k-1][x+1][v] - fs[s][m][z][y_k-1][x-1][v])  - (fs[s][m][z][y_k-1][x+2][v] - fs[s][m][z][y_k-1][x-2][v]))/(12.*dx) : 0.;

        const cmplxd fs_p1      = ((y_k+1) <= Nky-1) ? fs[s][m][z][y_k+1][x][v] : 0.;
        const cmplxd fs_m1      = ((y_k-1) >= 0   ) ? fs[s][m][z][y_k-1][x][v] :  0.;
         

        // not at the Nyquist frequency we have no coupling with higher frequencies
	
        // mode-mode connections
        register const cmplxd Island_A_F1 = ( - dMagIs_dx * (ky_m1 *fs_m1  + ky_p1 * fs_p1 )  +  MagIs  * ky_1 *  (dfs_dx_m1  - dfs_dx_p1 ))  ;
	
        /////////////// Finally the Vlasov equation calculate the time derivatve      //////////////////////
        cmplxd dg_dt = 
             // Island
             alpha * V[v] * (Island_A_F1 + Island_A_phi * F0) +   
*/


void VlasovCilk::initDataOutput(FileIO *fileIO) {
                Vlasov::initDataOutput(fileIO); 

                //check(H5LTset_attribute_string(vlasovGroup, ".", "CollisionModel","LennardBernstein"), DMESG("H5LTset_attribute"));
               //check(H5LTset_attribute_double(vlasovGroup, ".", "CollisionBeta"   ,  &collisionBeta, 1), DMESG("H5LTset_attribute"));
        
};


void VlasovCilk::Vlasov_EM(
                           cmplxd fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd vf0[NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd phi[NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           const cmplxd Ap [NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           const cmplxd Bp [NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           cmplxd    k2_phi[plasma->nfields][NzLD][NkyLD][NxLD],
                           cmplxd   dphi_dx[NzLB][NkyLD][NxLB],
                           cmplxd Xi       [NzLB][NkyLD][NxLB][NvLB],
                           cmplxd G        [NzLB][NkyLD][NxLB][NvLB],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           Fields *fields,
                           const double dt, const int rk_step, const double rk[3])
{ 

   
   const double B0 = plasma->B0;

   for(int s = NsLlD; s <= NsLuD; s++) {
        
      // small abbrevations
      const double w_n   = plasma->species(s).w_n;
      const double w_T   = plasma->species(s).w_T;
      const double alpha = plasma->species(s).alpha;
      const double sigma = plasma->species(s).sigma;
      const double Temp  = plasma->species(s).T0;
    
      const double sub = (plasma->species(s).doGyro) ? 3./2. : 1./2.;
      

      for(int m=NmLlD; m<= NmLuD;m++) { 
 
         setupXiAndG(fs, vf0, phi, Ap, Bp, Xi, G, V, M, m , s);

      //   if(nonLinear)   calculatePoissonBracket(Xi, _fs,m, s);
       
         // calculate for estimation of CFL condition
         for(int z=NzLlD; z<= NzLuD;z++) { omp_for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int x=NxLlD; x<= NxLuD;x++) { 
       
            const cmplxd phi_ = phi[s][m][z][y_k][x];

            dphi_dx[z][y_k][x] = (8.*(phi[s][m][z][y_k][x+1] - phi[s][m][z][y_k][x-1]) - (phi[s][m][z][y_k][x+2] - phi[s][m][z][y_k][x-2]))/(12.*dx)  ;  

            const cmplxd ky = cmplxd(0.,fft->ky(y_k));

             updateCFL(dphi_dx[z][y_k][x], ky*phi_, 0.);


     
      for(int v=NvLlD; v<= NvLuD;v++) {
        

          
        const cmplxd g    = fs[s][m][z][y_k][x][v];
        const cmplxd F0   = vf0[s][m][z][y_k][x][v];
        const cmplxd G_   = G[z][y_k][x][v];
        const cmplxd Xi_  = Xi[z][y_k][x][v];

        
        // Velocity derivaties for Lennard-Bernstein Collisional Model
        const cmplxd dfs_dv   = (8. *(fs[s][m][z][y_k][x][v+1] - fs[s][m][z][y_k][x][v-1]) - (fs[s][m][z][y_k][x][v+2] - fs[s][m][z][y_k][x][v-2]))/(12.*dv);
        const cmplxd ddfs_dvv = (16.*(fs[s][m][z][y_k][x][v+1] + fs[s][m][z][y_k][x][v-1]) - (fs[s][m][z][y_k][x][v+2] + fs[s][m][z][y_k][x][v-2]) - 30.*fs[s][m][z][y_k][x][v])/(12.*dv*dv);
        const double v2_rms = 1.;//pow2(alpha)
    
        
        /////////////// Finally the Vlasov equation calculate the time derivatve      //////////////////////

        // We use CD-4 (central difference fourth order for every variable)

        const cmplxd dG_dz   = (8.*(G[z+1][y_k][x][v] - G[z-1][y_k][x][v])    -1.*(G[z+2][y_k][x][v] - G[z-2][y_k][x][v]))/(12.*dz);
        const cmplxd dG_dx   = (8.*(G[z][y_k][x+1][v] - G[z][y_k][x-1][v])    -1.*(G[z][y_k][x+2][v] - G[z][y_k][x-2][v]))/(12.*dx);

        
        // magnetic prefactor defined as  $ \hat{B}_0 / \hat{B}_{0\parallel}^\star = \left[ 1 + \beta_{ref} \sqrt{\frac{\hat{m_\sigma T_{0\sigma}{2}}}}
        // note j0 is calculated and needs to be replaced, or ? no we calculate j1 ne ?!
        const double j0 = 0.;
        const double Bpre  = 1.; //1./(1. + plasma->beta * sqrt(m * T/2.) * j0 / (q * pow2(geo->B(x,y,z))) * V(v));
        const double CoJB = 1./geo->J(x,z);

        
        // Finally the Vlasov equation calculate the time derivatve      
        
        cmplxd dg_dt = 
            
          // driving term
          Bpre * (w_n + w_T * ((pow2(V[v])+ M[m] * B0)/Temp - sub)) * F0 * Xi_ * ky
          - Bpre * sigma * ((M[m] * B0 + 2.*pow2(V[v]))/B0) * (geo->Kx(x,z) * dG_dx - geo->Ky(x,z) * ky * G_)
          //- alpha * pow2(V[v]) * plasma->beta * plasma->w_p * G_ * ky
          -  CoJB *  alpha * V[v]* dG_dz
          + alpha  / 2. * M[m] * geo->dB_dz(x,z) * dfs_dv
          + Bpre *  sigma * (M[m] * B0 + 2. * pow2(V[v]))/B0 * geo->Kx(x,z) * ((w_n + w_T * (pow2(V[v]) + M[m] * B0)/Temp - sub) * dG_dx + sigma * dphi_dx[z][y_k][x] * F0);

          
        //////////////////////////// Vlasov End ////////////////////////////

        //  time-integrate the distribution function    
        ft [s][m][z][y_k][x][v] = rk[0] * ft[s][m][z][y_k][x][v] + rk[1] * dg_dt             ;
        fss[s][m][z][y_k][x][v] = f1[s][m][z][y_k][x][v]         + (rk[2] * ft[s][m][z][y_k][x][v] + dg_dt) * dt;
        
      }}} }}
   }
}


