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


#ifndef __cilk
#include <cilk/cilk_stub.h>
#endif
#include <cilk/cilk.h>


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

VlasovCilk::VlasovCilk(Grid *_grid, Parallel *_parallel, Setup *_setup, FileIO *fileIO, Geometry<HELIOS_GEOMETRY> *_geo, FFTSolver *fft)    
    : Vlasov(_grid, _parallel, _setup, fileIO, _geo, fft),
      dphi_dx(FortranArray<3>()),   dphi_dy(FortranArray<3>()),   
      dAp_dx(FortranArray<3>()),    dAp_dy(FortranArray<3>()),    k2p_phi(FortranArray<4>()), nonLinearTerms(HeliosStorage4)
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


int VlasovCilk::solve(std::string equation_type, Fields *fields, Array6z _fs, Array6z _fss, double dt, int rk_step, int user_boundary_type) 
{
   /// Enter Cilk Part here
  Xi_max = 0.;
  
  if((equation_type == "2D_ES")) Vlasov_2D((A6z) _fs.dataZero(), (A6z) _fss.dataZero(), (A6z) f0.dataZero(), (A6z) f.dataZero(), (A6z) ft.dataZero(), (A5z) fields->phi.dataZero(), (A4z) k2p_phi.dataZero(), (A4z) nonLinearTerms.dataZero(), X.dataZero(), V.dataZero(), M.dataZero(), fields, dt, rk_step, _fs);
  else if((equation_type == "2D_EM")) Vlasov_EM((A6z) _fs.dataZero(), (A6z) _fss.dataZero(), (A6z) f0.dataZero(), (A6z) f.dataZero(), (A6z) ft.dataZero(), (A5z) fields->phi.dataZero(), (A5z) fields->Ap.dataZero(), (A5z) fields->Bp.dataZero(), (A4z) k2p_phi.dataZero(), (A3z) dphi_dx.dataZero(), (A4z) Xi.dataZero(), (A4z) G.dataZero(), X.dataZero(), V.dataZero(), M.dataZero(), fields, dt, rk_step);
  else if((equation_type == "2DIsland")) Vlasov_2D_Island((A6z) _fs.dataZero(), (A6z) _fss.dataZero(), (A6z) f0.dataZero(), (A6z) f.dataZero(), (A6z) ft.dataZero(), (A5z) fields->phi.dataZero(), (A4z) k2p_phi.dataZero(), (A4z) nonLinearTerms.dataZero(), (A3z) dphi_dx.dataZero(), 
      X.dataZero(), V.dataZero(), M.dataZero(), fields, dt, rk_step, _fs);
  else if((equation_type == "2DLandauDamping")) Landau_Damping((A6z) _fs.dataZero(), (A6z) _fss.dataZero(), (A6z) f0.dataZero(), (A6z) f.dataZero(), (A6z) ft.dataZero(), (A5z) fields->phi.dataZero(), (A4z) k2p_phi.dataZero(), (A3z) dphi_dx.dataZero(), 
      X.dataZero(), V.dataZero(), M.dataZero(), fields, dt, rk_step);
  else   check(-1, DMESG("No Such Equation"));

  return HELIOS_SUCCESS;
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



Array4z VlasovCilk::calculatePoissonBracket(Array5z phi, Array6z f1, const int m, const int s) 
{

//   xy_phi(RxLD, RyLD, RzLD) = transform2Real(phi, m , s);
   transform2Real(phi, m , s); xy_phi(RxLD, RyLD, RzLD) = fft->rYOut(RxLD, RyLD, RzLD, Field::phi);
   setBoundaryXY(xy_phi,DIR_XY);

   // perform CD-4 derivative for dphi_dx , and dphi_dy
   for(int z = NzLlD; z <= NzLuD; z++) { for(int y=NyLlD; y<= NyLuD;y++) { for(int x=NxLlD; x<= NxLuD;x++)  {

     xy_dphi_dx(x, y, z) = (8.*(xy_phi(x+1, y,z) - xy_phi(x-1, y, z)) - (xy_phi(x+2,y,z) - xy_phi(x-2,y,z)))/(12.*dx); 
     xy_dphi_dy(x, y, z) = (8.*(xy_phi(x, y+1,z) - xy_phi(x, y-1, z)) - (xy_phi(x,y+2,z) - xy_phi(x,y-2,z)))/(12.*dy); 

   } } }

   setBoundaryXY(xy_dphi_dx,DIR_Y);
   setBoundaryXY(xy_dphi_dy,DIR_X);
   
   
   // phase space function & Poisson bracket
   const double fft_Norm = fft->Norm_Y_Backward * fft->Norm_Y_Backward * fft->Norm_Y_Forward;
   for(int v=NvLlD; v<=NvLuD;v++) { 
   
//      xy_f1(RxLD, RyLD, RzLD) = transform2Real(f1, v, m ,s);
      transform2Real(f1, v, m ,s); xy_f1(RxLD, RyLD, RzLD) = fft->rYOut(RxLD, RyLD, RzLD, Field::phi);
     setBoundaryXY(xy_f1,DIR_XY);
     
      
      ///////////////////////   calculate cross terms using Morinishi scheme    [ phi, F1]  //////////////////////////////////
      for(int z=NzLlD; z<=NzLuD;z++) { for(int y=NyLlD; y<=NyLuD;y++) { for(int x= NxLlD; x <= NxLuD; x++) {

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
                           cmplxd fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd vf0[NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd phi[NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           cmplxd k2_phi[plasma->nfields][NzLD][NkyLD][NxLD],
                           cmplxd nonLinear[NzLD][NkyLD][NxLD][NvLD],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           Fields *fields,
                           const double dt, const int rk_step, Array6z _fs)
{ 


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
       if(plasma->species(s).gyroModel == "Gyro-1") k2p_phi(RxLD, RkyLD, RzLD, RFields) = fields->gyroAverage(fields->Field0(RxLD, RkyLD, RzLD, RFields), 2, s,  Field::phi, true);
       if(calculate_nonLinear && (rk_step != 0)) calculatePoissonBracket(fields->phi, _fs,m, s);
       
       

       // calculate for estimation of CFL condition
       for(int z=NzLlD; z<= NzLuD;z++) {  omp_for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int x=NxLlD; x<= NxLuD;x++) { 
      


             const double Krook_nu = krook_nu * ( (X[x] > 0.8 * Lx/2.) || (X[x] < -0.8 * Lx/2.))  ?  0.1 * pow2(abs(X[x]) - 0.8 * Lx/2.): 0.;

             const cmplxd dphi_dx  = (8.*(phi[s][m][z][y_k][x+1] - phi[s][m][z][y_k][x-1]) - (phi[s][m][z][y_k][x+2] - phi[s][m][z][y_k][x-2]))/(12.*dx)  ;  
             const cmplxd phi_     = phi[s][m][z][y_k][x];
             
             const cmplxd ky = cmplxd(0.,fft->ky(y_k));
             const cmplxd kp = geo->get_kp(x, ky, z);



             // BUG : does it needs a mutex ?
             updateCFL(dphi_dx, ky*phi_, 0.);

            
            #pragma simd 
            for(int v=NvLlD; v<= NvLuD;v++) {
        
        const cmplxd g    = fs [s][m][z][y_k][x][v];
        const cmplxd F0   = vf0[s][m][z][y_k][x][v];

	     const cmplxd dfs_dv   = (8.  *(fs[s][m][z][y_k][x][v+1] - fs[s][m][z][y_k][x][v-1]) - (fs[s][m][z][y_k][x][v+2] - fs[s][m][z][y_k][x][v-2]))/(12.*dv);
        const cmplxd ddfs_dvv = (16. *(fs[s][m][z][y_k][x][v+1] + fs[s][m][z][y_k][x][v-1]) - (fs[s][m][z][y_k][x][v+2] + fs[s][m][z][y_k][x][v-2]) - 30.*fs[s][m][z][y_k][x][v])/(12.*pow2(dv));
        
       
        // hyperdiffusion terms
        const cmplxd d4fs_d4x = (-(fs[s][m][z][y_k][x-2][v] + fs[s][m][z][y_k][x+2][v]) + 4. * (fs[s][m][z][y_k][x-1][v] + fs[s][m][z][y_k][x+1][v]) - 6.*fs[s][m][z][y_k][x][v])/(16.*pow4(dx));
       
	    /////////////// Finally the Vlasov equation calculate the time derivatve      //////////////////////
        cmplxd dg_dt = 
             
	         // driving term (use dphi_dy instead of dXi_dy, because v * A does not vanish due to numerical errors)
             ky* (-(w_n + w_T * (((V[v]*V[v])+ M[m])/Temp  - sub)) * F0 * phi_
             // add first order gyro-average term (zero when full-gyro)
                  + 0.5 * w_T  * k2_phi[1][z][y_k][x] * F0)
      	     // Landau Damping term 
             - alpha  * V[v]* kp  * ( g + sigma * phi_ * F0)
             // collisional term (Lennard-Bernstein)
             + collisionBeta  * (g  + V[v] * dfs_dv + v2_rms * ddfs_dvv)
             // nonLinearTerm
	         + nonLinear[z][y_k][x][v]

             // hyperdiffusive terms
             + hyper_visc[DIR_X] * d4fs_d4x - krook_nu * g
          ;


        //////////////////////////// Vlasov End ////////////////////////////
        if(rk_step == 0) {
            fss[s][m][z][y_k][x][v] = dg_dt;
        } else {

            //  time-integrate the distribution function     
            if(rk_step == 1) ft[s][m][z][y_k][x][v] = dg_dt;
            else if((rk_step == 2) || (rk_step == 3)) ft[s][m][z][y_k][x][v] = ft[s][m][z][y_k][x][v] + 2.0*dg_dt;
            else    dg_dt = ft[s][m][z][y_k][x][v] + dg_dt;
        
            fss[s][m][z][y_k][x][v] = f1[s][m][z][y_k][x][v] + dg_dt*dt;
       }  

      }}} }
      
//       if(plasma->species(s).gyroModel != "Gyro") break;

      }
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
                           const double dt, const int rk_step)
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
        if(rk_step == 1) ft[s][m][z][y_k][x][v] = dg_dt;
        else if((rk_step == 2) || (rk_step == 3)) ft[s][m][z][y_k][x][v] = ft[s][m][z][y_k][x][v] + 2.0*dg_dt;
        else    dg_dt = ft[s][m][z][y_k][x][v] + dg_dt;
        
        fss[s][m][z][y_k][x][v] = f1[s][m][z][y_k][x][v] + dg_dt*dt;


      }}} }}
   }
}

void    VlasovCilk::Vlasov_2D_Global(
                           cmplxd fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd vf0[NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd phi[NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           cmplxd k2_phi[plasma->nfields][NzLD][NkyLD][NxLD],
                           Fields *fields,
                           const double dt, const int rk_step, Array6z _fs)
{ 
   /// Enter Cilk Part here

  Xi_max = 0.;

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
        if(rk_step == 0) {
            fss[s][m][z][y_k][x][v] = dg_dt;
        } else {

            if(rk_step == 1) ft[s][m][z][y_k][x][v] = dg_dt;
            else if((rk_step == 2) || (rk_step == 3)) ft[s][m][z][y_k][x][v] = ft[s][m][z][y_k][x][v] + 2.0*dg_dt;
            else    dg_dt = ft[s][m][z][y_k][x][v] + dg_dt;
        
            fss[s][m][z][y_k][x][v] = f1[s][m][z][y_k][x][v] + dg_dt*dt;
        }

      }}} }}
   }
//if(rk_step == 4)        std::cout << sum(pow2(fields->phi(RxLD, RkyLD, RzLD, 1, 1))) << " nl : " << sum(pow2(nonLinearTerms(RxLD, RkyLD, RzLD, RvLD))) << std::endl; 
}


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
                           const double dt, const int rk_step)
{ 

   

   Xi_max = 0.;


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

         // gyro-fluid model
       //if(nonLinear)   calculatePoissonBracket(Xi, _fs,m, s);
       
         // calculate for estimation of CFL condition
         for(int z=NzLlD; z<= NzLuD;z++) { omp_for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int x=NxLlD; x<= NxLuD;x++) { 
       
            const cmplxd phi_ = phi[s][m][z][y_k][x];
            dphi_dx[z][y_k][x] = (8.*(phi[s][m][z][y_k][x+1] - phi[s][m][z][y_k][x-1]) - (phi[s][m][z][y_k][x+2] - phi[s][m][z][y_k][x-2]))/(12.*dx)  ;  

            const cmplxd ky = cmplxd(0.,fft->ky(y_k));
            const cmplxd kp = geo->get_kp(x, ky, z);

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
        if(rk_step == 0) {
            fss[s][m][z][y_k][x][v] = dg_dt;
        } else {
   if(rk_step == 1) ft[s][m][z][y_k][x][v] = dg_dt;
   else if((rk_step == 2) || (rk_step == 3)) ft[s][m][z][y_k][x][v] = ft[s][m][z][y_k][x][v] + 2.0*dg_dt;
   else    dg_dt = ft[s][m][z][y_k][x][v] + dg_dt;
        
   fss[s][m][z][y_k][x][v] = f1[s][m][z][y_k][x][v] + dg_dt*dt;

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

  for(int z = NzLlB; z <= NzLuB; z++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { 
  for(int x = NxLlB; x <= NxLuB; x++) { for(int v   = NvLlB ;   v <= NvLuB ;   v++) { 

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


