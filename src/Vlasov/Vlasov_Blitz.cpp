#include "config.h"
#include "Vlasov_Blitz.h"


VlasovBlitz::VlasovBlitz(Grid *_grid, Parallel *_parallel, Setup *_setup, FileIO *fileIO, Geometry<HELIOS_GEOMETRY> *_geo, FFTSolver *fft)    : Vlasov_LenardBernstein(_grid, _parallel, _setup,fileIO,  _geo, fft),
    dphi_dx(FortranArray<3>()),    dphi_dy(FortranArray<3>()),    dAp_dx(FortranArray<3>()),    dAp_dy(FortranArray<3>()),    k2p_phi(FortranArray<3>())
{


           f1.resize(RxLB , RkyLD , RzLB, RvLB, RmLD, RsLD);  f1 = 0.e0;
          

           if(useAntiAliasing)  {
           xy_df1_dx.resize(RxLD , fft->Y_RyLD , RzLD);  xy_df1_dx = 0.e0;
           xy_df1_dy.resize(RxLD , fft->Y_RyLD , RzLD);  xy_df1_dy = 0.e0;

           xy_dphi_dx.resize(RxLD , fft->Y_RyLD , RzLD);  xy_dphi_dx = 0.e0;
           xy_dphi_dy.resize(RxLD , fft->Y_RyLD , RzLD);  xy_dphi_dy = 0.e0;
           

           } else {

           xy_f1.resize(RxLB , RyLB , RzLD);  xy_df1_dx = 0.e0;
           xy_df1_dx.resize(RxLB , RyLB , RzLD);  xy_df1_dx = 0.e0;
           xy_df1_dy.resize(RxLB , RyLB , RzLD);  xy_df1_dy = 0.e0;

           xy_phi.resize(RxLB , RyLB , RzLD);  xy_dphi_dx = 0.e0;
           xy_dphi_dx.resize(RxLB , RyLB , RzLD);  xy_dphi_dx = 0.e0;
           xy_dphi_dy.resize(RxLB , RyLB , RzLD);  xy_dphi_dy = 0.e0;
           
            }
           nonLinearTerms.resize(RxLD , RkyLD , RzLD, RvLD);  nonLinearTerms = 0.e0;
//        if(setup->get("Vlasov.Equation", "EM_ToroidalNonLinear") == "EM_2DLinear")  { 
		allocate(RxLB, RkyLD, RzLB, dphi_dy, dphi_dx, dAp_dy, dAp_dx);
		
		k2p_phi.resize(RxLD, RkyLD, RzLD) ; k2p_phi = 0.;
//	}
}

int VlasovBlitz::solve(std::string equation_tyoe, Fields *fields, Array6z _fs, Array6z _fss, double dt, int rk_step, int user_boundary_type) 
{

  if((equation_type == "EM_2DLinear")  && (plasma->nfields == 1)) Vlasov_2D          (_fs,_fss, ft, fields, dt, rk_step);
  else if((equation_type == "EM_2DLinearIsland")) Vlasov_2D_Island          (_fs,_fss, ft, fields, dt, rk_step);
  else   check(-1, DMESG("No Such Equation"));

  return HELIOS_SUCCESS;
}

void VlasovBlitz::setupXiAndG(Array6z g, Fields *fields) {

   for(int s = NsLlD; s <= NsLuD; s++) {
        
      // small abbrevations
      const double w_n = plasma->species(s).w_n;
      const double w_T = plasma->species(s).w_T;

      const double alpha = plasma->species(s).alpha;
      const double sigma = plasma->species(s).sigma;

      const double mass = plasma->species(s).m;
      const double Temp = plasma->species(s).T0;
      const double q_ch = plasma->species(s).q;
    
      const double B = plasma->B0;

      for(int m=NmLlB; m<= NmLuB;m++)  { for(int v=NvLlB; v<= NvLuB;v++) { for(int z=NzLlB; z<= NzLuB;z++){

         for(int y=NyLlB-2; y<= NyLuB+2;y++) { for(int x=NxLlB-2; x<= NxLuB+2;x++) {

             Xi(x,y,z,v,m,s) = fields->phi(x,y,z,m,s) - ((plasma->nfields >= 2) ? alpha* geo->eps_hat * plasma->beta * V(v)*fields->Ap(x,y,z,m,s) : 0.)
                                                      - ((plasma->nfields >= 3) ? alpha* geo->eps_hat * plasma->beta * M(m)*fields->Bp(x,y,z,m,s) : 0.);

         }}

          for(int y=NyLlB; y<= NyLuB;y++) { for(int x=NxLlB; x<= NxLuB;x++) { 

             G (x,y,z,v,m,s) = g(x,y,z,v,m,s) + sigma * Xi(x,y,z,v,m,s) * f0(x,y,z,v,m,s);
             f1(x,y,z,v,m,s) = g(x,y,z,v,m,s) - 
                          ((plasma->nfields >= 2) ? plasma->species(s).sigma * plasma->species(s).alpha * V(v)*geo->eps_hat 
                             * f0(x,y,z,v,m,s) * plasma->beta * fields->Ap(x,y,z,m,s) : 0.);
          }} 
      }}}

   }
};


void VlasovBlitz::Vlasov_2D(Array6z fs, Array6z fss, Array6z ft, Fields *fields, double dt, int rk_step) {

  const double w     = setup->get("Island.Width", 0.); 
  const double shear = setup->get("Geometry.Shear", 0.); 
  const bool useDorland1 = (fields->gyroAverageModel == "Dorland1");
  
  const double eps_v = setup->get("Vlasov.EpsV", 0.);

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
       if(useDorland1) k2p_phi(RxLD, RkyLD, RzLD) = fields->gyroAverage(fields->Field0(RxLD, RkyLD, RzLD, Field::phi), 2, s,  Field::phi, true);


       if(nonLinear) {
        nonLinearTerms = (useAntiAliasing) ? calculatePhiNonLinearityAA(fields->phi, fs, m, s): calculatePhiNonLinearity(fields->phi, fs, m, s);

       }

       // calculate for estimation of CFL condition
       for(int z=NzLlD; z<= NzLuD;z++) { 
        #pragma omp parallel for
        for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int x=NxLlD; x<= NxLuD;x++) { 
        
             dphi_dx(x,y_k,z)  = (8.*(fields->phi(x+1, y_k  , z, m, s) - fields->phi(x-1  , y_k, z, m, s)) - (fields->phi(x+2,   y_k, z, m, s) - fields->phi(x-2, y_k  , z, m, s)))/(12.*dx)  ;  

             const cmplxd d_ky = cmplxd(0.,-fft->ky(y_k));

             Xi_max(DIR_X) = max(Xi_max(DIR_X), abs(dphi_dx(x,y_k,z)));
             Xi_max(DIR_Y) = max(Xi_max(DIR_Y), abs(d_ky*fields->phi(x,y_k,z,m,s))); 
             Xi_max(DIR_Z) = 0.;


        }}}

      for(int v=NvLlD; v<= NvLuD;v++) {

      for(int z=NzLlD; z<= NzLuD;z++) { 
        #pragma omp parallel for
        for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int x=NxLlD; x<= NxLuD;x++) {  
        
                      
        const cmplxd d_ky  = cmplxd(0.,-fft->ky(y_k));

        const cmplxd g   = fs(x, y_k, z, v, m, s);
        const cmplxd F0  = f0(x, y_k, z, v, m, s);
        const cmplxd phi = fields->phi(x,y_k, z, m, s); 

        // Hyper diffusion terms
        const cmplxd d4g_dv    =  (-39. *(fs(x, y_k, z, v+1, m, s) - fs(x, y_k, z, v-1, m, s))  + 12. *(fs(x, y_k, z, v+2, m, s) - fs(x, y_k, z, v-2, m, s)) + 56. * fs(x,y_k,z, v, m, s));///pow4(dv);


        // collisions
        /*  
        const double df1_dv    = (8. *(f1(x, y, z, v+1, m, s) - f1(x, y, z, v-1, m, s))  -1. *(f1(x, y, z, v+2, m, s) - f1(x, y, z, v-2, m, s)))/(12.*dv);
        const cmplxd ddfs_dvv  = (16. *(fs(x, y_k, z, v+1, m, s) + fs(x, y_k, z, v-1, m, s))  -1. *(fs(x, y_k, z, v+2, m, s) + fs(x, y_k, z, v-2, m, s)) - 30.*fs(x,y_k,z,v,m,s))/(12.*pow2(dv));
        const double beta=2.e-3;
        const cmplxd v2_rms = 1.;//pow2(plasma->species(s).scale_v);
        */
        // finite difference stencil [ -8 1 0 -1 8 ]

        /////////////// Finally the Vlasov equation calculate the time derivatve      //////////////////////
 
        cmplxd dg_dt = 
          //- Nonlinearity (in case of Linear simulation term is zero
//          -  dXi_dx__dG_dy - dXi_dy__dG_dx
          // include electro-magnetic contributions
  //         + (((plasma->nfields >= 2) && !nonLinear) ? alpha * plasma->beta  * V(v) * (dA_dx__dfs_dy - dA_dy__dfs_dx) : 0.)
         nonLinearTerms(x,y_k,z,v)   + 
             // driving term (use dphi_dy instead of dXi_dy, because v * A does not vanish due to numerical errors)
             d_ky * (-(w_n + w_T * ((pow2(V(v))+ M(m))/Temp  - sub)) * F0 * phi

             // add first order gyro-average term (zero when full-gyro)
             + 0.5 * w_T  * k2p_phi(x,y_k,z) * F0

      	     // Landau Damping term and parallel ... ? - alpha  * V(v)* geo->getShear(x)  * ( g + sigma * phi * F0)) 
           - alpha  * V(v)* geo->getShear(x)  * ( g + sigma * phi * F0))
   //      + beta * (fs(x,y_k,z,v,m,s)  + V(v) * dfs_dv + v2_rms * ddfs_dvv)
          ;

          // simple pitch angle scattering in n without
          //+ Coll(x,y,z,v,m,s);
          //
            // Energy evolution term
//          + rhoOverLn * 1./mass*(shear * X(x) + theta) * dXi_dy *df1_dv;


        //////////////////////////// Vlasov End ////////////////////////////

        //  time-integrate the distribution function     
        if(rk_step == 1) ft(x, y_k, z, v, m, s) = dg_dt;
        else if((rk_step == 2) || (rk_step == 3)) ft(x, y_k, z, v, m, s) = ft(x, y_k, z, v, m, s) + 2.0*dg_dt;
        else    dg_dt = ft(x, y_k, z, v, m, s) + dg_dt;
        
        fss(x, y_k, z, v, m, s) = f(x, y_k, z, v, m, s) + dg_dt*dt;


      }}} }}
   }

}

void VlasovBlitz::printOn(ostream &output) const
{
	Vlasov::printOn(output);
            output << "Vlasov     |   Spatial : Morinishi        Time : Runge-Kutta 4 " << std::endl;
            output << "Vlasov     |   Hyper Viscosity : " << hyper_visc << std::endl;

};


int mapAAY(const int y_k) {
    int N = (Ny*3)/2;
    return (y_k <= Ny/2) ? y_k : N - (Ny-y_k);
};


void copyNyquiest(Array4z f, FFTSolver *fft) {
  int N = (Ny*3)/2;
  f(RxLD, N-Ny/2, RzLD, Field::phi) = -f(RxLD, Ny/2, RzLD, Field::phi); 


}


void setBoundary3(Array3z A) {

        A(Range(NxLlB, NxLlB+1), RyLD, RzLD) = A(Range(NxLuD-1, NxLuD),   RyLD, RzLD);
        A(Range(NxLuD+1, NxLuB), RyLD, RzLD) = A(Range(NxLlD, NxLlD+1),   RyLD, RzLD);
    
        A(RxLD, Range(NyLlB, NyLlB+1), RzLD) = A(RxLD, Range(NyLuD-1, NyLuD), RzLD);
        A(RxLD, Range(NyLuD+1, NyLuB), RzLD) = A(RxLD, Range(NyLlD, NyLlD+1), RzLD);

};

Array4z VlasovBlitz::calculatePhiNonLinearity(Array5z phi, Array6z fs, const int m, const int s) 
{
     
  

   // back transform dphi_dx(x,ky,z,m,s) -> dphi_dx(x,y,z,m,s) 
   fft->kYIn = 0.;
   for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) {
   for(int x=NxLlD; x<= NxLuD;x++)  fft->kYIn(x, y_k, RzLD, Field::phi) = (8.*(phi(x+1, y_k,RzLD, m,s) - phi(x-1, y_k, RzLD,m,s)) - (phi(x+2,   y_k, RzLD,m,s) - phi(x-2, y_k  , RzLD,m,s)))/(12.*dx); 
   }
   xy_dphi_dx(RxLD,RyLD,RzLD) = (fft->rYOut(RxLD, RyLD, RzLD, Field::phi));
   setBoundary3(xy_dphi_dx);
  

   // dphi_dy
   fft->kYIn = 0.;
   for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) fft->kYIn(RxLD, y_k, RzLD, Field::phi) = phi(RxLD, y_k, RzLD, m, s) ;// * cmplxd(0.,  -fft->ky(y_k));
   fft->solve(FFT_Y, FFT_BACKWARD, NxLD * NzLD);

   xy_dphi_dy(RxLD,RyLD,RzLD) = (fft->rYOut(RxLD, RyLD, RzLD, Field::phi));
   setBoundary3(xy_dphi_dy);
   
   for(int y=NyLlD; y<= NyLuD;y++)  xy_dphi_dy(RxLD, y, RzLD) = (8.*(xy_dphi_dy(RxLD, y+1,RzLD) - xy_dphi_dy(RxLD, y-1, RzLD)) - (xy_dphi_dy(RxLD,   y+2, RzLD) - xy_dphi_dy(RxLD, y-2  , RzLD)))/(12.*dy); 
   setBoundary3(xy_dphi_dy);
        
   
   // phase space function & Poisson bracket
   for(int v=NvLlD; v<=NvLuD;v++) { 
     

      fft->kYIn = 0.;
	  // ky * f1(x,ky,z,..) -> df1_dy(x,y,...)
      //#pragma omp parallel for
      for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) fft->kYIn(RxLD, y_k, RzLD, Field::phi) = fs(RxLD, y_k, RzLD, v, m, s);
      
        fft->solve(FFT_Y, FFT_BACKWARD, NxLD * NzLD);
        xy_f1(RxLD,RyLD,RzLD) = (fft->rYOut(RxLD, RyLD, RzLD, Field::phi));
        setBoundary3(xy_f1);

      
      ///////////////////////   calculate cross terms using Morinishi scheme    [ phi, F1]  //////////////////////////////////
      fft->rYIn = 0.;


      for(int z=NzLlD; z<=NzLuD;z++) { for(int y=NyLlD; y<=NyLuD;y++) { for(int x= NxLlD; x <= NxLuD; x++) {
                // use Morinishi scheme
// complex transformation but only real values are required ... why ?
                const cmplxd dXi_dy__dG_dx = (  8.*(real(xy_dphi_dy(x,y,z)+xy_dphi_dy(x+1,y,z))*real(xy_f1(x+1, y, z)) - real(xy_dphi_dy(x,y,z)+xy_dphi_dy(x-1,y,z))*real(xy_f1(x-1, y, z)))      
                                            - (     real(xy_dphi_dy(x,y,z)+xy_dphi_dy(x+2,y,z))*real(xy_f1(x+2, y, z)) - real(xy_dphi_dy(x,y,z)+xy_dphi_dy(x-2,y,z))*real(xy_f1(x-2, y, z))))/(24.*dx)     ;
                const cmplxd dXi_dx__dG_dy =  ( 8.*(real(xy_dphi_dx(x,y,z)+xy_dphi_dx(x,y+1,z))*real(xy_f1(x, y+1, z)) - real(xy_dphi_dx(x,y,z)+xy_dphi_dx(x,y-1,z))*real(xy_f1(x, y-1, z)))      
                                            - (     real(xy_dphi_dx(x,y,z)+xy_dphi_dx(x,y+2,z))*real(xy_f1(x, y+2, z)) - real(xy_dphi_dx(x,y,z)+xy_dphi_dx(x,y-2,z))*real(xy_f1(x, y-2, z))))/(24.*dy)     ;
        
            // we have always cross-terms thus FFT normalization is easy
	        fft->rYIn(x,y,z, Field::phi) =  (- dXi_dx__dG_dy - dXi_dy__dG_dx)/fft->Norm_Y;
	    
      }}}
     
      fft->solve(FFT_Y, FFT_FORWARD, NxLD * NzLD);

      for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) nonLinearTerms(RxLD,y_k,RzLD,v) = fft->kYOut(RxLD, y_k, RzLD, Field::phi);
    
   }
          
   return nonLinearTerms;

}




Array4z VlasovBlitz::calculatePhiNonLinearityAA(Array5z phi, Array6z fs, const int m, const int s) 
{
     
   // back transform dphi_dx(x,ky,z,m,s) -> dphi_dx(x,y,z,m,s) 
   fft->kYIn = 0.;
   for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) {
   for(int x=NxLlD; x<= NxLuD;x++)  fft->kYIn(x, mapAAY(y_k), RzLD, Field::phi) = (8.*(phi(x+1, y_k,RzLD, m,s) - phi(x-1, y_k, RzLD,m,s)) - (phi(x+2,   y_k, RzLD,m,s) - phi(x-2, y_k  , RzLD,m,s)))/(12.*dx); 
   }
   copyNyquiest(fft->kYIn, fft);
   fft->solve(FFT_Y, FFT_BACKWARD, NxLD * NzLD);
   xy_dphi_dx(RxLD,fft->Y_RyLD,RzLD) = fft->rYOut(RxLD, fft->Y_RyLD, RzLD, Field::phi);
            
   // dphi_dy
   fft->kYIn = 0.;
   //#-pragma omp parallel for
   for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) fft->kYIn(RxLD, mapAAY(y_k), RzLD, Field::phi) = phi(RxLD, y_k, RzLD, m, s) * cmplxd(0., -fft->ky(y_k));
   fft->solve(FFT_Y, FFT_BACKWARD, NxLD * NzLD);
   copyNyquiest(fft->kYIn, fft);
   xy_dphi_dy(RxLD,fft->Y_RyLD,RzLD) = fft->rYOut(RxLD, fft->Y_RyLD, RzLD, Field::phi);
            
   // phase space function & Poisson bracket
   for(int v=NvLlD; v<=NvLuD;v++) { 
       
      fft->kYIn = 0.;
      // FFT-back transformation for (ky->y)  : f(x,ky,z,v,m,s) -> f(x,y,z,v,m,s)
              
   for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) {
      for(int x=NxLlD; x<= NxLuD;x++)   fft->kYIn(x, mapAAY(y_k), RzLD, Field::phi) = (8.*(fs(x+1, y_k  , RzLD, v,m,s) - fs(x-1  , y_k, RzLD, v,m,s)) - (fs(x+2,   y_k, RzLD, v,m,s) - fs(x-2, y_k  , RzLD, v,m,s)))/(12.*dx);  
   }
      copyNyquiest(fft->kYIn, fft);
      fft->solve(FFT_Y, FFT_BACKWARD, NxLD * NzLD);
      xy_df1_dx(RxLD,fft->Y_RyLD,RzLD) = fft->rYOut(RxLD, fft->Y_RyLD, RzLD, Field::phi);
            
      fft->kYIn = 0.;
	  // ky * f1(x,ky,z,..) -> df1_dy(x,y,...)
      //#pragma omp parallel for
      for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) fft->kYIn(RxLD, mapAAY(y_k), RzLD, Field::phi) = fs(RxLD, y_k, RzLD, v, m, s) * cmplxd(0., -fft->ky(y_k));
      copyNyquiest(fft->kYIn, fft);
      fft->solve(FFT_Y, FFT_BACKWARD, NxLD * NzLD);
      xy_df1_dy(RxLD,fft->Y_RyLD,RzLD) = fft->rYOut(RxLD, fft->Y_RyLD, RzLD, Field::phi);

      
      ///////////////////////   calculate cross terms using Morinishi scheme    [ phi, F1]  //////////////////////////////////
      fft->rYIn = 0.;

      for(int z=NzLlD; z<=NzLuD;z++) { for(int y=fft->Y_NyLlD; y<=fft->Y_NyLuD;y++) { for(int x= NxLlD; x <= NxLuD; x++) {
	        const cmplxd dXi_dy__dG_dx = real(xy_dphi_dy(x,y,z)) * real(xy_df1_dx(x,y,z));  
		    const cmplxd dXi_dx__dG_dy = real(xy_dphi_dx(x,y,z)) * real(xy_df1_dy(x,y,z));  
                
            // we have always cross-terms thus FFT normalization is easy
	        fft->rYIn(x,y,z, Field::phi) = (dXi_dx__dG_dy - dXi_dy__dG_dx)/((fft->Norm_Y*3)/2);
	    
      }}}
     
      fft->solve(FFT_Y, FFT_FORWARD, NxLD * NzLD);

      for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) nonLinearTerms(RxLD,y_k,RzLD,v) = fft->kYOut(RxLD, mapAAY(y_k), RzLD, Field::phi);
    
   }
          
   return nonLinearTerms;

}



