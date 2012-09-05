#include "config.h"
#include "Vlasov_Blitz.h"


VlasovBlitz::VlasovBlitz(Grid *_grid, Parallel *_parallel, Setup *_setup, FileIO *fileIO, Geometry *_geo, FFTSolver *fft)    : Vlasov_LenardBernstein(_grid, _parallel, _setup,fileIO,  _geo, fft),
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

  return GKC_SUCCESS;
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



void VlasovBlitz::Vlasov_2D_Island(Array6z fs, Array6z fss, Array6z ft, Fields *fields, double dt, int rk_step) {

  const double w     = setup->get("Island.Width", 10.); 
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

        /////////////////////////////////// Non-Linear Terms ////////////////////////////////////////




       if(nonLinear) {
        nonLinearTerms = (useAntiAliasing) ? calculatePhiNonLinearityAA(fields->phi, fs, m, s): calculatePhiNonLinearity(fields->phi, fs, m, s);

       }


       for(int z=NzLlD; z<= NzLuD;z++) { 
        #pragma omp parallel for
        for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int x=NxLlD; x<= NxLuD;x++) { 
        

             dphi_dx(x,y_k,z)  = (8.*(fields->phi(x+1, y_k  , z, m, s) - fields->phi(x-1  , y_k, z, m, s)) - (fields->phi(x+2,   y_k, z, m, s) - fields->phi(x-2, y_k  , z, m, s)))/(12.*dx)  ;  
             const cmplxd d_ky = cmplxd(0.,-fft->ky(y_k));

             Xi_max(DIR_X) = max(Xi_max(DIR_X), abs(dphi_dx(x,y_k,z)));
             //Xi_max(DIR_Y) = max(Xi_max(DIR_Y), max(abs(d_ky*fields->phi(x,y_k,z,m,s)),abs(alpha*Lv*geo->getShear(x))); 
             //Xi_max(DIR_Z) = 0.;
             Xi_max(DIR_Y) = max(Xi_max(DIR_Y), abs(d_ky*fields->phi(x,y_k,z,m,s))); 
             Xi_max(DIR_Z) = max(Xi_max(DIR_Z), 0.5 * abs(d_ky * geo->getShear(x))); 
             



        }}}

      for(int v=NvLlD; v<= NvLuD;v++) {

      for(int z=NzLlD; z<= NzLuD;z++) { 
        #pragma omp parallel for
        for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int x=NxLlD; x<= NxLuD;x++) {  
        
                      
        const cmplxd d_ky  = cmplxd(0.,-fft->ky(y_k));

        const cmplxd g   = fs(x, y_k, z, v, m, s);
        const cmplxd F0  = f0(x, y_k, z, v, m, s);
        const cmplxd phi = fields->phi(x,y_k, z, m, s); 

        // Hyper diffusion terms (fourth order deriviative)
        const cmplxd d4g_dv    =  (-39. *(fs(x, y_k, z, v+1, m, s) - fs(x, y_k, z, v-1, m, s))  + 12. *(fs(x, y_k, z, v+2, m, s) - fs(x, y_k, z, v-2, m, s)) + 56. * fs(x,y_k,z, v, m, s));///pow4(dv);



        /////////////////////////////////////////////////// Magnetic Island Contribution    /////////////////////////////////////////
                   
        // using ky = ky' +- ky'  and  using ky = -ky' +- ky''
        // what about gyro-average ? this is in gyro-coordinates but island is in particle coordinates
        const double zeta=0.1;
        const double MagIs     = - 0.5 * w*w * shear / 16. * cos(zeta * X(x));

        const double dMagIs_dx = 0.5 * w*w*shear/16. * sin(zeta * X(x)) * zeta;

        // F1 
        // const cmplxd dfs_dx_p1 = ((y_k+1) <= NkyLuD) ? (8. *(fs(x+1, y_k+1, z, v, m, s) - f1(x-1, y_k+1, z, v, m, s))  - (f1(x+2, y_k+1, z, v, m, s) - f1(x-2, y_k+1, z, v, m, s)))/(12.*dx): 0.;
        //  const cmplxd dfs_dx_m1 = ((y_k-1) >= NkyLlD) ? (8. *(fs(x+1, y_k-1, z, v, m, s) - f1(x-1, y_k-1, z, v, m, s))  - (f1(x+2, y_k-1, z, v, m, s) - f1(x-2, y_k-1, z, v, m, s)))/(12.*dx): 0.;

        //const cmplxd ky_fs_m1  = ((y_k-1) >= NkyLlD) ?  cmplxd(0., -fft->ky(y_k-1)) * fs(x,y_k-1,z,v,m,s) : 0.;
        //const cmplxd ky_fs_p1  = ((y_k+1) <= NkyLuD) ?  cmplxd(0., -fft->ky(y_k+1)) * fs(x,y_k+1,z,v,m,s) : 0.;
        
        //const cmplxd ky_fs_m1  = ((y_k-1) >= NkyLlD) ?  cmplxd(0.,  -fft->ky(y_k)) * fs(x,y_k-1,z,v,m,s) : 0.;
        //const cmplxd ky_fs_p1  = ((y_k+1) <= NkyLuD) ?  cmplxd(0.,  -fft->ky(y_k)) * fs(x,y_k+1,z,v,m,s) : 0.;
       

        //const cmplxd ky_fs_m1  = cmplxd(0., -fft->ky(y_k)) * ( ((y_k-1) >= 0   ) ? (((y_k-1) != Ny/2  ) ? fs(x,y_k-1,z,v,m,s) : 0.): fs(x,Ny-1,z,v,m,s) );


        // how does they couple with Nyquist frequencty ??
        // take care of boundaries !!
       cmplxd dfs_dx_p1         =  ((y_k+1) <= Ny-1) ?
                              	   (8. *(fs(x+1, y_k+1, z, v, m, s) - fs(x-1, y_k+1, z, v, m, s))  - (fs(x+2, y_k+1, z, v, m, s) - fs(x-2, y_k+1, z, v, m, s)))/(12.*dx) :
                                   (8. *(fs(x+1,     0, z, v, m, s) - fs(x-1,     0, z, v, m, s))  - (fs(x+2,     0, z, v, m, s) - fs(x-2,     0, z, v, m, s)))/(12.*dx);

        const cmplxd dfs_dx_m1  =  ((y_k-1) >= 0   ) ? 
                                   (8. *(fs(x+1, y_k-1, z, v, m, s) - fs(x-1, y_k-1, z, v, m, s))  - (fs(x+2, y_k-1, z, v, m, s) - fs(x-2, y_k-1, z, v, m, s)))/(12.*dx) :
                                   (8. *(fs(x+1, Ny -1, z, v, m, s) - fs(x-1, Ny -1, z, v, m, s))  - (fs(x+2, Ny -1, z, v, m, s) - fs(x-2, Ny -1, z, v, m, s)))/(12.*dx);

        cmplxd ky_fs_p1         = cmplxd(0., -fft->ky(y_k)) * ( ((y_k+1) <= Ny-1) ? fs(x,y_k+1,z,v,m,s) : fs(x,   0,z,v,m,s) );
        const cmplxd ky_fs_m1   = cmplxd(0., -fft->ky(y_k)) * ( ((y_k-1) >= 0   ) ? fs(x,y_k-1,z,v,m,s) : fs(x,Ny-1,z,v,m,s) );





        // not at the Nyquist frequency we have no coupling with higher frequencies
	    if(y_k == Ny/2) ky_fs_p1 =  0.;
	    if(y_k == Ny/2) dfs_dx_p1 =  0.;


        // phi
        //const cmplxd ky_phi_p1 =  ((y_k+1) <= NkyLuD) ? cmplxd(0., -fft->ky(y_k+1)) * fields->phi(x,y_k+1,z,m,s) : 0.;

              cmplxd dphi_dx_p1 = ((y_k+1) <= Ny-1) ? dphi_dx(x,y_k+1, z): dphi_dx(x,  0, z);
        const cmplxd dphi_dx_m1 = ((y_k-1) >=    0) ? dphi_dx(x,y_k-1, z): dphi_dx(x,Ny-1,z);
              cmplxd ky_phi_p1 =  cmplxd(0., -fft->ky(y_k)) * (((y_k+1) <= Ny-1) ? fields->phi(x,y_k+1,z,m,s) : fields->phi(x,    0,z,m,s));
        const cmplxd ky_phi_m1 =  cmplxd(0., -fft->ky(y_k)) * (((y_k-1) >=    0) ? fields->phi(x,y_k-1,z,m,s) : fields->phi(x, Ny-1,z,m,s));
        
        // not at the Nyquist frequency we have no coupling with higher frequencies
	    if(y_k == Ny/2) dphi_dx_p1 =  0.;
	    if(y_k == Ny/2)  ky_phi_p1 =  0.;
          
               const cmplxd dky_1 = cmplxd(0., -fft->ky(1));//1.;//cmplxd(0.,-fft->ky(1.));
	
	// mode-mode connections
        const cmplxd Island =
                           //                                     -  alpha * V(v) * (   -  cmplxd(0., 1.) * MagIs  * (dfs_dx_m1  - dfs_dx_p1 )) ; 
                             - alpha * V(v) * ( dMagIs_dx * (ky_fs_m1  + ky_fs_p1 )  -  dky_1 * MagIs  * (dfs_dx_p1  - dfs_dx_m1 ))  
                             + alpha * V(v) * ( dMagIs_dx * (ky_phi_m1 + ky_phi_p1) -   dky_1 * MagIs *  (dphi_dx_p1 + dphi_dx_m1)) * F0;

             // (((y_k-1) < NkyLlG) ? 0. : cmplxd(0.,-fft->ky(y_k-1)) * fs(x,y_k-1,z,v,m,s)) + (((y_k+1) > NkyLuG) ?  0. : cmplxd(0.,-fft->ky(y_k+1)) * fs(x,-y_k+1,z,v,m,s)))); 
//            -   Island_m1 * Ax *  dfs_dx));
             // [A_\parallel, F_1] + [A_\paralell, \phi]
 /*      (
              (dA_dx * (
              (((y_k-1) < NkyLlD) ? 0. : dky_p_m1 * fs(x,y_k-1,z,v,m,s)) - (((-y_k+1) > NkyLuD) ?  0. : dky_m_p1 * fs(x,-y_k+1,z,v,m,s)) 
            + (((y_k+1) > NkyLuD) ? 0. : dky_p_p1 * fs(x,y_k+1,z,v,m,s)) - (((-y_k-1) < NkyLlD) ?  0. : dky_m_m1 * fs(x,-y_k-1,z,v,m,s)) )
            -   Island_m1 * zeta * Ax *  dfs_dx)
            + (dA_dx * (
              (((y_k-1) < NkyLlD) ? 0. : dky_p_m1 * fields->phi(x,y_k-1,z,m,s)) - (((-y_k+1) > NkyLuD) ?  0. : dky_m_p1 * fields->phi(x,-y_k+1,z,m,s)) 
            + (((y_k+1) > NkyLuD) ? 0. : dky_p_p1 * fields->phi(x,y_k+1,z,m,s)) - (((-y_k-1) < NkyLlD) ?  0. : dky_m_m1 * fields->phi(x,-y_k-1,z,m,s)) )
            -   Island_m1 * zeta * Ax *  dphi_dx(x,y_k,z)) );
*/
        /////////////// Finally the Vlasov equation calculate the time derivatve      //////////////////////
 
        cmplxd dg_dt = 
          //- Nonlinearity (in case of Linear simulation term is zero  -  dXi_dx__dG_dy - dXi_dy__dG_dx
        //nonLinearTerms(x,y_k,z,v) +
         Island + 
          // include electro-magnetic contributions
  //         + (((plasma->nfields >= 2) && !nonLinear) ? alpha * plasma->beta  * V(v) * (dA_dx__dfs_dy - dA_dy__dfs_dx) : 0.)
             
             // driving term (use dphi_dy instead of dXi_dy, because v * A does not vanish due to numerical errors)
             d_ky * (-(w_n + w_T * ((pow2(V(v))+ M(m))/Temp  - sub)) * F0 * phi

             // add first order gyro-average term (zero when full-gyro)
             + 0.5 * w_T  * k2p_phi(x,y_k,z) * F0

      	     // Landau Damping term and parallel ... ? 
             - alpha  * V(v)* geo->getShear(x)  * ( g + sigma * phi * F0)) - 1.e-6 * d4g_dv 
        // Collisional Term
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

/*  
void VlasovBlitz::Vlasov_2D_Island(Array6z fs, Array6z fss, Array6z ft, Fields *fields, double dt, int rk_step) {

  const double w     = setup->get("Island.Width", 0.); 
  const double shear = setup->get("Geometry.Shear", 0.); 
  const bool useDorland1 = (fields->gyroAverageModel == "Dorland1");
   
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
	
       // back transform f(x,ky,z,v,m,s) -> f(x,y,z,v)
       for(int z=NzLlD; z<= NzLuD;z++) {  for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int v=NvLlD; v<= NvLuD;v++) { 

            fft->k1InY(fft->Rk1yL) = fs(x, RkyLD, z, v, m, s);
            fft->solve(FFT_1DY, FFT_BACKWARD, Field::phi);
            f(x,RyLD,v,m,s) = fft->r1OutY(RyLD);
    
        }}}

       // calculate non-linear terms

        if(nonLinear) {
 
            // we use aliasing, that so our input array is 3/2 wider
            phi_xy = fft->solve(FFT_1D, phi, FFT::ALIASING);
            phi_xy = fft->solve(FFT_1D, phi, FFT::ALIASING);




/// To calculate non-linear terms we use Mornishi-stencil
    //       dXi_dy__dG_dx =  ( 8.*((dphi_dy(x,y,z)+dphi_dy(x+1,y,z))*fs(x+1, y, z, v, m, s) - (dphi_dy(x,y,z)+dphi_dy(x-1,y,z))*fs(x-1, y, z, v, m, s))      
    //                                - ((dphi_dy(x,y,z)+dphi_dy(x+2,y,z))*fs(x+2, y, z, v, m, s) - (dphi_dy(x,y,z)+dphi_dy(x-2,y,z))*fs(x-2, y, z, v, m, s)))/(24.*dx)     ;

     //      dXi_dx__dG_dy =  ( 8.*((dphi_dx(x,y,z)+dphi_dx(x,y+1,z))*fs(x, y+1, z, v, m, s) - (dphi_dx(x,y,z)+dphi_dx(x,y-1,z))*fs(x, y-1, z, v, m, s))      
     //                                  - ((dphi_dx(x,y,z)+dphi_dx(x,y+2,z))*fs(x, y+2, z, v, m, s) - (dphi_dx(x,y,z)+dphi_dx(x,y-2,z))*fs(x, y-2, z, v, m, s)))/(24.*dy)     ;
        }






	   // take derivatives for dphi_dx, dphi_dy
       for(int z=NzLlD; z<= NzLuD;z++) { for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int x=NxLlB; x<= NxLuB;x++) { 
        
             dphi_dx(x,y_k,z)  = (8.*(fields->phi(x+1, y_k  , z, m, s) - fields->phi(x-1  , y_k, z, m, s)) - (fields->phi(x+2,   y_k, z, m, s) - fields->phi(x-2, y_k  , z, m, s)))/(12.*dx)  ;  

             const cmplxd d_ky = cmplxd(0.,-fft->ky(y_k));

             Xi_max(DIR_X) = max(Xi_max(DIR_X), abs(dphi_dx(x,y_k,z)));
             Xi_max(DIR_Y) = max(Xi_max(DIR_Y), abs(d_ky*fields->phi(x,y_k,z))); 
             Xi_max(DIR_Z) = 0.;


        }}}
      
      for(int v=NvLlD; v<= NvLuD;v++) { 
      for(int z=NzLlD; z<= NzLuD;z++) { for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) {  for(int x=NxLlD; x<= NxLuD;x++) { 
                      
        const cmplxd d_ky  = cmplxd(0.,fft->ky(y_k));

        const cmplxd g   = fs(x, y_k, z, v, m, s);
        const cmplxd F0  = f0(x, y_k, z, v, m, s);
        const cmplxd phi = fields->phi(x,y_k, z, m, s); 

        
        double dXi_dy__dG_dx = 0.;
        double dXi_dx__dG_dy = 0.;


        /////////////////////////////////// Non-Linear Terms ////////////////////////////////////////




        // collisions
//        const double df1_dv    = (8. *(f1(x, y, z, v+1, m, s) - f1(x, y, z, v-1, m, s))  -1. *(f1(x, y, z, v+2, m, s) - f1(x, y, z, v-2, m, s)))/(12.*dv);
//        const cmplxd ddfs_dvv  = (16. *(fs(x, y_k, z, v+1, m, s) + fs(x, y_k, z, v-1, m, s))  -1. *(fs(x, y_k, z, v+2, m, s) + fs(x, y_k, z, v-2, m, s)) - 30.*fs(x,y_k,z,v,m,s))/(12.*pow2(dv));
//        const double beta=2.e-3;
//        const cmplxd v2_rms = 1.;//pow2(plasma->species(s).scale_v);
// finite difference stencil [ -8 1 0 -1 8 ]

        // magnetic island
        const double zeta=0.1;
        const double Ax = - w*w * shear / 16. * cos(zeta * x);

        const double dA_dx = w*w*shear/16. * sin(zeta * x) * zeta;
        const cmplxd dA_dy = Ax;
        const cmplxd dfs_dx  = (8. *(fs(x+1, y_k, z, v, m, s) - f1(x-1, y_k, z, v, m, s))  - (f1(x+2, y_k, z, v, m, s) - f1(x-2, y_k, z, v, m, s)))/(12.*dx);

        /////////////// Finally the Vlasov equation calculate the time derivatve      //////////////////////
 
        cmplxd dg_dt = 
          //- Nonlinearity (in case of Linear simulation term is zero
//          -  dXi_dx__dG_dy - dXi_dy__dG_dx
          // include electro-magnetic contributions
  //         + (((plasma->nfields >= 2) && !nonLinear) ? alpha * plasma->beta  * V(v) * (dA_dx__dfs_dy - dA_dy__dfs_dx) : 0.)
             
             // driving term (use dphi_dy instead of dXi_dy, because v * A does not vanish due to numerical errors)
             d_ky * (-(w_n + w_T * ((pow2(V(v))+ M(m))/Temp  - sub)) * F0 * phi

             // add first order gyro-average term (zero when full-gyro)
             + 0.5 * w_T  * k2p_phi(x,y_k,z) * F0

      	     // Landau Damping term and parallel ... ? 
             - alpha  * V(v)* geo->getShear(x)  * ( g + sigma * phi * F0))

            //   Magnetic Island
//            	Mode coupling table
  
//	                A_1\parallel  = ..
//
                   //dA_dx * dfs_dyi					 + dA_dy * dfs_dx
//  
// 
//                    const cmplxd ky = -cmplxd(0.,1.) * fft->ky(y_k);
//
                   // using ky = ky' +- ky'  and  using ky = -ky' +- ky''
                   // what about gyro-average ? this is in gyro-coordinates but island is in particle coordinates
             + alpha * V(v) * (dA_dx * (
              (((y_k-1) < NkyLlD) ? 0. : cmplxd(0., -fft->ky(y_k-1)) * fs(x,y_k-1,z,v,m,s)) - (((-y_k+1) > NkyLuD) ?  0. : cmplxd(0.,-fft->ky(-y_k+1)) * fs(x,-y_k+1,z,v,m,s)) 
            +     (((y_k+1) > NkyLuD) ? 0. : cmplxd(0., -fft->ky(y_k+1)) * fs(x,y_k+1,z,v,m,s)) - (((-y_k-1) < NkyLlD) ?  0. : cmplxd(0.,-fft->ky(-y_k-1)) * fs(x,-y_k-1,z,v,m,s)) )
            -   (cmplxd(0.,-fft->ky(1)) * zeta * dA_dy *  dfs_dx))
;           + alpha * V(v) * (dA_dx * (
             (((y_k-1) < NkyLlD) ? 0. : cmplxd(0., -fft->ky(y_k-1)) * fields->phi(x,y_k-1,z,m,s)) - (((-y_k+1) > NkyLuD) ?  0. : cmplxd(0.,-fft->ky(-y_k+1)) * fields->phi(x,-y_k+1,z,m,s)) 
            +     (((y_k+1) > NkyLuD) ? 0. : cmplxd(0., -fft->ky(y_k+1)) * fields->phi(x,y_k+1,z,m,s)) - (((-y_k-1) < NkyLlD) ?  0. : cmplxd(0.,-fft->ky(-y_k-1)) * fields->phi(x,-y_k-1,z,m,s)) )
            -   (cmplxd(0.,-fft->ky(1)) * zeta * dA_dy *  dphi_dx))
     ;
        // Collisional Term
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



  */

  /*  
  void VlasovBlitz::Vlasov_EM_Linear   (Array6d fs, Array6d fss, Array6d ft, Fields *fields, double dt, int rk_step) {
       
  
  const double B0 = plasma->B0;

       for(int s = NsLlD; s <= NsLuD; s++) {
        
         // small abbrevations
         const double w_n = plasma->species(s).w_n;
         const double w_T = plasma->species(s).w_T;

        const double a = plasma->species(s).alpha;
        const double sigma = plasma->species(s).sigma;

        const double m = plasma->species(s).m;
        const double T = plasma->species(s).T0;
        const double q = plasma->species(s).q;
    
        const double B = plasma->B0;

        // what is eps ?
        const double eps  = 1.;
        const double sub = (plasma->species(s).doGyro) ? 3./2. : 1./2.;

          for(int m=NmLlB; m<= NmLuB;m++)  { for(int v=NvLlB; v<= NvLuB;v++) { for(int z=NzLlB; z<= NzLuB;z++){
          
            for(int y=NyLlB-2; y<= NyLuB+2;y++) { for(int x=NxLlB-2; x<= NxLuB+2;x++) { 
                            Xi(x,y,z,v,m,s) = fields->phi(x,y,z,m,s);//  - a* eps * plasma->beta * V(v)*fields->Ap(x,y,z,m,s);
            }}

          for(int y=NyLlB; y<= NyLuB;y++) { for(int x=NxLlB; x<= NxLuB+2;x++) { 

                            G(x,y,z,v,m,s)  = fs(x,y,z,v,m,s) - sigma * Xi(x,y,z,v,m,s) * f0(x,y,z,v,m,s);
          }}
        
          }}} 


          for(int m=NmLlD; m<= NmLuD;m++)  {
        for(int v=NvLlD; v<= NvLuD;v++) { 
          for(int z=NzLlD; z<= NzLuD;z++){
          for(int y=NyLlD; y<= NyLuD;y++) { 
        for(int x=NxLlD; x<= NxLuD;x++) { 
                      

                        // We use CD-4 (central difference fourth order for every variable)
                       const double dXi_dz     = (8. *(Xi(x, y  , z+1, v, m, s) - Xi(x  , y  , z-1, v, m, s))  -1. *(Xi(x  , y  , z+2, v, m, s) - Xi(x  , y  , z-2, v, m, s)))/(12.*dz); 
                       const double dXi_dy     = (8. *(Xi(x, y+1, z  , v, m, s) - Xi(x  , y-1, z  , v, m, s))  -1. *(Xi(x  , y+2, z  , v, m, s) - Xi(x  , y-2, z  , v, m, s)))/(12.*dy)  ;  
                      
                       const double dfs_dv   = (8. *(fs(x, y, z, v+1, m, s) - fs(x, y, z, v-1, m, s))  -1. *(fs(x, y, z, v+2, m, s) - fs(x, y, z, v-2, m, s)))/(12.*dv);

                       const double dG_dx   = (8.*(G(x, y, z+1, v, m, s) - G(x, y, z-1, v, m, s))    -1.*(G(x, y, z+2, v, m, s) - G(x, y, z-2, v, m, s)))/(12.*dx);
                       const double dG_dy   = (8.*(G(x, y, z+1, v, m, s) - G(x, y, z-1, v, m, s))    -1.*(G(x, y, z+2, v, m, s) - G(x, y, z-2, v, m, s)))/(12.*dy);
                       const double dG_dz   = (8.*(G(x, y, z+1, v, m, s) - G(x, y, z-1, v, m, s))    -1.*(G(x, y, z+2, v, m, s) - G(x, y, z-2, v, m, s)))/(12.*dz);

                       const double F0      =  f0(x,y,z,v,m,s);

                        // magnetic prefactor defined as  $ \hat{B}_0 / \hat{B}_{0\parallel}^\star = \left[ 1 + \beta_{ref} \sqrt{\frac{\hat{m_\sigma T_{0\sigma}{2}}}}
                        //                                                                           \frac{\hat{j}_{0\parallel}}{\hat{q}_sigma \hat{B}_0^2 v_\parallel  \right]^{-1}
                        //
                        // note j0 is calculated and needs to be replaced, or ? no we calculate j1 ne ?!
                       const double j0 = 0.;
                       const double Bpre  = 1.; //1./(1. + plasma->beta * sqrt(m * T/2.) * j0 / (q * pow2(geo->B(x,y,z))) * V(v));
                   
                       const double CoJB = 1./geo->J(x,y,z);

                        // Finally the Vlasov equation calculate the time derivatve      
                        double fts = 
                        // driving term
                             + Bpre * (w_n + w_T * ((pow2(V(v))+ M(m) * B0)/T - sub)) * F0 * dXi_dy
                             - Bpre * sigma * ((M(m) * B0 + 2.*pow2(V(v)))/B0 * geo->Kx(x,y,z)) * dG_dx
                             - Bpre * sigma * ((M(m) * B0 + 2.*pow2(V(v)))/B0 * geo->Ky(x,y,z) - a * pow2(V(v)) * plasma->beta * plasma->w_p) * dG_dy
                             -  CoJB * ( a * V(v)* dG_dz)
                             + a * V(v) / 2. * M(m) * geo->dB_dz(x,y,z) * dfs_dv
                             + Bpre *  sigma * (M(m) * B0 + 2. * pow2(V(v)))/B0 * geo->Kx(x,y,z) * (w_n + w_T * (pow2(V(v)) + M(m) * B0)/T - 3./2.) * F0;
                             // - Nonlinearity - (dXi_dx * dG_dy - dXi_dy * dG_dx)
                             ;
// add hyperdiffusive terms, collisions, heating 
                        
                             
                         //    const double fts = alpha * Vel(v)* dH_dz
                         //    - dXi_dy * (plasma->species(s).w_n + plasma->species(s).w_T * (Vel(v)**2.e0+Mu(m)* B0  - sub))  * f0(x, y, z, v, m, s ) 

        //  time-integrate the distribution function     
        if(rk_step == 1) ft(x, y, z, v, m, s) = fts;
        else if((rk_step == 2) || (rk_step == 3)) ft(x, y, z, v, m, s) = ft(x, y, z, v, m, s) + 2.0*fts;
        else    fts = ft(x, y, z, v, m, s) + fts;
        
        fss(x, y, z, v, m, s) = f(x, y, z, v, m, s) + fts*dt;


        }}} }}}


       
       }


void VlasovBlitz::printOn(ostream &output) const
        {

            output << "Vlasov     |   Spatial : Morinishi        Time : Runge-Kutta 4 " << std::endl;
            output << "Vlasov     |   Hyper Viscosity : " << hyper_visc << std::endl;
            output << " Vlasov    | " << (calcTerms & LINEAR ? "Linear," : "") << 
                                        (calcTerms & NON_LINEAR ? "NonLinear," : "") << 
                                        (calcTerms & GEO_TERM ? "Geo," : "") << 
                                        (calcTerms & PTRAP_TERM ? "PTrap," : "") << 
                                        (calcTerms & VPARALLEL_TERM ? "VParallel," : "") << std::endl; 


        };




*/




