/*
 * =====================================================================================
 *
 *       Filename:  Init.cpp
 *
 *    Description: Definition of Initial conditions 
 *
 *         Author:  Paul P. Hilscher (2009-), 
 *
 * =====================================================================================
 */

#include "Init.h"
#include "Plasma.h"
#include "Special/HermitePoly.h"

#include "Reader/ReaderXYV.h"
#include "System.h"

Init::Init(Parallel *parallel, Grid *grid, Setup *setup, Vlasov *vlasov, Fields *fields, Geometry *_geo) : geo(_geo) {

   epsilon_0          = setup->get("Init.Epsilon0", 1.e-14); 
   sigma              = setup->get("Init.Sigma"   , 3.e-1); 


   ////////////////////////////////////////////////////  Initial Condition Maxwellian f0 = (...) ///////////////
   for(int s = NsLlD; s <= NsLuD; s++) {
   
     // Initialize Form of f, Ap, phi, and g, we need superposition between genereal f1 pertubration and species dependent
   std::string perturb_f1s_str = setup->get("Plasma.Species" + Setup::number2string(s) + ".InitF1", "0.");
   std::string perturb_f1_str  = setup->get("Init.F1", "0.");
   FunctionParser f1s_parser = setup->getFParser();
   check(((f1s_parser.Parse(perturb_f1s_str + "+" +  perturb_f1_str, "x,y,z,v,m") == -1) ? 1 : -1), DMESG("Parsing error of Initial condition n(x)"));
                                 
   	const double VOff = setup->get("Plasma.Species" + Setup::number2string(s) + ".VelocityOffset", 0.);
   
   	for(int m = NmLlD; m <= NmLuD; m++) { for(int v = NvLlD; v <= NvLuD; v++) { 
   	for(int z = NzLlD; z <= NzLuD; z++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) {  for(int x = NxLlD; x <= NxLuD; x++) {
       
                    const double n = plasma->species(s).n(x);
                    const double T = plasma->species(s).T(x);

                    // included rotation, not additionally the flux surfacedensity is also scaled withing n = n_0 exp(e/t) 
                    const double w = 0., r = 0;

		            // although only F0(x,k_y=0,...) is not equal zero, we perturb all modes, as F0 in Fourier space "acts" like a nonlinearity, which couples modes together
                    vlasov->f0(x, y_k, z, v, m, s)  =  (n / pow( M_PI*T, 1.5) * exp(-pow2(V(v) - w*r + VOff)/T) * ((plasma->species(s).doGyro == true) ?   exp(- M(m)    * plasma->B0/T) :  T/(plasma->B0)));

                   // for df
                   double pos[] = {X(x), 0., Z(z), V(v), M(m) };
                   vlasov->f (x,y_k,z,v,m,s) = f1s_parser.Eval(pos)*vlasov->f0(x,y_k,z,v,m,s);

   }}} }}

   if(plasma->global == false)  vlasov->f(RxLB, RkyLD, RzLB, RvLB, RmLB, RsLD) = 0.e0;
   else                         vlasov->f = vlasov->f0;
   
   // check for predefined perturbations 
   PerturbationMethod = setup->get("Init.Perturbation", "");
   if     (PerturbationMethod == "PSFEqualModePower") PerturbationPSFMode (vlasov, s, plasma->global ? 1. : 0.); 
   else if(PerturbationMethod == "PSFNoise")          PerturbationPSFNoise(vlasov, s, plasma->global ? 1. : 0.);
   else if(PerturbationMethod == "PSFExp")            PerturbationPSFExp  (vlasov, s, plasma->global ? 1. : 0.);
   else if(PerturbationMethod == "PSFHermitePoly")    PerturbationHermitePolynomial(vlasov, s, plasma->global ? 1. : 0., setup->get("Init.HermitePolynomial", 0));
   else if(PerturbationMethod == "NoPerturbation") ;
   else check(-1, DMESG("No such Perturbation Method"));  


   }
   
   ////////////////////////////////////////////// Recheck Charge Neutrality for Perturabtion and correct  /////////////////////////

/*
   if(Ns > 1 && setup->get("Init.ChargeNeutral", 0)) {
     check(-1, DMESG("Buggy"));
     double total_charge = 0.;
   	for(int m = NmLlD; m <= NmLuD; m++) total_charge = parallel->collect((NkyLlD == 0) ? sum(abs(fields->calculateChargeDensity((A6zz) vlasov->f0.dataZero(), (A6zz) vlasov->f.dataZero(), (A4zz) fields->Field0.dataZero(), m, 1)(RxLD, 0, RzLD, Q::rho))) : 0.); 
    
        for(int s = NsLlD; s <= NsLuD; s++) vlasov->f(RxLB, RkyLD, RzLB, RvLB, RmLB, s) *= (1. - total_charge/((double) Ns));
   }
*/

   /////////////////////////////////////// Initialize dynamic Fields for Ap (we use Canonical Momentum Method) ////////////////////////////
   if(plasma->nfields >= 2) {
   	FunctionParser Ap_parser = setup->getFParser();
   	FunctionParser phi_parser = setup->getFParser();
   	check(((phi_parser.Parse(setup->get("Init.phi", "0."), "x,y_k,z") == -1) ? 1 : -1), DMESG("Parsing error of Initial condition n(x)"));
   	check((( Ap_parser.Parse(setup->get("Init.Ap ", "0."), "x,y_k,z") == -1) ? 1 : -1), DMESG("Parsing error of Initial condition n(x)"));

   	for(int s = NsLlD; s <= NsLuB; s++) { for(int m = NmLlD; m <= NmLuD; m++) {
   
     	for(int z = NzLlD; z <= NzLuD; z++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) {  for(int x = NxLlD; x <= NxLuD; x++) {

         	//const Complex pos[3] = { X(x), Y(y_k), Z(z) };
         	const double pos[3] = { X(x), y_k, Z(z) };
         //	fields->Field0 (x,y_k,z,Field::Ap ) = (y_k == 1) ? Complex(sqrt(2.),sqrt(2.)) * Ap_parser.Eval(pos) : 0.;
         //	fields->Field0(x,y_k,z, Field::phi) = Complex(sqrt(2.), sqrt(2.)) * phi_parser.Eval(pos);
  
     	} } }

       	fields->gyroAverage(fields->Field0(RxLD,RkyLD, RzLD, RFields), fields->Field(RxLD, RkyLD, RzLD, m , s, RFields), m, s, true);

    	} }
           
   }






/*
   /////////////////////////////////////// Initialize dynamic Fields for Bp (we use Canonical Momentum Method) ////////////////////////////
   if(plasma->nfields >= 3) {

   	FunctionParser Bp_parser = setup->getFParser();
   	check(((Bp_parser.Parse(setup->get("Init.Bp", "0."), "x,y,z") == -1) ? 1 : -1), DMESG("Parsing error of Initial condition n(x)"));

   	for(int s = NsLlD; s <= NsLuB; s++) { for(int m = NmLlD; m <= NmLuB; m++) {
   
     	for(int z = NzLlD; z <= NzLuB; z++) { for(int y = NkyLlD; y <= NkyLuD; y++) {  for(int x = NxLlD; x <= NxLuD; x++) {

         	const double pos[3] = { X(x), Y(y), Z(z) };
         	fields->Bp(x,y,z,m,s ) = Bp_parser.Eval(pos);
  
     	}}}

        fields->Bp (RxLD, RyLD, RzLD, m, s) = fields->gyroAverage(fields->Bp (RxLD,RyLD, RzLD,  m, s), m, s, Field::phi);
   
   ` 	}}
           
   }
*/

   //////////////////////////////// construct g from f : g = f + .... .//////////////////////////////////////////////////// 
   // we defined g = f + sigma_j * alpha_j * vp * F_j0 eps * berta * Ap
   
   for(int s = NsLlD; s <= NsLuD; s++) { for(int m   = NmLlD ; m   <= NmLuD ; m++  ) {  for(int v = NvLlD; v <= NvLuD; v++) { 
   for(int z = NzLlD; z <= NzLuD; z++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) {  for(int x = NxLlD; x <= NxLuD; x++) {
 
                   vlasov->f(x, y_k, z, v, m, s)  = vlasov->f(x,y_k,z,v,m,s)
                             + ((plasma->nfields >= 2) ? plasma->species(s).sigma * plasma->species(s).alpha * V(v)*geo->eps_hat 
                                                         * vlasov->f0(x,y_k,z,v,m,s) * plasma->beta * fields->Ap(x,y_k,z,m,s) : 0.);
   }}} }}}


   ////////////////////////////////////////////////////////    Set Fixed Fields  phi, Ap, Bp //////////////////
   // set fixed fields if requested,  initialize Fields, Ap

   if (setup->get("Init.FixedPhi", "0.").substr(0,4) == "File") setFieldFromDataFile(setup, fields->Field0, Field::phi, setup->get("Init.FixedPhi", "0."));
   else if (plasma->nfields >= 1) setFieldFromFunction(setup, fields->Field0, Field::phi, setup->get("Init.FixedPhi", "0."));

   if (setup->get("Init.FixedAp" , "0.").substr(0,4) == "File") setFieldFromDataFile(setup, fields->Field0, Field::Ap, setup->get("Init.FixedAp", "0."));
   else if (plasma->nfields >= 2) setFieldFromFunction(setup, fields->Field0, Field::Ap, setup->get("Init.FixedAp", "0."));

   if (setup->get("Init.FixedBp", "0.").substr(0,4) == "File") setFieldFromDataFile(setup, fields->Field0, Field::Bp, setup->get("Init.FixedBp", "0."));
   else if (plasma->nfields >= 3) setFieldFromFunction(setup, fields->Field0, Field::Bp, setup->get("Init.FixedBp", "0."));


  ////////////////////////////////////////////////////////////// 
   // Perform gyro-average of the fields
   for(int s = NsLlD; s <= NsLuD; s++) { for(int m = NmLlD; m <= NmLuD; m++) {
        
       	fields->gyroAverage(fields->Field0, fields->Field(RxLD, RkyLD, RzLD, m , s, RFields), m, s, true);
   
   }}
   
   //////////////////////////////////////// Done ///////////////
   // Boundaries and Done   
   vlasov->setBoundary( vlasov->f0   );
   vlasov->setBoundary( vlasov->f    );
   fields->updateBoundary(); // fields->Field);


};

int Init::PerturbationPSFNoise(Vlasov *vlasov, int s, double pre) {
    
    // add s to initialization of RNG due to fast iteration over s (which time is not resolved)
    ranlib::UniformClosed<double> wNoise;
    wNoise.seed(System::getTime() + System::getProcessID() + s);

   
    for(int s = NsLlD; s <= NsLuD; s++) { for(int m = NmLlD; m <= NmLuD; m++) { for(int v = NvLlD; v<=NvLuD; v++) {
    for(int x = NxLlD; x <= NxLuD; x++) { for(int y = NkyLlD; y <= NkyLuD; y++) { for(int z = NzLlD; z<=NzLuD; z++) {
      vlasov->f(x,y,z,v,m,s) = vlasov->f(x,y,z,v,m,s) + epsilon_0*(wNoise.random()-0.5e0) * vlasov->f0(x,y,z,v,m,s);
    }}}}}}

   return GKC_SUCCESS;
}   


int Init::PerturbationPSFExp(Vlasov *vlasov, int s, double pert) {
   
   for(int s = NsLlD; s <= NsLuD; s++) { for(int m = NmLlD; m <= NmLuD; m++) { for(int v = NvLlD; v<=NvLuD; v++) {
   for(int x = NxLlD; x <= NxLuD; x++) { for(int y = NkyLlD; y <= NkyLuD; y++) { for(int z = NzLlD; z<=NzLuD; z++) {
      vlasov->f(x,y,z,v,m,s) += vlasov->f0(x,y,z,v,m,s) * (pert + Perturbation(x,y,z,1.e-10,0.1));
   }}}}}}
  
   return GKC_SUCCESS;
}


int Init::PerturbationPSFMode(Vlasov *vlasov, int s, double pre) {
 
   //  Callculates the phase (of what ??!)
   auto Phase = [=] (const int q, const int N)  -> double { return 2.*M_PI*((double) (q-1)/N); };
   // check if value is reasonable
   if(!((s >= NsLlD) && (s <= NsLuD))) return GKC_SUCCESS; 
       
   for(int x = NxLlD; x <= NxLuD; x++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { for(int z = NzLlD; z<=NzLuD; z++) {
      
      double pert_x=0., pert_y=0., pert_z=0.;
      
      pert_z = (Nz == 1) ? 1. : 0.;
      //for(int m = 0; m <= Ny/2; m++) pert_y += 1.e-6 * cos(m*(2.*M_PI*Y(y)/Ly   +Phase(m,Ny)));
      for(int r = 1; r <= Nx/2; r++) pert_x += cos(r*(2.*M_PI*X(x)/Lx)+Phase(r,Nx));//*exp(pow2(M_PI*X(x)/Lx));
      for(int m = 1; m <= Nky-1; m++) pert_y += 0.;
      //for(int m = 1; m <= Ny/2; m++) pert_y += exp(pow2(2.*M_PI/Ly*y_k));                                            //     cos(m*(2.*M_PI*Y(y)/Ly)+Phase(m,Ny))*exp(pow2(M_PI*Y(y)/Ly));
      for(int n = 1; n <= Nz/2; n++) pert_z += cos(n*(2.*M_PI*Z(z)/Lz)+Phase(n,Nz));//*exp(pow2(M_PI*Z(z)/Lz));
      
      if(pert_x == 0.) pert_x = 1.;
      if(pert_y == 0.) pert_y = 1.;
      if(pert_z == 0.) pert_z = 1.;
        // for(int r = 1; r <= Nx/2; r++) pert_x +=  cos(r*(2.*M_PI*X(x)/Lx)+Phase(r,Nx));
        // for(int n = 1; n <= Nz/2; n++) pert_z +=  cos(n*(2.*M_PI*Z(z)/Lz)+Phase(n,Nz));
        // vlasov->f(x, y, z, RvLD, RmLD, s)  = (pre + (pert_x*pert_y*pert_z)) * vlasov->f0(x,y,z,RvLD, RmLD, s) / (Nx * Ny * Nz);
        vlasov->f(x, y_k, z, RvLD, RmLD, s)  =  epsilon_0 * (pert_x*pert_y*pert_z) * vlasov->f0(x,y_k,z,RvLD, RmLD, s) * (1./ (Nx * Nky * Nz));

    }}}

   return GKC_SUCCESS;
}



double inline Init::Perturbation(int x,int y,int z, const double epsilon_0, const double sigma) {

            return  epsilon_0*exp(-(  pow2(X(x)/Lx) + pow2(Y(y)/Ly - 0.5) + pow2(Z(z)/Lz - 0.5))/(2.*pow2(sigma)));
}




void Init::setFieldFromDataFile(Setup *setup, Array4C Field0, int n, std::string path) {

  
   std::vector<std::string> token = Setup::split(path, ":");
   std::cout << "Use Reader : " << token[1] << std::endl;
   
   Reader *reader = new ReaderXYV(token[2],setup);

   for(int z = NzLlD; z <= NzLuD; z++) { for(int y = NkyLlD; y <= NkyLuD; y++) {  for(int x = NxLlD; x <= NxLuD; x++) {
        
          Field0(x,y,z,n) = reader->getValue(X(x),Y(y),Z(z));

   }}}

    delete(reader);

}
   

void Init::setFieldFromFunction(Setup *setup, Array4C Field0, int n , std::string func) {

   FunctionParser parser = setup->getFParser();

   check(((parser.Parse (func, "x,y_k,z") == -1) ? 1 : -1), DMESG("Parsing error of Initial condition"));


   for(int z = NzLlD; z <= NzLuD; z++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) {  for(int x = NxLlD; x <= NxLuD; x++) {

      const double pos[3] = { X(x), y_k, Z(z) };

      if(plasma->nfields >= n) Field0(x,y_k,z,n) = parser.Eval(pos);

   }}}
   
};



int Init::PerturbationHermitePolynomial(Vlasov *vlasov, int s, double pert, int l) {
  
     const double a =  pow2(0.4); 
     
   for(int x = NxLlD; x <= NxLuD; x++) {
     const double phi_x = Special::HermiteFunction(l, X(x));
     for(int y = NkyLlD; y <= NkyLuD; y++) { for(int z = NzLlD; z<=NzLuD; z++) {
   for(int s = NsLlD; s <= NsLuD; s++) { for(int m = NmLlD; m <= NmLuD; m++) { for(int v = NvLlD; v<=NvLuD; v++) {
      vlasov->f(x,y,z,v,m,s) += epsilon_0 * phi_x *vlasov->f0(x,y,z,v,m,s);
   }}}}}}
  
   return GKC_SUCCESS;
}

