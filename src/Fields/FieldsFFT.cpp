/*
 * =====================================================================================
 *
 *       Filename: FieldsFFT.cpp
 *
 *    Description: Fields Solver in Fourier space
 *
 *         Author: Paul P. Hilscher (2009-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include "FieldsFFT.h"
#include "SpecialMath.h"

FieldsFFT::FieldsFFT(Setup *setup, Grid *grid, Parallel *parallel, FileIO *fileIO, Geometry<HELIOS_GEOMETRY> *geo, FFTSolver *fft) 
: Fields(setup, grid, parallel, fileIO,  geo), Fourier3D(setup, grid, fft)
{
	phi_yz.resize(fft->Rk1xL); phi_yz = 0.;

   screenNyquist = setup->get("Fields.screenNyquist", 1);
} 


FieldsFFT::~FieldsFFT() {

}

Array4z FieldsFFT::solveFieldEquations(Array4z Q, Timing timing) {
   
  
   // need to set other Q to zero to avoid accumulation of errors
   // is this necessary ? avoid accumulation of numerical errors in FFT 
   fft->rXOut = 0.;

   // transform to fourier space
   fft->rXIn(RxLD, RkyLD, RzLD, RFields) = Q(RxLD, RkyLD, RzLD, RFields);
   fft->solve(FFT_X, FFT_FORWARD, FFT_FIELDS);
    
   // for 3 fields phi and B_perp are coupled and solved together
   if     ( solveEq &  Field::phi & Field::Bpp ) solveBParallelEquation(Q(RxLD, RkyLD, RzLD, Field::Bp ), timing);
   else if( solveEq &  Field::phi              ) solvePoissonEquation  (Q(RxLD, RkyLD, RzLD, Field::phi), timing);
   if     ( solveEq &  Field::Ap               ) solveAmpereEquation   (Q(RxLD, RkyLD, RzLD, Field::Ap ), timing);

    
   // suppresses modes in all fields
   suppressModes(fft->kXIn, Field::phi);

   fft->solve(FFT_X, FFT_BACKWARD, FFT_FIELDS);
  
   // replace calcultated field with fixed one if option is set Don't need - or ?
   //if(!(solveEq & Field::phi) && (plasma->nfields >= 1)) fft->rXOut(RxLD, RkyLD, RzLD, Field::phi) = Field0(RxLD, RkyLD, RzLD, Field::phi);
   //if(!(solveEq & Field::Ap ) && (plasma->nfields >= 2)) fft->rXOut(RxLD, RkyLD, RzLD, Field::Ap ) = Field0(RxLD, RkyLD, RzLD, Field::Ap );
   //if(!(solveEq & Field::Bpp) && (plasma->nfields >= 3)) fft->rXOut(RxLD, RkyLD, RzLD, Field::Bp ) = Field0(RxLD, RkyLD, RzLD, Field::Bp );

    
   return fft->rXOut(RxLD, RkyLD, RzLD, RFields);


}



Array3z FieldsFFT::solvePoissonEquation(Array3z rho, Timing timing)
{
    // Calculate flux-surface averaging
    // Note : how to deal with FFT normalization here ?
    if(plasma->species(0).doGyro) phi_yz = calcFluxSurfAvrg(fft->kXOut);
    
    // adiabatic response term (if no adiabatic species included n0 = 0)
    const double adiab = plasma->species(0).n0 * pow2(plasma->species(0).q)/plasma->species(0).T0;
    
    // solve PoissonEq. in Fourier space, using periodic boundary conditions
    for(int z=NzLlD;z<=NzLuD;z++) { omp_for(int y_k=NkyLlD; y_k<=NkyLuD;y_k++) { for(int x_k=fft->K1xLlD; x_k<= fft->K1xLuD;x_k++) {

          // Set kx=0/ky=0 component to zero (gauge freedom)
          if((x_k == 0) && (y_k == 0)) { fft->kXIn(x_k,y_k,z, Field::phi) = (cmplxd) 0.e0 ; continue; }
          
          const double k2_p = fft->k2_p(x_k,y_k,z);
          
          const double lhs  = plasma->debye2 * k2_p + sum_qqnT_1mG0(k2_p) + adiab;
          const cmplxd rhs  =  (fft->kXOut(x_k,y_k,z, Q::rho) + adiab * phi_yz(x_k))/fft->Norm_X; 
          
          fft->kXIn(x_k,y_k,z, Field::phi) = rhs/lhs;
         
    } } }
        
    // note we requre this one !! maybe from Analysis or so ?
    return rho;    
}


Array3z FieldsFFT::solveAmpereEquation(Array3z j, Timing timing)
{
    
  for(int z=NzLlD; z<=NzLuD;z++) {  omp_for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int x_k=fft->K1xLlD; x_k<= fft->K1xLuD;x_k++) {
 
     if((x_k == 0) && (y_k == 0)) { fft->kXIn(x_k,y_k,z, Field::Ap) = (cmplxd) 0.e0; continue; }
          
     const double k2_p = fft->k2_p(x_k,y_k,z);
           
     const double lhs =  - k2_p - Yeb * sum_sa2qG0(k2_p);
     const cmplxd rhs =  fft->kXOut(x_k,y_k,z, Q::jp)/fft->Norm_X;
          
     fft->kXIn(x_k,y_k,z, Field::Ap ) = rhs/lhs;

  } } }
      
 return j;

}            

// Note : what about adiabtic/flux-surface averaging ? Or are kinetic electrons required ?
Array3z FieldsFFT::solveBParallelEquation(Array3z phi, Timing timing) 
{

  const double adiab = plasma->species(0).n0 * pow2(plasma->species(0).q)/plasma->species(0).T0;
    
  for(int z=NzLlD; z<=NzLuD;z++) { omp_for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int x_k=fft->K1xLlD; x_k<= fft->K1xLuD;x_k++) {
 
     if((x_k == 0) && (y_k == 0)) { fft->kXIn(x_k,y_k,z, Field::Bp) = (cmplxd) 0.e0; continue; }
         
     const double k2_p = fft->k2_p(x_k,y_k,z);

     const double C_1 = plasma->debye2 * k2_p + sum_qqnT_1mG0(k2_p) + adiab;
     const double C_2 = - sum_qnB_Delta(k2_p);
     const double C_3 = 2./plasma->beta - sum_2TnBB_Delta(k2_p);
          
     const cmplxd M_00 = fft->kXOut(x_k, y_k, z, Q::rho);
     const cmplxd M_01 = fft->kXOut(x_k, y_k, z, Q::jo );

     fft->kXIn(x_k, y_k, z, Field::phi) = (C_3 * M_00 - C_2 * M_01) / (C_1 * C_3 - pow2(C_2)) / fft->Norm_X;
     fft->kXIn(x_k, y_k, z, Field::Bp ) = (C_1 * M_01 - C_2 * M_00) / (C_1 * C_3 - pow2(C_2)) / fft->Norm_X;

    } } }
    
    return phi;
}


// Note : This is very unpolished and is it also valid in toroidal case ? 
Array1z FieldsFFT::calcFluxSurfAvrg(Array4z kXOut) 
{

  if(parallel->Coord(DIR_VMS) != 0) check(-1, DMESG("calcFluxAverg should only be called by XYZ-main nodes"));
   
  phi_yz = 0.; 

  // Note : In FFT the ky=0 components carries the offset over y (integrated value), thus
  //        by divinding through the number of points we get the averaged valued
  for(int z=NzLlD; z<=NzLuD;z++) { if((NkyLlD <= 0) && (NkyLuD >= 0)) { for(int x_k=fft->K1xLlD; x_k<= fft->K1xLuD; x_k++) {
            
     if(x_k == 0) { phi_yz(x_k) = (cmplxd) 0.e0 ; continue; }
          
     const double k2_p = fft->k2_p(x_k,0,z);
     
     const double lhs  = plasma->debye2 * k2_p + sum_qqnT_1mG0(k2_p);
     // A(x_k, 0) is the sum over y, thus A(x_k, 0)/Ny is the average 
     const cmplxd rhs  =  kXOut(x_k,0,z, Field::phi)/((double) (grid->NyGD*Nz));
     phi_yz(x_k) = rhs/lhs;
    
  } } }

  // average over z-direction & distribute over y 
  parallel->collect(phi_yz, OP_SUM, DIR_XYZ);  

  return phi_yz;

}


Array4z FieldsFFT::gyroFull(Array4z A4, int m, int s, const int nField2, bool gyroField)  
{
   if(!(do_gyro && (parallel->Coord(DIR_V) == 0) && plasma->species(s).doGyro)) return A4;
        
   fft->rXIn(RxLD, RkyLD, RzLD, RFields) = A4(RxLD, RkyLD, RzLD, RFields);
   fft->solve(FFT_X, FFT_FORWARD, FFT_FIELDS);
        
   for(int nField = 1; nField <= plasma->nfields; nField++) {

      
      // get therma gyro-radius^2 of species and lambda =  2 x b 
     	const double rho_t2 = plasma->species(s).T0 * plasma->species(s).m / (pow2(plasma->species(s).q) * plasma->B0); 
     	const double lambda2 = 2. * M(m) * rho_t2;

      // perform gyro-average in Fourier space for rho/phi field
      for(int z=NzLlD; z<=NzLuD;z++) { omp_for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) {  for(int x_k=fft->K1xLlD; x_k<= fft->K1xLuD;x_k++) {
    
         const double k2_p = fft->k2_p(x_k,y_k,z);
         fft->kXIn(x_k,y_k,z, nField) = fft->kXOut(x_k,y_k,z, nField)/fft->Norm_X *  ((nField != 3) ? SpecialMath::BesselJ0(sqrt(lambda2 * k2_p)) 
                                                                                                    : SpecialMath::BesselI1(sqrt(lambda2 * k2_p)));
           
         // According to GENE code the Nyquiest frequency in x is unphysicall and thus needs to be screened out * geo->B0k(k_x,k_y,z)
         // highest modes explodes in kinetic simulations.
         if((y_k == Nky-1) && screenNyquist) fft->kXIn(x_k,y_k,z, nField) = 0.;

      } } }
   }
        
   fft->solve(FFT_X, FFT_BACKWARD, FFT_FIELDS);
       
   return fft->rXOut(RxLD, RkyLD, RzLD, RFields);

}

// This is depracated
Array4z FieldsFFT::gyroFirst(Array4z A4, int m, int s, const int nField2, bool gyroField)  
{

  if(gyroField==false) return A4;
        
  fft->rXIn(RxLD, RkyLD, RzLD, RFields) = A4(RxLD, RkyLD, RzLD, RFields);
  fft->solve(FFT_X, FFT_FORWARD, FFT_FIELDS);
           
  for(int nField = 1; nField <= plasma->nfields; nField++) {

    // get therma gyro-radius^2 of species and lambda =  2 x b 
    const double rho_t2 = plasma->species(s).T0 * plasma->species(s).m / (pow2(plasma->species(s).q) * plasma->B0); 

    for(int z=NzLlD; z<=NzLuD;z++) { for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) {  for(int x_k=fft->K1xLlD; x_k<= fft->K1xLuD;x_k++) {
    
      // perform gyro-average in Fourier space for rho/phi field
      const double k2p_rho2 = fft->k2_p(x_k, y_k, z) * rho_t2;

      if(m==2) fft->kXIn(x_k,y_k,z, nField) = k2p_rho2 * fft->kXOut(x_k,y_k,z, nField)/fft->Norm_X *  exp(-k2p_rho2);
      else     fft->kXIn(x_k,y_k,z, nField) =            fft->kXOut(x_k,y_k,z, nField)/fft->Norm_X *  exp(-k2p_rho2);

      // we set Nyquiest frequency to zero, as the ordering of the fftw with only + k_Ny, and non-linearity required both +k_Ny AND -k_Ny
      if((y_k == Nky-1) && screenNyquist) fft->kXIn(x_k,y_k,z, nField) = 0.;
            
    } } }
	
  }

  fft->solve(FFT_X, FFT_BACKWARD, FFT_FIELDS);
  return fft->rXOut(RxLD, RkyLD, RzLD, RFields);

};


Array4z FieldsFFT::gyroAverage(Array4z A4, int m, int s, const int nField, bool gyroField) 
{

  // check if we should do gyroAvearge
  if     (plasma->species(s).gyroModel == "Drift" ) return A4;
  else if(plasma->species(s).gyroModel == "Gyro"  ) return gyroFull (A4, m, s, nField, gyroField);  
  else if(plasma->species(s).gyroModel == "Gyro-1") return gyroFirst(A4, m, s, nField, gyroField); 
  else   check(-1, DMESG("No such gyro-average Model"));
  
  return A4;

}

// @todo also show solver
void FieldsFFT::printOn(ostream &output) const
{
  Fields::printOn(output);
  output   << "Fields     |  FFT "  << std::endl;
}


void FieldsFFT::calculateFieldEnergy(Array4z Q, double& phiEnergy, double& ApEnergy, double& BpEnergy)
{
      phiEnergy = 0.; ApEnergy = 0.; BpEnergy = 0.;

      if(parallel->Coord(DIR_VMS) == 0) {
        

        fft->rXIn(RxLD, RkyLD, RzLD, RFields) = Field0(RxLD,RkyLD, RzLD, RFields);
        fft->solve(FFT_X, FFT_FORWARD, FFT_FIELDS);

        // Add only kinetic contributions
        for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int z=NzLlD; z<=NzLuD;z++) { for(int x_k=fft->K1xLlD; x_k<= fft->K1xLuD;x_k++) {
              
              const double k2_p = fft->k2_p(x_k,y_k,z);
              
              if(plasma->nfields >= 1) phiEnergy += (plasma->debye2*k2_p + sum_qqnT_1mG0(k2_p)) * pow2(abs(fft->kXOut(x_k,y_k,z,Field::phi)))/fft->Norm_X;
              if(plasma->nfields >= 2) ApEnergy  += (k2_p - Yeb * sum_sa2qG0(k2_p))     * pow2(abs(fft->kXOut(x_k,y_k,z,Field::Ap )))/fft->Norm_X;
              if(plasma->nfields >= 3) BpEnergy  += 0.;
           
        } } }
      


        // Add Polarization term , see Y.Idomura et al., ...., Kinetic siulations of turnbulence plasmas, Eq. (27)
     // modify field energy term to add adiabatic contributions

        // short-circuit contribution for ITG mode 
        Array1z phi_yz(RxLD); phi_yz = 0.;
       /* 
        if(plasma->species(0).doGyro) {
          // muss ich wirklich ueber yz mitteln ? nicht summieren ??
	    for(int x = NxLlD; x <= NxLuD; x++) phi_yz(x) = sum(fields->phi(x, RyLD, RzLD, 1, 1))/((double) Ny*Nz);
//	    for(int x = NxLlD; x <= NxLuD; x++) phi_yz(x) = sum(fields->phi(x, RyLD, RzLD, 1, 1));
            parallel->collect(phi_yz, OP_SUM, DIR_YZ);
        }
      
        * */ 
        // add Adiabatic contributions
        const double adiab = plasma->species(0).n0 * pow2(plasma->species(0).q)/plasma->species(0).T0;
        for(int x = NxLlD; x <= NxLuD; x++)  phiEnergy += adiab * abs(sum(pow2(phi(x,RkyLD,RzLD, 1, 1) - phi_yz(x))));
     
      }
        // still HDF-5 is a bit buggy and we need to distribute the value over all nodes
       //phiEnergy =  parallel->collect(4. * M_PI * phiEnergy * scaleXYZ / (8.e0 * M_PI) / (initialEkin(TOTAL) == 0. ? 1. : initialEkin(TOTAL)), OP_SUM, DIR_XYZ);
       //ApEnergy  =  parallel->collect(4. * M_PI *  ApEnergy * scaleXYZ / (8.e0 * M_PI) / (initialEkin(TOTAL) == 0. ? 1. : initialEkin(TOTAL)), OP_SUM, DIR_XYZ);
       //phiEnergy =  parallel->collect(4. * M_PI * phiEnergy * scaleXYZ / (8.e0 * M_PI), OP_SUM, DIR_XYZ);
       //ApEnergy  =  parallel->collect(4. * M_PI *  ApEnergy * scaleXYZ / (8.e0 * M_PI), OP_SUM, DIR_XYZ);
       //
       phiEnergy =  abs( parallel->collect(4. * M_PI * phiEnergy * grid->dXYZ / (8.e0 * M_PI), OP_SUM, DIR_ALL) ); 
       ApEnergy  =  abs( parallel->collect(4. * M_PI *  ApEnergy * grid->dXYZ / (8.e0 * M_PI), OP_SUM, DIR_ALL) );
       BpEnergy  = 0.;
      
}


//////////////// Some Helper functions ////////////////////


double FieldsFFT::sum_qqnT_1mG0(const double k2_p) 
{
   double g0 = 0.;
   for(int s = NsGlD; s <= NsGuD; s++) {

      const double qqnT   = plasma->species(s).n0 * pow2(plasma->species(s).q)/plasma->species(s).T0;
      const double rho_t2 = plasma->species(s).T0  * plasma->species(s).m / pow2(plasma->species(s).q * plasma->B0);
       
      // doGyro = true
      if(plasma->species(s).doGyro == true) g0 += qqnT * SpecialMath::_1mGamma0( rho_t2 * k2_p);
      else g0 += qqnT * (rho_t2 * k2_p);
        
      //g0 += qqnT * _1mGamma0( rho_t2 * k2_p, plasma->species(s).doGyro);
   }
   return g0;
};
  

double FieldsFFT::sum_sa2qG0(const double kp_2) 
{
    double g0 = 0.;
    for(int s = NsGlD; s <= NsGuD; s++) {
       const double sa2q   = plasma->species(s).sigma * pow2(plasma->species(s).alpha) * plasma->species(s).q;
       const double rho_t2 = plasma->species(s).T0 * plasma->species(s).m / (pow2(plasma->species(s).q) * plasma->B0);

       const double b = rho_t2 * kp_2;
       //g0 += sa2q * (1.e0 - _1mGamma0( rho_t2 * kp_2, plasma->species(s).doGyro));
       //if(plasma->species(s).doGyro == true) g0 += sa2q * I0(b)*exp(-b);
       if(plasma->species(s).doGyro == true) g0  += sa2q * (1.e0 - SpecialMath::_1mGamma0(b));
       else g0 += sa2q;
    }
    return g0;
};


double FieldsFFT::sum_qnB_Delta(const double k2_p) 
{

   double g0 = 0.;
   for(int s = NsGlD; s <= NsGuD; s++) {

      const double qn    = plasma->species(s).q * plasma->species(s).n0;
      const double rho_t2 = plasma->species(s).T0 * plasma->species(s).m / (pow2(plasma->species(s).q) * plasma->B0);

      g0 += qn * SpecialMath::Delta(rho_t2 * k2_p);
   }
   return g0/plasma->B0;
};

double FieldsFFT::sum_2TnBB_Delta(const double k2_p) 
{
   double g0 = 0.;
   for(int s = NsGlD; s <= NsGuD; s++) {

      const double Tn   = plasma->species(s).T0 * plasma->species(s).n0;
      const double rho_t2 = plasma->species(s).T0 * plasma->species(s).m / (pow2(plasma->species(s).q) * plasma->B0);

      g0 += Tn * SpecialMath::Delta(rho_t2 * k2_p);
   }

   return 2. * g0 /pow2(plasma->B0) ;
};
