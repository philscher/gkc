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
#include "Special/SpecialMath.h"



FieldsFFT::FieldsFFT(Setup *setup, Grid *grid, Parallel *parallel, FileIO *fileIO, Geometry *geo, FFTSolver *_fft) 
: Fields(setup, grid, parallel, fileIO,  geo), fft(_fft)
{
   
   screenNyquist = setup->get("Fields.screenNyquist", 1);
} 


FieldsFFT::~FieldsFFT() 
{

}

   
void FieldsFFT::solveFieldEquations(CComplex Q     [Nq][NzLD][NkyLD][NxLD],
                                    CComplex Field0[Nq][NzLD][NkyLD][NxLD]) 
{

   // Transform to fourier space (x,ky) -> (kx,ky)
   fft->solve(FFT_X_FIELDS, FFT_FORWARD, FieldsFFT::Q.data());

   // if (phi,Ap,Bp) is solved together, (phi,Bp) are coupled and have to be solved together
   if     ( solveEq &  Field::phi & Field::Bpp ) solveBParallelEquation((A4zz) fft->kXOut.dataZero(), (A4zz) fft->kXIn.dataZero());
   else if( solveEq &  Field::phi              ) solvePoissonEquation  ((A4zz) fft->kXOut.dataZero(), (A4zz) fft->kXIn.dataZero());
   if     ( solveEq &  Field::Ap               ) solveAmpereEquation   ((A4zz) fft->kXOut.dataZero(), (A4zz) fft->kXIn.dataZero());

   // suppresses modes in all fields (needs work)
   //fft->suppressModes(fft->kXIn, Field::phi);

   // replace calcultated field with fixed one if option is set Don't need - or ?
   //if(!(solveEq & Field::phi) && (Nq >= 1)) fft->rXOut(RxLD, RkyLD, RzLD, Field::phi) = Field0(RxLD, RkyLD, RzLD, Field::phi);
   //if(!(solveEq & Field::Ap ) && (Nq >= 2)) fft->rXOut(RxLD, RkyLD, RzLD, Field::Ap ) = Field0(RxLD, RkyLD, RzLD, Field::Ap );
   //if(!(solveEq & Field::Bpp) && (Nq >= 3)) fft->rXOut(RxLD, RkyLD, RzLD, Field::Bp ) = Field0(RxLD, RkyLD, RzLD, Field::Bp );

   // transform back to real-space (kx,ky) -> (x,ky)
   fft->solve(FFT_X_FIELDS, FFT_BACKWARD, FieldsFFT::Field0.data());
   
   return;

}

void FieldsFFT::solvePoissonEquation(CComplex kXOut[Nq][NzLD][NkyLD][FFTSolver::X_NkxL],
                                     CComplex kXIn [Nq][NzLD][NkyLD][FFTSolver::X_NkxL])
{
    // Calculate flux-surface averaging  Note : how to deal with FFT normalization here ?
    CComplex phi_yz[FFTSolver::X_NkxL];
    if(plasma->species[0].doGyro) calcFluxSurfAvrg(kXOut, phi_yz);
    
    // adiabatic response term (if no adiabatic species included n0 = 0)
    const double adiab = plasma->species[0].n0 * pow2(plasma->species[0].q)/plasma->species[0].T0;
    
    omp_for_C2(int z=NzLlD;z<=NzLuD;z++) { for(int y_k=NkyLlD; y_k<=NkyLuD; y_k++) { simd_for(int x_k=fft->K1xLlD; x_k<= fft->K1xLuD; x_k++) {

          // Set kx=0/ky=0 component to zero (gauge freedom)
          if((x_k == 0) && (y_k == 0)) { kXIn[Field::phi][z][y_k][x_k] = 0.e0 ; continue; }
         
          // BUG (Optimize) : need to vectorize this one 
          const double k2_p = fft->k2_p(x_k,y_k,z);
          
          const double lhs    = plasma->debye2 * k2_p + sum_qqnT_1mG0(k2_p) + adiab;
          const CComplex rhs  = (kXOut[Q::rho][z][y_k][x_k])/fft->Norm_X; 
          
          kXIn[Field::phi][z][y_k][x_k] = rhs/lhs;
    
    } } } // z, y_k, x_k
        
    return;    
}



void FieldsFFT::solveAmpereEquation(CComplex kXOut[Nq][NzLD][NkyLD][FFTSolver::X_NkxL],
                                    CComplex kXIn [Nq][NzLD][NkyLD][FFTSolver::X_NkxL])

{
    
  omp_for_C2(int z=NzLlD; z<=NzLuD;z++) {  for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { simd_for(int x_k=fft->K1xLlD; x_k<= fft->K1xLuD;x_k++) {
 
     // Set kx=0/ky=0 component to zero (gauge freedom)
     if((x_k == 0) && (y_k == 0)) { kXIn[Field::Ap][z][y_k][x_k] = 0.e0 ; continue; }
          
     const double k2_p = fft->k2_p(x_k,y_k,z);
           
     const double lhs   =  - k2_p - Yeb * sum_sa2qG0(k2_p);
     const CComplex rhs =  kXOut[Q::jp][z][y_k][x_k]/fft->Norm_X; 
     
     kXIn[Field::Ap][z][y_k][x_k] = rhs/lhs;

  } } } // z, y_k, x_k
      
 return;

}            

void FieldsFFT::solveBParallelEquation(CComplex kXOut[Nq][NzLD][NkyLD][FFTSolver::X_NkxL],
                                       CComplex kXIn [Nq][NzLD][NkyLD][FFTSolver::X_NkxL])
{

  const double adiab = plasma->species[0].n0 * pow2(plasma->species[0].q)/plasma->species[0].T0;
    
  omp_for_C2(int z=NzLlD; z<=NzLuD;z++) { for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { simd_for(int x_k=fft->K1xLlD; x_k<= fft->K1xLuD;x_k++) {
 
     if((x_k == 0) && (y_k == 0)) { kXIn[Field::phi][z][y_k][x_k] =0.; continue; }
         
     const double k2_p = fft->k2_p(x_k,y_k,z);

     const double C_1 = plasma->debye2 * k2_p + sum_qqnT_1mG0(k2_p) + adiab;
     const double C_2 = - sum_qnB_Delta(k2_p);
     const double C_3 = 2./plasma->beta - sum_2TnBB_Delta(k2_p);
          
     const CComplex M_00 = kXOut[Q::rho][z][y_k][x_k];
     const CComplex M_01 = kXOut[Q::jo ][z][y_k][x_k];

     kXIn[Field::phi][z][y_k][x_k] = (C_3 * M_00 - C_2 * M_01) / (C_1 * C_3 - pow2(C_2)) / fft->Norm_X;
     kXIn[Field::Bp ][z][y_k][x_k] = (C_1 * M_01 - C_2 * M_00) / (C_1 * C_3 - pow2(C_2)) / fft->Norm_X;

    } } } // z, y_k, x_k
    
    return;
}


// Note : This is very unpolished and is it also valid in toroidal case ? 
void FieldsFFT::calcFluxSurfAvrg(CComplex kXOut[Nq][NzLD][NkyLD][FFTSolver::X_NkxL],
                                 CComplex phi_yz[FFTSolver::X_NkxL])
{

  if(parallel->Coord[DIR_VMS] != 0) check(-1, DMESG("calcFluxAverg should only be called by XYZ-main nodes"));
   
  phi_yz[:] = 0.; 

  // Note : In FFT the ky=0 components carries the offset over y (integrated value), thus
  //        by divinding through the number of points we get the averaged valued
  for(int z=NzLlD; z<=NzLuD;z++) { if((NkyLlD <= 0) && (NkyLuD >= 0)) { for(int x_k=fft->K1xLlD; x_k<= fft->K1xLuD; x_k++) {
            
     if(x_k == 0) { phi_yz[x_k] = (CComplex) 0.e0 ; continue; }
          
     const double k2_p = fft->k2_p(x_k,0,z);
     
     const double lhs  = plasma->debye2 * k2_p + sum_qqnT_1mG0(k2_p);
     
     // A(x_k, 0) is the sum over y, thus A(x_k, 0)/Ny is the average 
     const CComplex rhs  =  kXOut[Field::phi][z][0][x_k]/((double) (grid->NyGD*Nz));
     
     phi_yz[x_k] = rhs/lhs;
    
  } } }

  // average over z-direction & distribute over y 
//  parallel->collect(phi_yz, OP_SUM, DIR_XYZ);  

  return;

}


// Pretty sure needs to be turned around
void FieldsFFT::gyroFull(Array4C In, Array4C Out,
                         CComplex kXOut[Nq][NzLD][NkyLD][FFTSolver::X_NkxL],
                         CComplex kXIn [Nq][NzLD][NkyLD][FFTSolver::X_NkxL],
                         const int m, const int s, const bool gyroFields)  
{
   
   // Test if gyro-field is OK !     
   fft->solve(FFT_X_FIELDS, FFT_FORWARD, In.data());
   
   // solve for all fields at once
   for(int n = 1; n <= Nq; n++) {
      
        // get therma gyro-radius^2 of species and lambda =  2 x b 
        const double rho_t2  = plasma->species[s].T0 * plasma->species[s].m / (pow2(plasma->species[s].q) * plasma->B0); 
        const double lambda2 = 2. * M[m] * rho_t2;

        // perform gyro-average in Fourier space for rho/phi field
        for(int z=NzLlD; z<=NzLuD;z++) { omp_for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) {  simd_for(int x_k=fft->K1xLlD; x_k<=fft->K1xLuD;x_k++) {
    
         const double k2_p = fft->k2_p(x_k,y_k,z);
          
          kXIn[n][z][y_k][x_k] = kXOut[n][z][y_k][x_k]/fft->Norm_X * ((n != 3) ? j0(sqrt(lambda2 * k2_p)) 
                                                                               : SFL::i1(sqrt(lambda2 * k2_p)));

         // According to GENE code the Nyquiest frequency in x is unphysicall and thus needs to be screened out
         // highest modes explodes in kinetic simulations.
         if((y_k == Nky-1) && screenNyquist) kXIn[n][z][y_k][x_k]  = 0.;
      
       } } }
   }
   
   fft->solve(FFT_X_FIELDS, FFT_BACKWARD, Out.data());
       
   return;

}


void FieldsFFT::gyroFirst(Array4C In, Array4C Out,
                         CComplex kXOut[Nq][NzLD][NkyLD][FFTSolver::X_NkxL],
                         CComplex kXIn [Nq][NzLD][NkyLD][FFTSolver::X_NkxL],
                         const int s, const bool gyroFields)  
{

  if(gyroFields==false) { Out = In; return; };

  fft->solve(FFT_X_FIELDS, FFT_FORWARD, In.data());
           
   // solve for all fields at once
   for(int n = 1; n <= Nq; n++) {

    // get therma gyro-radius^2 of species and lambda =  2 x b 
    const double rho_t2 = plasma->species[s].T0 * plasma->species[s].m / (pow2(plasma->species[s].q) * plasma->B0); 

    for(int z=NzLlD; z<=NzLuD;z++) { omp_for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) {  simd_for(int x_k=fft->K1xLlD; x_k<= fft->K1xLuD;x_k++) {
    
      // perform gyro-average in Fourier space for rho/phi field
      const double k2p_rhoth2 = fft->k2_p(x_k, y_k, z) * rho_t2;
          
      kXIn[n][z][y_k][x_k] = kXOut[n][z][y_k][x_k]/fft->Norm_X * exp(-k2p_rhoth2); 

      if((y_k == Nky-1) && screenNyquist) kXIn[n][z][y_k][x_k]  = 0.;
            
    } } }
   
   }


  fft->solve(FFT_X_FIELDS, FFT_BACKWARD, Out.data());
  
  return;

};


// note : back gyro-average goes over only one field not all !
void FieldsFFT::gyroAverage(Array4C In, Array4C Out, const int m, const int s, const bool gyroFields) 
{
  // check if we should do gyroAvearge
  if     (plasma->species[s].gyroModel == "Drift" ) Out=In;
  else if(plasma->species[s].gyroModel == "Gyro"  ) gyroFull (In, Out, (A4zz) fft->kXOut.dataZero(), (A4zz) fft->kXIn.dataZero(), m, s, gyroFields);  
  else if(plasma->species[s].gyroModel == "Gyro-1") gyroFirst(In, Out, (A4zz) fft->kXOut.dataZero(), (A4zz) fft->kXIn.dataZero(),    s, gyroFields);  
  else   check(-1, DMESG("No such gyro-average Model"));
  
  return;

}


void FieldsFFT::printOn(ostream &output) const
{
  Fields::printOn(output);
  output   << "Fields     |  FFT Solver"  << std::endl;
}


void FieldsFFT::getFieldEnergy(double& phiEnergy, double& ApEnergy, double& BpEnergy)
{
      phiEnergy = 0.; ApEnergy = 0.; BpEnergy = 0.;

      if(parallel->Coord[DIR_VMS] == 0) {
        

        //fft->rXIn(RxLD, RkyLD, RzLD, RFields) = Field0(RxLD,RkyLD, RzLD, RFields);
        fft->solve(FFT_X_FIELDS, FFT_FORWARD, Field0.data());

        // Add only kinetic contributions
        for(int z=NzLlD; z<=NzLuD;z++) { omp_for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { simd_for(int x_k=fft->K1xLlD; x_k<= fft->K1xLuD;x_k++) {
              
              const double k2_p = fft->k2_p(x_k,y_k,z);
              
              if(Nq >= 1) phiEnergy += (plasma->debye2*k2_p + sum_qqnT_1mG0(k2_p)) * pow2(abs(fft->kXOut(x_k,y_k,z,Field::phi)))/fft->Norm_X;
              if(Nq >= 2) ApEnergy  += (k2_p - Yeb * sum_sa2qG0(k2_p))     * pow2(abs(fft->kXOut(x_k,y_k,z,Field::Ap )))/fft->Norm_X;
              if(Nq >= 3) BpEnergy  += 0.;
           
        } } }
      


        // Add Polarization term , see Y.Idomura et al., ...., Kinetic siulations of turnbulence plasmas, Eq. (27)
     // modify field energy term to add adiabatic contributions
        CComplex phi_yz[NxLD]; // RxLD

       /* 
        // short-circuit contribution for ITG mode 
        if(plasma->species[0].doGyro) {
          // muss ich wirklich ueber yz mitteln ? nicht summieren ??
       for(int x = NxLlD; x <= NxLuD; x++) phi_yz(x) = sum(fields->phi(x, RyLD, RzLD, 1, 1))/((double) Ny*Nz);
//       for(int x = NxLlD; x <= NxLuD; x++) phi_yz(x) = sum(fields->phi(x, RyLD, RzLD, 1, 1));
            parallel->collect(phi_yz, OP_SUM, DIR_YZ);
        }
      
        // add Adiabatic contributions
        const double adiab = plasma->species[0].n0 * pow2(plasma->species[0].q)/plasma->species[0].T0;
        for(int x = NxLlD; x <= NxLuD; x++)  phiEnergy += adiab * abs(__sec_reduce(pow2(phi(x,RkyLD,RzLD, 1, 1) - phi_yz[x-NxLlD])));
        * */ 
     
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


///////////////////////////////////////// Some Helper functions ///////////////////////////////////////////



double FieldsFFT::sum_qqnT_1mG0(const double k2_p) 
{
   double g0 = 0.;

   for(int s = NsGlD; s <= NsGuD; s++) {

      const double qqnT   = plasma->species[s].n0 * pow2(plasma->species[s].q)/plasma->species[s].T0;
      const double rho_t2 = plasma->species[s].T0  * plasma->species[s].m / pow2(plasma->species[s].q * plasma->B0);
       
      //g0 += qqnT * SpecialMath::_1mGamma0_Pade( rho_t2 * k2_p);
      g0 += qqnT * SpecialMath::_1mGamma0( rho_t2 * k2_p);
        
   }
   return g0;
};
  
double FieldsFFT::sum_sa2qG0(const double kp_2) 
{
    double g0 = 0.;
    
    #pragma simd
    for(int s = NsGlD; s <= NsGuD; s++) {

       const double sa2q   = plasma->species[s].sigma * pow2(plasma->species[s].alpha) * plasma->species[s].q;
       const double rho_t2 = plasma->species[s].T0 * plasma->species[s].m / (pow2(plasma->species[s].q) * plasma->B0);

       const double b = rho_t2 * kp_2;
       
       g0  += sa2q * (1.e0 - SpecialMath::_1mGamma0_Pade(b));
    }
    return g0;
};


double FieldsFFT::sum_qnB_Delta(const double k2_p) 
{

   double g0 = 0.;
   
   for(int s = NsGlD; s <= NsGuD; s++) {

      const double qn    = plasma->species[s].q * plasma->species[s].n0;
      const double rho_t2 = plasma->species[s].T0 * plasma->species[s].m / (pow2(plasma->species[s].q) * plasma->B0);

      g0 += qn * SpecialMath::Delta(rho_t2 * k2_p);
   }
   return g0/plasma->B0;
};

double FieldsFFT::sum_2TnBB_Delta(const double k2_p) 
{
   double g0 = 0.;
   
   for(int s = NsGlD; s <= NsGuD; s++) {

      const double Tn   = plasma->species[s].T0 * plasma->species[s].n0;
      const double rho_t2 = plasma->species[s].T0 * plasma->species[s].m / (pow2(plasma->species[s].q) * plasma->B0);

      g0 += Tn * SpecialMath::Delta(rho_t2 * k2_p);
   }

   return 2. * g0 /pow2(plasma->B0) ;
};
