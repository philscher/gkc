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

  return;
} 


FieldsFFT::~FieldsFFT() 
{

}

   
void FieldsFFT::solveFieldEquations(const CComplex Q     [Nq][NzLD][Nky][NxLD],
                                          CComplex Field0[Nq][NzLD][Nky][NxLD]) 
{
  // Transform to fourier space (x,ky) -> (kx,ky)
  #pragma omp single
  fft->solve(FFT_Type::X_FIELDS, FFT_Sign::Forward, ArrayField0.data((CComplex *) Q));

  // if (phi,Ap,Bp) is solved together, (phi,Bp) are coupled and have to be solved together
  if     ( solveEq & Field::Iphi & Field::IBp ) solveBParallelEquation((A4zz) fft->kXOut, (A4zz) fft->kXIn);
  else if( solveEq & Field::Iphi              ) solvePoissonEquation  ((A4zz) fft->kXOut, (A4zz) fft->kXIn);
  if     ( solveEq & Field::IAp               ) solveAmpereEquation   ((A4zz) fft->kXOut, (A4zz) fft->kXIn);
  #pragma omp barrier

  // suppresses modes in all fields (needs work)
  // fft->suppressModes(fft->kXIn, Field::phi);

  // replace calcultated field with fixed one if option is set Don't need - or ?
  //if(!(solveEq & Field::phi) && (Nq >= 1)) fft->rXOut(RxLD, RkyLD, RzLD, Field::phi) = Field0(RxLD, RkyLD, RzLD, Field::phi);
  //if(!(solveEq & Field::Ap ) && (Nq >= 2)) fft->rXOut(RxLD, RkyLD, RzLD, Field::Ap ) = Field0(RxLD, RkyLD, RzLD, Field::Ap );
  //if(!(solveEq & Field::Bpp) && (Nq >= 3)) fft->rXOut(RxLD, RkyLD, RzLD, Field::Bp ) = Field0(RxLD, RkyLD, RzLD, Field::Bp );

  // transform back to real-space (kx,ky) -> (x,ky)
  #pragma omp single
  fft->solve(FFT_Type::X_FIELDS, FFT_Sign::Backward, ArrayField0.data((CComplex *) Field0));


  return;

}

void FieldsFFT::solvePoissonEquation(CComplex kXOut[Nq][NzLD][Nky][FFTSolver::X_NkxL],
                                     CComplex kXIn [Nq][NzLD][Nky][FFTSolver::X_NkxL])
{
  // Calculate flux-surface averaging  Note : how to deal with FFT normalization here ?
  // OpenMP how to set this variable to be shared ? copyprivate(phi_yz)
  CComplex phi_yz[Nx]; phi_yz[:] = 0.;

  #pragma omp single
  {
    if(species[0].doGyro) calcFluxSurfAvrg(kXOut, phi_yz);
  }

  // adiabatic response term (if no adiabatic species included n0 = 0)
  const double adiab = species[0].n0 * pow2(species[0].q)/species[0].T0;
    
  #pragma omp for collapse(2) nowait
  for(int z = NzLlD; z <= NzLuD; z++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { 
    
  simd_for(int x_k = fft->K1xLlD; x_k <= fft->K1xLuD; x_k++) {

    // Set kx=0/ky=0 component to zero (gauge freedom)
    if((x_k == 0) && (y_k == 0)) { kXIn[Field::phi][z][y_k][x_k] = 0.e0 ; continue; }
         
    // BUG (Optimize) : need to vectorize this one 
    const double k2_p = fft->k2_p(x_k,y_k,z);
         
    // adiabatic term \adia ( \phi - <\phi>_{yz}), we shift flux averaging term <\phi>_{yz}
    // to rhs due to FFT normalization
    const double lhs    = plasma->debye2 * k2_p + sum_qqnT_1mG0(k2_p) + adiab;
    const CComplex rhs  = (kXOut[Field::phi][z][y_k][x_k] + adiab*phi_yz[x_k])/fft->Norm_X; 
          
    kXIn[Field::phi][z][y_k][x_k] = rhs/lhs;
         
    // where to remove ? Here or @ gyroFull ? ... I guess better here ... !
    if(( (y_k == Nky-1) || (x_k == Nx/2)) && screenNyquist) kXIn[Field::phi][z][y_k][x_k]  = 0.;
    
  } } } // z, y_k, x_k
        
  return;    

}



void FieldsFFT::solveAmpereEquation(CComplex kXOut[Nq][NzLD][Nky][FFTSolver::X_NkxL],
                                    CComplex kXIn [Nq][NzLD][Nky][FFTSolver::X_NkxL])

{
  #pragma omp for collapse(2) nowait
  for(int z = NzLlD; z <= NzLuD; z++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { 
    
  simd_for(int x_k = fft->K1xLlD; x_k <= fft->K1xLuD; x_k++) {
 
    // Set kx=0/ky=0 component to zero (gauge freedom)
    if((x_k == 0) && (y_k == 0)) { kXIn[Field::Ap][z][y_k][x_k] = 0.e0 ; continue; }
          
    const double k2_p = fft->k2_p(x_k,y_k,z);
           
    const double lhs   =  - k2_p - Yeb * sum_sa2qG0(k2_p);
    const CComplex rhs =  kXOut[Field::Ap][z][y_k][x_k]/fft->Norm_X; 
     
    kXIn[Field::Ap][z][y_k][x_k] = rhs/lhs;

  } } } // z, y_k, x_k
      
  return;

}            

void FieldsFFT::solveBParallelEquation(CComplex kXOut[Nq][NzLD][Nky][FFTSolver::X_NkxL],
                                       CComplex kXIn [Nq][NzLD][Nky][FFTSolver::X_NkxL])
{

  const double adiab = species[0].n0 * pow2(species[0].q)/species[0].T0;
    
  #pragma omp for collapse(2) nowait
  for(int z = NzLlD; z <= NzLuD; z++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { 
   
  simd_for(int x_k = fft->K1xLlD; x_k <= fft->K1xLuD; x_k++) {
 
    if((x_k == 0) && (y_k == 0)) { kXIn[Field::phi][z][y_k][x_k] =0.; continue; }
         
    const double k2_p = fft->k2_p(x_k,y_k,z);

    const double C_1 = plasma->debye2 * k2_p + sum_qqnT_1mG0(k2_p) + adiab;
    const double C_2 = - sum_qnB_Delta(k2_p);
    const double C_3 = 2./plasma->beta - sum_2TnBB_Delta(k2_p);
          
    const CComplex M_00 = kXOut[Field::phi][z][y_k][x_k];
    const CComplex M_01 = kXOut[Field::Bp ][z][y_k][x_k];

    kXIn[Field::phi][z][y_k][x_k] = (C_3 * M_00 - C_2 * M_01) / (C_1 * C_3 - pow2(C_2)) / fft->Norm_X;
    kXIn[Field::Bp ][z][y_k][x_k] = (C_1 * M_01 - C_2 * M_00) / (C_1 * C_3 - pow2(C_2)) / fft->Norm_X;

  } } } // z, y_k, x_k
    
  return;
}


// Note : is it also valid in toroidal case ? 
void FieldsFFT::calcFluxSurfAvrg(CComplex kXOut[Nq][NzLD][Nky][FFTSolver::X_NkxL],
                                 CComplex phi_yz[Nx])
{
  const double _kw_NxNy = 1./((2.*Nky-2.)*Nx*Nz)  ; // Number of poloidal points in real space

  // Note : In FFT the ky=0 components carries the offset over y (integrated value), thus
  //        by divinding through the number of points we get the averaged valued
  for(int z = NzLlD; z <= NzLuD; z++) { if(NkyLlD == 0) { 
    
  for(int x_k = fft->K1xLlD; x_k <= fft->K1xLuD; x_k++) {
            
    if(x_k == 0) { phi_yz[x_k] = 0.e0 ; continue; }
         
    // k_y < A >_y = 0 (thus no need to include)
    const double k2_p  = fft->k2_p(x_k, 0, z);
     
    const double lhs   = plasma->debye2 * k2_p + sum_qqnT_1mG0(k2_p);
    // A(x_k, 0) is the sum over y, thus A(x_k, 0)/Ny is the average 
    const CComplex rhs =  kXOut[Field::phi][z][0][x_k] * _kw_NxNy;
     
    phi_yz[x_k] = rhs/lhs;
    
  } } }

  // average over z-direction
  parallel->reduce(phi_yz, Op::sum, DIR_Z, Nx);  

  return;

}

void FieldsFFT::doubleGyroExp(const CComplex In [NzLD][Nky][NxLD], 
                                    CComplex Out[NzLD][Nky][NxLD], const int m, const int s)
{
       
  fft->solve(FFT_Type::X_FIELDS, FFT_Sign::Forward, (void *) &In[0][0][0]);

  [=](CComplex kXOut[Nq][NzLD][Nky][FFTSolver::X_NkxL],
      CComplex kXIn [Nq][NzLD][Nky][FFTSolver::X_NkxL])
  {
       
    const double qqnT   = species[s].n0 * pow2(species[s].q)/species[s].T0;
    const double rho_t2 = species[s].T0  * species[s].m / pow2(species[s].q * plasma->B0);
   
    for(int q = 0; q < Nq; q++) {
      
    //#pragma omp for collapse(2)
    for(int z = NzLlD; z <= NzLuD; z++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { 
      
    for(int x_k = fft->K1xLlD; x_k <= fft->K1xLuD; x_k++) {
    
      const double k2_p = fft->k2_p(x_k,y_k,z);
          
      kXIn[q][z][y_k][x_k] = kXOut[q][z][y_k][x_k]/fft->Norm_X * ((m == 0) ? SFL::i0e(rho_t2 * k2_p) 
                                                                           : SpecialMath::Delta_1(rho_t2 * k2_p));

    } } } // z, y_k, x
   
    } // q

  } ((A4zz) fft->kXOut, (A4zz) fft->kXIn);
       
  fft->solve(FFT_Type::X_FIELDS, FFT_Sign::Backward, &Out[0][0][0]);

  return;

}

void FieldsFFT::gyroFull(const CComplex In   [Nq][NzLD][Nky][NxLD             ], 
                               CComplex Out  [Nq][NzLD][Nky][NxLD             ],
                               CComplex kXOut[Nq][NzLD][Nky][FFTSolver::X_NkxL],
                               CComplex kXIn [Nq][NzLD][Nky][FFTSolver::X_NkxL],
                         const int m, const int s, bool stack)  
{
  
  #pragma omp single
  fft->solve(FFT_Type::X_FIELDS, FFT_Sign::Forward, stack ? (void *) &In[0][0][0][0] : (void *) &In[0][NzLlD][0][NxLlD]);

  // get therma gyro-radius^2 of species and lambda =  2 x b 
  const double rho_t2  = species[s].T0 * species[s].m / (pow2(species[s].q) * plasma->B0); 
  const double lambda2 = 2. * M[m] * rho_t2;

  // solve for all fields at once
  for(int q = 0; q < Nq; q++) {
      
  #pragma omp for collapse(2)
  for(int z = NzLlD; z <= NzLuD; z++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { 

  for(int x_k = fft->K1xLlD; x_k <= fft->K1xLuD; x_k++) {
    
    const double k2_p = fft->k2_p(x_k,y_k,z);
          
    kXIn[q][z][y_k][x_k] = kXOut[q][z][y_k][x_k]/fft->Norm_X * ((q != 3) ? j0(sqrt(lambda2 * k2_p)) 
                                                                         : SFL::i1(sqrt(lambda2 * k2_p)));
                                                                         //?: SFL::j1(sqrt(lambda2 * k2_p)));

  } } }

  } // q
   
  #pragma omp single
  fft->solve(FFT_Type::X_FIELDS, FFT_Sign::Backward, stack ? &Out[0][0][0][0] : &Out[0][NzLlD][0][NxLlD]);
 
  return;
}


void FieldsFFT::gyroFirst(const CComplex In   [Nq][NzLD][Nky][NxLD], 
                                CComplex Out  [Nq][NzLD][Nky][NxLD],
                                CComplex kXOut[Nq][NzLD][Nky][FFTSolver::X_NkxL],
                                CComplex kXIn [Nq][NzLD][Nky][FFTSolver::X_NkxL],
                          const int s, const bool gyroFields)  
{

  if(gyroFields==false) {  

     Out[0:Nq][NzLlD:NzLD][NkyLlD:Nky][NxLlD:NxLD] = In[0:Nq][NzLlD:NzLD][NkyLlD:NkyLD][NxLlD:NxLD];
     return; 
  };
   
  fft->solve(FFT_Type::X_FIELDS, FFT_Sign::Forward, (void *) &In[0][NzLlD][NkyLlD][NxLlD]);
           
  // solve for all fields at once
  for(int q = 0; q < Nq; q++) {

  // get therma gyro-radius^2 of species and lambda =  2 x b 
  const double rho_t2 = species[s].T0 * species[s].m / (pow2(species[s].q) * plasma->B0); 

  #pragma omp for collapse(2)
  for(int z=NzLlD; z<=NzLuD;z++) { for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) {  simd_for(int x_k=fft->K1xLlD; x_k<= fft->K1xLuD;x_k++) {
    
    // perform gyro-average in Fourier space for rho/phi field
    const double k2p_rhoth2 = fft->k2_p(x_k, y_k, z) * rho_t2;
          
    kXIn[q][z][y_k][x_k] = kXOut[q][z][y_k][x_k]/fft->Norm_X * exp(-k2p_rhoth2); 

    if(((y_k == Nky-1) || (x_k == Nx/2)) && screenNyquist) kXIn[q][z][y_k][x_k]  = 0.;
            
  } } }
   
  }


  fft->solve(FFT_Type::X_FIELDS, FFT_Sign::Backward, &Out[0][NzLlD][NkyLlD][NxLlD]);
  
  return;

};


// note : back gyro-average goes over only one field not all !
void FieldsFFT::gyroAverage(const CComplex In [Nq][NzLD][Nky][NxLD], CComplex Out[Nq][NzLD][NkyLD][NxLD],
                            const int m, const int s, const bool gyroFields, const bool stack)
{
  if     (species[s].gyroModel == "Drift" ) {
       
      Out[0:Nq][NzLlD:NzLD][NkyLlD:Nky][NxLlD:NxLD]
   =   In[0:Nq][NzLlD:NzLD][NkyLlD:Nky][NxLlD:NxLD];
  }
  else if(species[s].gyroModel == "Gyro"  ) gyroFull (In, Out, (A4zz) fft->kXOut, (A4zz) fft->kXIn, m, s, stack);  
  else if(species[s].gyroModel == "Gyro-1") gyroFirst(In, Out, (A4zz) fft->kXOut, (A4zz) fft->kXIn,    s, gyroFields);  
  else   check(-1, DMESG("No such gyro-average Model"));
  
  return;

}


void FieldsFFT::printOn(std::ostream &output) const
{
  Fields::printOn(output);
  output   << "Fields     |  FFT Solver"  << std::endl;
}


void FieldsFFT::getFieldEnergy(double& phiEnergy, double& ApEnergy, double& BpEnergy)
{

  phiEnergy = ApEnergy = BpEnergy = 0.;

  if(parallel->Coord[DIR_VMS] == 0) {
        
    [&](CComplex kXIn  [Nq][NzLD][Nky][FFTSolver::X_NkxL],
        CComplex kXOut [Nq][NzLD][Nky][FFTSolver::X_NkxL],
        CComplex Field0[Nq][NzLD][Nky][NxLD])
    {

      fft->solve(FFT_Type::X_FIELDS, FFT_Sign::Forward, ArrayField0.data((CComplex *) Field0));

      // Calculate |∇ φ|^2 (check if FFT Norm is correct)
      //
      //#pragma omp parallel for, collapse(2), reduce(+,phiEnergy,ApEnergy,BpEnergy)
      for(int z = NzLlD; z <= NzLuD; z++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { 
      
      for(int x_k = fft->K1xLlD; x_k <= fft->K1xLuD; x_k++) {
              
        const double k2_p = fft->k2_p(x_k,y_k,z);
        if(Nq >= 1) phiEnergy += k2_p * pow2(cabs(kXOut[Field::phi][z][y_k][x_k])) * (1. + plasma->debye2);
        if(Nq >= 2) ApEnergy  += cabs(k2_p * pow2(kXOut[Field::Ap ][z][y_k][x_k]))/fft->Norm_X;
        if(Nq >= 3) BpEnergy  += cabs(       pow2(kXOut[Field::Bp ][z][y_k][x_k]))/fft->Norm_X;
        //if(Nq >= 1) phiEnergy += (plasma->debye2*k2_p + sum_qqnT_1mG0(k2_p)) * pow2(cabs(kXOut[Field::phi][z][y_k][x_k]))/fft->Norm_X;
           
      } } }
      
      ////////////// add contribution to field energy arising from adiabatic term ( including flux-averaging correction ) ///////////////
      CComplex phi_yz[Nx]; phi_yz[:] = 0.;
      if(species[0].doGyro) calcFluxSurfAvrg(kXOut, phi_yz);
        
      const double adiab = species[0].n0 * pow2(species[0].q)/species[0].T0;
      // add contributions from adiabatic species (eventuellay with flux-averaging term)
      for(int z = NzLlD; z <= NzLuD; z++) { simd_for(int x_k = fft->K1xLlD; x_k <= fft->K1xLuD; x_k++) {
        
        phiEnergy += adiab * __sec_reduce_add(pow2(cabs((kXOut[Field::phi][z][:][x_k] - phi_yz[x_k])))) ;

      } } // z, x_k

      phiEnergy =  parallel->reduce(phiEnergy * grid->dXYZ / (8. * M_PI), Op::sum, DIR_XYZ) ; 
      ApEnergy  =  parallel->reduce( ApEnergy * grid->dXYZ / (8. * M_PI), Op::sum, DIR_XYZ) ;
      BpEnergy  =  parallel->reduce( BpEnergy * grid->dXYZ / (8. * M_PI), Op::sum, DIR_XYZ) ;
      
    } ((A4zz) fft->kXIn, (A4zz) fft->kXOut, (A4zz) Field0); 

  }     
}


///////////////////////////////////////// Some Helper functions ///////////////////////////////////////////



double FieldsFFT::sum_qqnT_1mG0(const double k2_p) 
{
   
  double g0 = 0.;

  for(int s = NsGlD; s <= NsGuD; s++) {

    const double qqnT   = species[s].n0 * pow2(species[s].q)/species[s].T0;
    const double rho_t2 = species[s].T0  * species[s].m / pow2(species[s].q * plasma->B0);
       
    g0 += qqnT * SpecialMath::_1mGamma0_Pade( rho_t2 * k2_p);
  //g0 += qqnT * SpecialMath::_1mGamma0( rho_t2 * k2_p);
        
  }
  return g0;

};
  
double FieldsFFT::sum_sa2qG0(const double kp_2) 
{
   
  double g0 = 0.;
    
  #pragma simd
  for(int s = NsGlD; s <= NsGuD; s++) {

    const double sa2q   = species[s].sigma * pow2(species[s].alpha) * species[s].q;
    const double rho_t2 = species[s].T0 * species[s].m / (pow2(species[s].q) * plasma->B0);

    const double b = rho_t2 * kp_2;
       
    g0  += sa2q * (1.e0 - SpecialMath::_1mGamma0_Pade(b));
   
  }
  return g0;

};


double FieldsFFT::sum_qnB_Delta(const double k2_p) 
{

  double g0 = 0.;
   
  for(int s = NsGlD; s <= NsGuD; s++) {

    const double qn    = species[s].q * species[s].n0;
    const double rho_t2 = species[s].T0 * species[s].m / (pow2(species[s].q) * plasma->B0);

    g0 += qn * SpecialMath::Delta(rho_t2 * k2_p);
  }
  return g0/plasma->B0;

};

double FieldsFFT::sum_2TnBB_Delta(const double k2_p) 
{
   
  double g0 = 0.;
   
  for(int s = NsGlD; s <= NsGuD; s++) {

    const double Tn     = species[s].T0 * species[s].n0;
    const double rho_t2 = species[s].T0 * species[s].m / (pow2(species[s].q) * plasma->B0);

    g0 += Tn * SpecialMath::Delta(rho_t2 * k2_p);
  }

  return 2. * g0 /pow2(plasma->B0) ;
};


