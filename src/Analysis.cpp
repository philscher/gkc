/*
 * =====================================================================================
 *
 *       Filename: Analysis.cpp
 *
 *    Description: Analysis functions for gkc
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include "Analysis.h"

extern "C" double cabs (CComplex z);
extern "C" double creal(CComplex z);
extern "C" double cimag(CComplex z);
extern "C" double carg (CComplex z);



enum SpecDir   {SPEC_NO=-1, SPEC_XY=0, SPEC_XZ=1, SPEC_YZ=2};


Analysis::Analysis(Parallel *_parallel, Vlasov *_vlasov, Fields *_fields, Grid *_grid, Setup *_setup, FFTSolver *_fft, FileIO *fileIO, Geometry *_geo) : 
  parallel(_parallel),setup(_setup), vlasov(_vlasov), grid(_grid), fields(_fields), geo(_geo),  fft(_fft),
     A4(FortranArray<4>()), A4_z(FortranArray<4>())

{

       // set initial energy
       initialEkin.resize(Range(0, NsGuD)); initialEkin = 0.e0;
       for(int s=NsLlD; s<= NsLuD; s++)  initialEkin(s) = getKineticEnergy((A6zz) vlasov->f0.dataZero(), V.dataZero(), M.dataZero(), s);
       initialEkin(0) = sum(initialEkin(RsLD));
       
       // Spectrum
//       if(setup->dirSpectrumAvrg & SPEC_XZ)  spectrumXZ.resize(fft->RkxL, fft->RkzL); spectrumXZ = 0.0;
//       if(setup->dirSpectrumAvrg & SPEC_YZ)  spectrumYZ.resize(fft->RkyL, fft->RkzL); spectrumYZ = 0.0;
//       if(setup->dirSpectrumAvrg & SPEC_XY)  spectrumXY.resize(fft->RkxL, fft->RkyL); spectrumXY = 0.0;

    
     A_xyz.resize(RxLD, RkyLD, RzLD); A_xyz = 0.;

     A4.resize(RxLD, RkyLD, RzLD, RsLD);
     A4_z.resize(RxLD, RkyLD, RzLD, RsLD);

       
     initDataOutput(setup, fileIO);
}



Analysis::~Analysis() {

   closeData();
    

}


void Analysis::getPowerSpectrum(CComplex  kXOut  [Nq][NzLD][NkyLD][FFTSolver::X_NkxL], 
                                CComplex  Field0 [Nq][NzLD][NkyLD][NxLD]      ,
                                double    pSpecX [plasma->nfields][Nx/2], double pSpecY [plasma->nfields][Nky],
                                double    pPhaseX[plasma->nfields][Nx/2], double pPhaseY[plasma->nfields][Nky])
{

  double   pSpec[plasma->nfields][Nx];
  CComplex pFreq[plasma->nfields][Nx];

  // Note :  take care that the FFT output is of for e.g. an 8 number sequence [0, 1, 2, 3, 4,-3,-2,-1]
  // Note : atan2 is not a linear function, thus phase has to be calculated at last
  if(parallel->Coord[DIR_VMS] == 0) {

    // Note : We have domain decomposition in X but not in Y
    fft->solve(FFT_X_FIELDS, FFT_FORWARD, fields->Field0.data());
         
    for(int n = 1; n <= plasma->nfields; n++) {
                
      // Mode power & Phase shifts for X (domain decomposed)
      omp_for(int x_k=fft->K1xLlD; x_k<= fft->K1xLuD; x_k++) {

        pSpec[n-1][x_k] = sqrt(__sec_reduce_add(cabs(kXOut[n][NzLlD:NzLD][NkyLlD:NkyLD][x_k]))/fft->Norm_X); 
        pFreq[n-1][x_k] =      __sec_reduce_add(     kXOut[n][NzLlD:NzLD][NkyLlD:NkyLD][x_k] )/fft->Norm_X;
        
      }
            
      
      // Mode power & Phase shifts for Y (not decomposed)
      omp_for(int y_k = NkyLlD; y_k <=  NkyLuD ; y_k++) { 
        
        pSpecY [n-1][y_k] = sqrt(__sec_reduce_add(cabs(Field0[n][NzLlD:NzLD][y_k][NxLlD:NxLD]))); 
        pPhaseY[n-1][y_k] = carg(__sec_reduce_add(     Field0[n][NzLlD:NzLD][y_k][NxLlD:NxLD]));
      }
      
    }
             
    // Finalized calculations for x-domain (decomposed and negative frequency modes)

    // Sum up with X-values from other processes
    parallel->collect(&pSpec[0][0], OP_SUM, DIR_XYZ, Nq * Nx);        
    parallel->collect(&pFreq[0][0], OP_SUM, DIR_XYZ, Nq * Nx);        
         
    // map back from fftw [k0, k1, k2, ..., k_Ny, -k_(N/y-2), ... -k_1]
    for(int x_k = 1; x_k < Nx/2; x_k++) pFreq [:][x_k] = pFreq[:][Nx-x_k];
    for(int x_k = 1; x_k < Nx/2; x_k++) pSpec [:][x_k] = pSpec[:][Nx-x_k];

    pPhaseX[:][0:Nx/2] = carg(pFreq[:][0:Nx/2]);
    pSpecX [:][0:Nx/2] =      pSpec[:][0:Nx/2];
         
  } // if DIR_XYZ

         return;
};


//////////////////////// Calculate scalar values ///////////////////////////



// note the hear flux rate have to be calculted after getkinetic energy
// this should go hand-in-hand with the temperature calculation
double Analysis::getTotalHeatFlux(int s) 
{
    double heatFlux = 0.e0;
  // for(int s = ((sp == TOTAL) ? NsLlD : sp); s <= NsLuD && ((sp != TOTAL) ? s == sp : true)  ; s++) heatFlux += sum(getHeatFluxKy(sp));
 //   if((s >= NsLlD) && (s <= NsLuD)) {
        heatFlux = sum(getHeatFluxKy(s));
 //   } 
    return parallel->collect(heatFlux);

};


double Analysis::getTotalParticleFlux(int s) 
{
      
    double particleFlux = 0.e0;
   
    //for(int s = ((sp == TOTAL) ? NsLlD : sp); s <= NsLuD && ((sp != TOTAL) ? s == sp : true)  ; s++) particleFlux += sum(getParticleFluxKy(sp)); 
    //for(int s = ((sp == TOTAL) ? NsLlD : sp); s <= NsLuD && ((sp != TOTAL) ? s == sp : true)  ; s++)

//    Array3R G(RFields, RsGD, RkyLD); G = 0.;
    
//    if((s >= NsLlD) && (s <= NsLuD)) {

        //G(RFields, RsLD, RkyLD) = getParticleFluxKy(s);
        particleFlux = sum(getParticleFluxKy(s));

//    }
    return parallel->collect(particleFlux);

};


double Analysis::getKineticEnergy(const CComplex f[NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB], 
                                  const double V[NvGB], const double M[NmGB], const int s) 
{
      
       double kineticEnergy=0.e0;
     
       const double v2_d6Z = M_PI * plasma->species(s).n0 * plasma->species(s).T0 * plasma->B0 * dv * dm * grid->dXYZ * plasma->species(s).scale_v ;
              
       #pragma omp parallel for reduction(+:kineticEnergy) collapse (3)
       for(int m=NmLlD; m<= NmLuD;m++) { for(int z=NzLlD; z<= NzLuD;z++) {  for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { 
              
              for(int v=NvLlD; v<=NvLuD; v++) { 
                 kineticEnergy += cabs(__sec_reduce_add(f[s][m][NzLlD:NzLD][y_k][NxLlD:NxLD][v]) * pow2(V[v]) * v2_d6Z);
              }
              kineticEnergy += cabs(__sec_reduce_add(f[s][m][NzLlD:NzLD][y_k][NxLlD:NxLD][NvLlD:NvLD]) * M[m]) * v2_d6Z;

       } } } 

       // substract initial kientic energy or not ?
      // return  (parallel->collect(kineticEnergy, OP_SUM, DIR_ALL) - initialEkin(sp))/((initialEkin(sp) == 0.) ? 1. : initialEkin(sp));
      return  parallel->collect(kineticEnergy, OP_SUM, DIR_ALL);

};

   


double Analysis::getEntropy(int sp) 
{
  return 0.; 
   double entropy = 0.e0;
   for(int s = NsLlD; s <= NsLuD ; s++) { 
                //entropy += grid->dXYZV * abs(pow2(sum(vlasov->f(RxLD, RkyLD, RzLD, RvLD, RmLD, s) 
                //                            - vlasov->f0 (RxLD, RkyLD, RzLD, RvLD, RmLD, s)))/ sum(vlasov->f0(RxLD, RkyLD, RzLD, RvLD, RmLD, s)));
                entropy += grid->dXYZV * abs(pow2(sum(vlasov->f(RxLD, RkyLD, RzLD, RvLD, RmLD, s) 
                                            - vlasov->f0 (RxLD, RkyLD, RzLD, RvLD, RmLD, s)))/ sum(vlasov->f0(RxLD, RkyLD, RzLD, RvLD, RmLD, s)));
   }
        
   return parallel->collect(entropy);
};



double Analysis::getParticleNumber(const CComplex f[NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB], const int s)
{ 
   double number = 0.e0;
        
   for(int m = NmLlD; m <= NmLuD; m++) { 
   
     const double pn_d6Z = M_PI  * dv * grid->dm(m) * grid->dXYZV;
   
     #pragma omp parallel for reduction(+:number)
     for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) {

               number +=  __sec_reduce_add(f[NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][y_k][NxLlD:NxLD][NvLlD:NvLD]) * pn_d6Z;

     }  }

     return parallel->collect(number, OP_SUM);
};


///////////////////////////////// Calculate Moments  /////////////////////////////////////////////////

void Analysis::getNumberDensity(const  CComplex f[NsLB][NmLB][NzLB][NkyLB][NxLB][NvLB],
                                       CComplex D[NsLD][NzLD][NkyLD][NxLD], 
                                const double V[NvGB], const double M[NmGB])
{

    for(int s = NsLlD; s <= NsLuD; s++) { for(int m = NmLlD; m <= NmLuD; m++) { 

       const double pn_d6Z = M_PI  * dv * grid->dm(m) * grid->dXYZV;

    for(int x=NxLlD; x<= NxLuD;x++) { for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int z=NzLlD; z<= NzLuD;z++){

          D[s][z][y_k][x] =  __sec_reduce_add(f[s][m][z][y_k][x][NvLlD:NvLD]) * pn_d6Z;

    }}} }} 

    parallel->collect(&D[NsLlD][NzLlD][NkyLlD][NxLlD], OP_SUM, NsLD * NzLD * NkyLD * NxLD, DIR_VM);

    return;
};


Array4R Analysis::getMomentumParallel() {
/* 
    const double d6Z = dXYZV * plasma->B0 * M_PI;
    const double alpha_s = plasma->species(s).q / plasma->species(s).T0;

       for(int v=NvLlD; v<= NvLuD;v++)  A4(x,y,z,s) = alpha_s * V(v) * F(x, y, z, v, RmLD, s) * d6Z;
 * */

    return parallel->collect(A4, OP_SUM, DIR_VM);
}


void Analysis::getTemperatureParallel(const  CComplex f[NsLB][NmLB][NzLB][NkyLB][NxLB][NvLB],
                                             CComplex A4_z[NsLD][NzLD][NkyLD][NxLD], 
                                      const double V[NvGB], const double M[NmGB])
{

  for(int s = NsLlD; s <= NsLuD; s++) {

      const double d6Z = M_PI * plasma->species(s).n0 * plasma->species(s).T0 * plasma->B0 * dv * dm * grid->dXYZ;

      omp_for_C3(int m=NmLlD; m<=NmLuD; m++) { for(int z=NzLlD; z<= NzLuD; z++) { for(int y_k=NkyLlD; y_k<= NkyLuD; y_k++) { 
            
        const Complex ky=Complex(0.,-fft->ky(y_k));
        
        for(int x=NxLlD; x<=NxLuD ; x++) { 

            // calculate temperature
            A4_z[s][z][y_k][x] = __sec_reduce_add((pow2(V[NvLlD:NvLD]) + M[m] * plasma->B0) * f[s][m][z][y_k][x][NvLlD:NvLD]) * d6Z;
         
        } 
      
      } } } // y_k, z, m
  
    } // s


    parallel->collect(&A4_z[NsLlD][NzLlD][NkyLlD][NxLlD], OP_SUM, NsLD * NzLD * NkyLD * NxLD, DIR_VM);

    return;
};



Array4R Analysis::getTemperatureOthogonal() {
/* 
    const double d6Z = dXYZV * plasma->B0 * M_PI;
       for(int m=NmLlD; m<= NmLuD;m++)  A4(x,y,z,s) = (M(m) * plasma->B0 - 1.) * F(x, y, z, RvLD, m, s) * d6Z;
 * */
    return parallel->collect(A4, OP_SUM, DIR_VM);
}


//////////////////////// Calculate additional heat fluxes /////////////////////

Array4R Analysis::getHeatFluxParallel() {
/* 
    for(int s = NsLlD; s <= NsLuD; s++) {

    const double d6Z = dXYZV * plasma->B0 * M_PI;
    const double alpha_s = plasma->species(s).q / plasma->species(s).T0;

    for(int x=NxLlD; x<= NxLuD;x++) { for(int y=NyLlD; y<= NyLuD;y++) { for(int z=NzLlD; z<= NzLuD;z++){


       for(int v=NvLlD; v<= NvLuD;v++)  A4(x,y,z,s) = alpha_s * V(v) * (pow2(V(v)) - 3./2.) * abs(F(x, y, z, v, RmLD, s)) * d6Z;

    }}} }

 * */
    return parallel->collect(A4, OP_SUM, DIR_VM);
}

Array4R Analysis::getHeatFluxOrthogonal() {
/* 
    for(int s = NsLlD; s <= NsLuD; s++) {

    const double alpha_s = plasma->species(s).q / plasma->species(s).T0;
    const double d6Z = dXYZV * plasma->B0 * M_PI;

    for(int x=NxLlD; x<= NxLuD;x++) { for(int y=NyLlD; y<= NyLuD;y++) { for(int z=NzLlD; z<= NzLuD;z++){
    for(int v=NvLlD; v<= NvLuD;v++) { for(int m=NmLlD; m<= NmLuD;m++)  {


         A4(x,y,z,s) = alpha_s * V(v) * (M(m) * plasma->B0 - 1.) * F(x, y, z, v, RmLD, s) * d6Z;

    }} }}} }

 * */
    return parallel->collect(A4, OP_SUM, DIR_VM);
}


// note the hear flux rate have to be calculted after getkinetic energy
// this should go hand-in-hand with the temperature calculation
Array4C Analysis::getHeatFlux(int sp) 
{
    A4_z = 0.;

    Array3C V2(RxLD, RkyLD, RzLD); 
    Array3C W2(RxLD, RkyLD, RzLD); 
   
   //for(int s = ((sp == TOTAL) ? NsLlD : sp); s <= NsLuD && ((sp != TOTAL) ? s == sp : true)  ; s++) {
  // for(int s = ((sp == TOTAL) ? NsLlD : sp); (s <= ((sp == TOTAL) ? NsLuD : sp)) && ( (sp >= NsLlD) && (sp <= NsLuD))  ; s++) { 
  for(int s = NsLlD; s <= NsLuD; s++) {
      
      const double d6Z = M_PI * plasma->species(s).n0 * plasma->species(s).T0 * plasma->B0 * dv * dm * grid->dXYZ;

      for(int m=NmLlD; m<= NmLuD;m++ ) { 
      
        V2 = 0.; W2 = 0.;
        // multiply in real-space 
      for(int z=NzLlD; z<= NzLuD;z++){ for(int x=NxLlD; x  <= NxLuD ;  x++) { for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { 

            const Complex ky=Complex(0.,-fft->ky(y_k));
              
            V2(x,y_k,z) = ky * fields->phi(x,y_k,z,m,s);
        
             // integrate over velocity space
             for(int v=NvLlD; v<= NvLuD;v++) W2(x,y_k,z) += (pow2(V(v)) +  M(m) * plasma->B0) * vlasov->f(x, y_k, z, v, m, s) * d6Z;

         }}}

            A4_z(RxLD, RkyLD, RzLD,s) = fft->multiply(V2, W2, A_xyz);
     
      } }

    return parallel->collect(A4_z, OP_SUM, DIR_VM);


};



// note the hear flux rate have to be calculted after getkinetic energy
// this should go hand-in-hand with the temperature calculation
Array3R Analysis::getHeatFluxKy(int sp) 
{
    Array3R Q(RFields, RkyLD, RsLD); Q = 0.;
      
    getHeatFlux(sp);
    //
    for(int s = NsLlD; s <= NsLuD; s++) {
        
       for(int z=NzLlD; z<= NzLuD;z++) { for(int x=NxLlD; x  <= NxLuD ;  x++) { Q(Field::phi, RkyLD, s) +=  abs(A4_z(x, RkyLD,z,s)) ; }}
                 
    }
   return Q; 

};





Array3R  Analysis::getParticleFluxKy(int sp) 
{
    Array3R G(RFields, RkyLD, RsLD); G = 0.;
    
    Array3C V2(RxLD, RkyLD, RzLD); V2 = 0.;
    Array3C W2(RxLD, RkyLD, RzLD); W2 = 0.;
   
  for(int s = NsLlD; s <= NsLuD; s++) {
      
      const double d6Z = M_PI * plasma->species(s).n0 * plasma->species(s).T0 * plasma->B0 * dv * dm * grid->dXYZ;

        for(int m=NmLlD; m<= NmLuD;m++ ) { for(int z=NzLlD; z<= NzLuD;z++){ 
      
        // multiply in real-space 
        for(int x=NxLlD; x  <= NxLuD ;  x++) { for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { 

            const Complex ky=Complex(0.,-fft->ky(y_k));
              
            V2(x,y_k,z) = ky * fields->phi(x,y_k,z,m,s);
        
             // integrate over velocity space
             W2(x,y_k,z) = sum(vlasov->f(x, y_k, z, RvLD, m, s)) * d6Z;
         }} 

           // V2(RxLD, RkyLD, RzLD) = 
           fft->multiply(V2, W2, A_xyz);
      
            // sum over x and z
            for(int x=NxLlD; x  <= NxLuD ;  x++) { G(Field::phi, RkyLD, s) +=  abs(A_xyz(x, RkyLD,z)) ; }
                 
        } }
   }


    
   return parallel->collect(G);

};

      


 
        
int Analysis::updateSpectrum(unsigned int dir) {
       /* 
        firstIndex  i_idx;
        secondIndex j_idx;
        thirdIndex  k_idx;
       if     (dir == SPEC_XZ)  {
            double avrg_factor = 1.e0 / ((double) (setup->spectrumAvrg[SPEC_XZ][SPEC_END] - setup->spectrumAvrg[SPEC_XZ][SPEC_START]));
            spectrumXZ = spectrumXZ + avrg_factor/Ny * sum(phi_k, j_idx);
      }
       else if(dir == SPEC_YZ)  {
            double avrg_factor = 1.e0 / ((double) (setup->spectrumAvrg[SPEC_YZ][SPEC_END] - setup->spectrumAvrg[SPEC_YZ][SPEC_START]));
//spectrumYZ = spectrumYZ + avrg_factor/Nx * sum(phi_k, i_idx);
            spectrumYZ = spectrumYZ + avrg_factor * phi_k(5, fft->RkyL, fft->RkzL);
      }
       else if(dir == SPEC_XY)  {
            double avrg_factor = 1.e0 / ((double) (setup->spectrumAvrg[SPEC_XY][SPEC_END] - setup->spectrumAvrg[SPEC_XY][SPEC_START]));
            spectrumXY = spectrumXY + avrg_factor/Nz * sum(phi_k, k_idx);
      }
      else check(-1, DMESG("No such direction for spectrum")); 
  

      // set to zero
        * */ 
        return GKC_SUCCESS;
   };

////////////////////////     Calculate x-dependent values     /////////////////////

Array2C Analysis::getSpectrum(unsigned int dir) {
       Array2C spectrum;
       if     (dir == SPEC_XZ)  spectrum.reference(spectrumXZ);
       else if(dir == SPEC_YZ)  spectrum.reference(spectrumXZ);
       else if(dir == SPEC_XY)  spectrum.reference(spectrumXY);
       else   check(-1, DMESG("No such direction for get spectrum"));
    
       return spectrum;
   };

  

void Analysis::initDataOutput(Setup *setup, FileIO *fileIO) {
        
     analysisGroup = fileIO->newGroup(fileIO->getFileID(), "/Analysis");
        
     hsize_t offset0[] = { 0, 0, 0, 0, 0, 0, 0 };
    
     //###################################### Analysis - Heat fluxes ################################
     
     // Heat Flux ky and Particle FluxKy ( per species) 
     hid_t fluxGroup = fileIO->newGroup(analysisGroup, "Flux");
     
     hsize_t FSky_dim[]       = { plasma->nfields, Nky, Ns  , 1 }; 
     hsize_t FSky_maxdim[]    = { plasma->nfields, Nky, Ns  , H5S_UNLIMITED} ;
     hsize_t FSky_chunkdim[]  = { plasma->nfields, Nky, NsLD, 1 };
     hsize_t FSky_chunkBdim[] = { plasma->nfields, Nky, NsLD, 1 };
     hsize_t FSky_offset[]    = { 0, 0  , NsLlD-1, 0  };

     FA_heatKy      = new FileAttr("Heat"   , fluxGroup, 4, FSky_dim, FSky_maxdim, FSky_chunkdim, offset0,  FSky_chunkBdim, FSky_offset, parallel->Coord[DIR_XYZVM] == 0);
     FA_particleKy  = new FileAttr("Density", fluxGroup, 4, FSky_dim, FSky_maxdim, FSky_chunkdim, offset0,  FSky_chunkBdim, FSky_offset, parallel->Coord[DIR_XYZVM] == 0);
    
     H5Gclose(fluxGroup);
     
     //###################################### Moments - Heat fluxes ################################
     hid_t momentGroup = check(H5Gcreate(fileIO->getFileID(), "/Moments",H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), DMESG("Error creating group file for Phi : H5Gcreate"));

     hsize_t moment_dim[]       =  { grid->NzGD , NkyLD      , grid->NxGD , NsLD , 1};
     hsize_t moment_maxdim[]    =  { grid->NzGD , NkyLD      , grid->NxGD , Ns   , H5S_UNLIMITED};
     hsize_t moment_chunkBdim[] =  { NzLD       , NkyLD      , NxLD       , NsLD , 1          };
     hsize_t moment_chunkdim[]  =  { NzLD       , NkyLD      , NxLD       , NsLD , 1};
     hsize_t moment_moffset[]   =  { 0, 0, 0, 0, 0 };
     hsize_t moment_offset[]    =  { NzLlD-3, 0     , NxLlD-3, 0     ,  0  };
     
     bool momWrite = (parallel->Coord[DIR_VM] == 0);
     
     
     FA_Mom_Tp        = new FileAttr("Temperature_v", momentGroup, 5, moment_dim , moment_maxdim   , moment_chunkdim   , moment_moffset    ,  moment_chunkBdim  , moment_offset, momWrite, fileIO->complex_tid);
     FA_Mom_HeatFlux  = new FileAttr("HeatFlux"     , momentGroup, 5, moment_dim , moment_maxdim   , moment_chunkdim   , moment_moffset    ,  moment_chunkBdim  , moment_offset, momWrite, fileIO->complex_tid);
     FA_Mom_Density   = new FileAttr("Density"      , momentGroup, 5, moment_dim , moment_maxdim   , moment_chunkdim   , moment_moffset    ,  moment_chunkBdim  , moment_offset, momWrite, fileIO->complex_tid);
     FA_Mom_Time  = fileIO->newTiming(momentGroup);
        
     H5Gclose(momentGroup);
      
     dataOutputMoments   = Timing(setup->get("DataOutput.Moments.Step", -1)       , setup->get("DataOutput.Moments.Time", -1.));


     //###################################### Power Spectrum  ################################################
     // X-scalarValue
      hid_t growGroup = fileIO->newGroup(analysisGroup, "PowerSpectrum");

     hsize_t grow_x_dim[]       = { plasma->nfields, Nx/2+1, 1 }; 
     hsize_t grow_x_maxdim[]    = { plasma->nfields, Nx/2+1, H5S_UNLIMITED} ;
     hsize_t grow_x_chunkdim[]  = { plasma->nfields, Nx/2+1, 1 };
     hsize_t grow_x_chunkBdim[] = { plasma->nfields, Nx/2+1, 1 };
     FA_grow_x  = new FileAttr("X", growGroup, 3, grow_x_dim, grow_x_maxdim, grow_x_chunkdim, offset0,  grow_x_chunkBdim, offset0, parallel->myRank == 0);

     // Y-scalarValue
     hsize_t grow_y_dim[]       = { plasma->nfields, Nky, 1 };
     hsize_t grow_y_maxdim[]    = { plasma->nfields, Nky, H5S_UNLIMITED };
     hsize_t grow_y_chunkdim[]  = { plasma->nfields, Nky, 1 };
     hsize_t grow_y_chunkBdim[] = { plasma->nfields, Nky, 1 };
   
     FA_grow_y  = new FileAttr("Y", growGroup, 3, grow_y_dim, grow_y_maxdim, grow_y_chunkdim, offset0,  grow_y_chunkBdim, offset0, parallel->myRank == 0);

     FA_grow_t  = fileIO->newTiming(growGroup);

     H5Gclose(growGroup);
     
     
     hid_t freqGroup = fileIO->newGroup(analysisGroup, "PhaseShift");
     FA_freq_x  = new FileAttr("X", freqGroup, 3, grow_x_dim, grow_x_maxdim, grow_x_chunkdim, offset0,  grow_x_chunkBdim, offset0, parallel->myRank == 0);
     FA_freq_y  = new FileAttr("Y", growGroup, 3, grow_y_dim, grow_y_maxdim, grow_y_chunkdim, offset0,  grow_y_chunkBdim, offset0, parallel->myRank == 0);
     FA_freq_t  = fileIO->newTiming(freqGroup);
     H5Gclose(freqGroup);

      //////////////////////////////////////////////////////////////// Setup Table for scalar data ////////////////////////////////////////////////////////
              
      ScalarValues_t scalarValues;
     
      size_t SV_offset[] = { HOFFSET( ScalarValues_t, timestep       ), HOFFSET( ScalarValues_t, time     ), HOFFSET( ScalarValues_t, phiEnergy       ),
                             HOFFSET( ScalarValues_t, ApEnergy       ), HOFFSET( ScalarValues_t, BpEnergy ), HOFFSET( ScalarValues_t, particle_number ),
                             HOFFSET( ScalarValues_t, kinetic_energy ), HOFFSET( ScalarValues_t, entropy  ), HOFFSET( ScalarValues_t, heat_flux       ),
                             HOFFSET( ScalarValues_t, particle_flux  ) };

      size_t SV_sizes[] = { sizeof(scalarValues.timestep), sizeof(scalarValues.time    ), sizeof(scalarValues.phiEnergy), 
                            sizeof(scalarValues.ApEnergy), sizeof(scalarValues.BpEnergy), Ns * sizeof(scalarValues.particle_number[0]), Ns * sizeof(scalarValues.kinetic_energy[0]), 
                           Ns * sizeof(scalarValues.entropy[0]), Ns * sizeof(scalarValues.heat_flux[0]), Ns * sizeof(scalarValues.particle_flux[0])};

      hid_t SV_types[] = { H5T_NATIVE_INT, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, 
                           fileIO->species_tid, fileIO->species_tid, fileIO->species_tid, fileIO->species_tid, fileIO->species_tid } ;
  
      const char *SV_names[] = { "Timestep", "Time", "phiEnergy", "ApEnergy", "BpEnergy", "ParticleNumber", "KineticEnergy", "Entropy", "HeatFlux", "ParticleFlux" };

      SVTable = new TableAttr(analysisGroup, "scalarValues", 10, SV_names, SV_offset, SV_types, SV_sizes, &scalarValues); 

      dataOutputStatistics  = Timing(setup->get("DataOutput.Statistics.Step", -1), setup->get("DataOutput.Statistics.Time", -1.));

}
 


void Analysis::getFieldEnergy(double& phiEnergy, double& ApEnergy, double& BpEnergy)
{

   fields->getFieldEnergy(phiEnergy, ApEnergy, BpEnergy);

};
  
int Analysis::writeData(Timing timing, double dt)

{

  // important - include X-dependence values

  if (timing.check(dataOutputMoments, dt)       )   {

    getTemperatureParallel((A6zz ) vlasov->f.dataZero(), (A4zz) A4_z.dataZero(), V.dataZero(), M.dataZero());
//           FA_Mom_Tp->write(A4_z.data());
//           FA_Mom_HeatFlux->write(getHeatFlux().data());
//           FA_Mom_Density->write(getNumberDensity().data());
//           FA_Mom_Density->write(A4_z.data());
//           FA_Mom_Time->write(&timing);
    writeMessage("Data I/O : Moments output");
    
  }

  if (timing.check(dataOutputStatistics, dt)       )   {
  
    // Stack allocation (size ok ?) and alignement ? (well, speed not critical here...)
    double pSpecX [plasma->nfields] [Nx/2], pSpecY [plasma->nfields] [Nky ],
           pPhaseX[plasma->nfields] [Nx/2], pPhaseY[plasma->nfields] [Nky ];
     
    getPowerSpectrum((A4zz) fft->kXOut.dataZero(), (A4zz) fields->Field0.dataZero(), pSpecX, pSpecY, pPhaseX, pPhaseY);
    
    // Seperatly writing ? Hopefully it is buffered ... (passing stack pointer ... OK ?)
    FA_grow_x->write( &pSpecX[0][0] ); FA_grow_y->write(&pSpecY[0][0]); FA_grow_t->write(&timing);
    FA_freq_x->write( &pPhaseX[0][0]); FA_freq_y->write(&pPhaseY[0][0]); FA_freq_t->write(&timing);
    
    // Heat Flux skip this crap
    //Array3R heatKy; heatKy.reference(getHeatFluxKy());
    //FA_heatKy->write(heatKy.data());
    //Array3R particleKy; particleKy.reference(getParticleFluxKy());
    //FA_particleKy->write(particleKy.data());
    
    ScalarValues scalarValues;
    
    // calculate kinetic Energy first, need for initial_e ! sum over domain
    scalarValues.timestep = timing.step;
    scalarValues.time     = timing.time;
    
    fields->getFieldEnergy(scalarValues.phiEnergy, scalarValues.ApEnergy, scalarValues.BpEnergy);
    
    //  Get scalar Values for every species ( this is bad calculate them alltogether)
    for(int s = NsGlD; s <= NsGuD; s++) {
    
      scalarValues.particle_number[s-1]  = getParticleNumber((A6zz) vlasov->f.dataZero(), s);
      scalarValues.entropy        [s-1]  = getEntropy(s)                               ;
      scalarValues.kinetic_energy [s-1]  = getKineticEnergy((A6zz) vlasov->f.dataZero(), V.dataZero(), M.dataZero(), s) ;
      scalarValues.particle_flux  [s-1]  = getTotalParticleFlux(s)                     ;
      scalarValues.heat_flux      [s-1]  = getTotalHeatFlux(s)                         ;
    }
      
    
    SVTable->append(&scalarValues);
    
    // write out to Terminal/File
    
    std::stringstream messageStream;
    messageStream << std::endl << std::endl << "Analysis | " << std::setprecision(3);
    messageStream << "Field Energy : (phi) " << scalarValues.phiEnergy  << "  (Ap) " << scalarValues.ApEnergy  <<  "  (Bp) " << scalarValues.BpEnergy << std::endl; 
    double charge = 0., kinetic_energy=0.;
    
    for(int s = NsGlD; s <= NsGuD; s++) {
    
      messageStream << "         | " << 
        plasma->species(s).name << "   N : " << scalarValues.particle_number[s-1]  << "  Kinetic Energy: " << scalarValues.kinetic_energy[s-1] ;
      messageStream << "   Particle Flux :" << scalarValues.particle_flux[s-1]    << "  Heat Flux : " << scalarValues.heat_flux[s-1] << std::endl;
        charge += plasma->species(s).q  * scalarValues.particle_number[s-1];
      kinetic_energy += scalarValues.kinetic_energy[s-1];
    }
    
    messageStream << //"------------------------------------------------------------------" <<
      "         | Total Energy " << kinetic_energy+scalarValues.phiEnergy + scalarValues.ApEnergy + scalarValues.BpEnergy << "    Total Charge = " << ((plasma->species(0).n0 != 0.) ? 0. : charge) 
      << std::endl;  
    parallel->print(messageStream.str());
    
  }
  
  return GKC_SUCCESS;
  
}

     //###################################### Spectrum  ################################################
     // YZ-Spectrum
/*
     hsize_t spec_yz_dim[]      = { Ny/2+1, Nz/2+1, 1 };
     hsize_t spec_yz_maxdim[]   = { Ny/2+1, Nz/2+1, H5S_UNLIMITED};
     hsize_t spec_yz_chunkdim[] = { kNyL, kNzL, 1};

     FA_spec_yz  = new FileAttr("Unnamed",3, spec_yz_dim, spec_yz_maxdim, spec_yz_chunkdim, offset0,  spec_yz_chunkdim, parallel->isFFTGroup);
     
     // XY-Spectrum
     hsize_t spec_xy_dim[]      = { Nx/2+1, Ny/2+1, 1};
     hsize_t spec_xy_maxdim[]   = { Nx/2+1, Ny/2+1, H5S_UNLIMITED};
     hsize_t spec_xy_chunkdim[] = { kNxL, kNyL, 1};

     FA_spec_xy  = new FileAttr("Unnamed",3, spec_xy_dim, spec_xy_maxdim, spec_xy_chunkdim, offset0,  spec_xy_chunkdim, parallel->isFFTGroup);
     
     // XZ-Spectrum
     hsize_t spec_xz_dim[]      = { Nx/2+1, Nz/2+1, 1 };
     hsize_t spec_xz_maxdim[]   = { Nx/2+1, Nz/2+1, H5S_UNLIMITED };
     hsize_t spec_xz_chunkdim[] = { kNyL, kNzL, 1 };
     
     FA_spec_xz  = new FileAttr("Unnamed",3, spec_xz_dim, spec_xz_maxdim, spec_xz_chunkdim, offset0,  spec_xz_chunkdim, parallel->isFFTGroup);
*/

        
/* 
  int FileIO::writeSpectrum(Array3C phi_k, Analysis *analysis, Timing timing)
  {
      return GKC_SUCCESS;
    if((parallel->isFFTGroup) && (setup->spectrumAvrg[SPEC_XZ][SPEC_START] <= timing.step) && (timing.step <= setup->spectrumAvrg[SPEC_XZ][SPEC_END])) 
      analysis->updateSpectrum(SPEC_XZ);
    
    if((parallel->isFFTGroup) && (setup->spectrumAvrg[SPEC_YZ][SPEC_START] <= timing.step) && (timing.step <= setup->spectrumAvrg[SPEC_YZ][SPEC_END])) 
      analysis->updateSpectrum(SPEC_YZ);


    if((parallel->isFFTGroup) && (setup->spectrumAvrg[SPEC_XY][SPEC_START] <= timing.step) && (timing.step <= setup->spectrumAvrg[SPEC_XY][SPEC_END])) 
      analysis->updateSpectrum(SPEC_XY);


//      if(timeStep == setup->spectrumAvrg[SPEC_YZ][SPEC_END]) FA_spec_yz->write(analysis->getSpectrum(SPEC_YZ).data());
//      if(timeStep == setup->spectrumAvrg[SPEC_XZ][SPEC_END]) FA_spec_xz->write(analysis->getSpectrum(SPEC_XZ).data());
//      if(timeStep == setup->spectrumAvrg[SPEC_XY][SPEC_END]) FA_spec_xy->write(analysis->getSpectrum(SPEC_XY).data());

      return GKC_SUCCESS;
        

  }
 */


// separately we set MPI struct
void Analysis::setMPIStruct()
{
/*  
   long long int 

   MPI_Aint mpi_addr_SV;
   MPI_Aint mpi_addr_SV_timestep;

   addr_struct = (long long int) &_ScalarValues;
   addr_struct = (long long int) &_ScalarValues.timestep;

   MPI_Get_address(&_ScalarValues         , &mpi_addr_SV         );
   MPI_Get_address(&_ScalarValues.timestep, &mpi_addr_SV_timestep);

   types[] = { MPI_DOUBLE };

   blocklengths[] = { 1 };

   displacements[] = { mpi_addr_SV - mpi_addr_SV };

   MPI_Type_create_struct(1, blocklengths, displacements, types, &mpi_SV_t);
*/

};



void Analysis::closeData() {
  delete FA_heatKy; 
  delete FA_particleKy;
  delete FA_grow_x; delete FA_grow_y; delete FA_grow_t;     
  delete FA_freq_x; delete FA_freq_y; delete FA_freq_t;    
       
  delete FA_Mom_Tp;
  delete FA_Mom_HeatFlux;
  delete FA_Mom_Density;
  delete FA_Mom_Time;

  delete SVTable;
  H5Gclose(analysisGroup);


}
    
void Analysis::printOn(ostream &output)  const
{ 



};
