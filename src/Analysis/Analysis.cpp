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



Analysis::Analysis(Parallel *_parallel, Vlasov *_vlasov, Fields *_fields, Grid *_grid, Setup *_setup, FFTSolver *_fft, FileIO *fileIO, Geometry *_geo) : 
 
parallel(_parallel),setup(_setup), vlasov(_vlasov), grid(_grid), fields(_fields), geo(_geo),  fft(_fft)

{

     initData(setup, fileIO);
}



Analysis::~Analysis() 
{

   closeData();
    

}


void Analysis::getPowerSpectrum(CComplex  kXOut  [Nq][NzLD][NkyLD][FFTSolver::X_NkxL], 
                                CComplex  Field0 [Nq][NzLD][NkyLD][NxLD]      ,
                                double    pSpecX [Nq][Nx/2+1], double pSpecY [Nq][Nky],
                                double    pPhaseX[Nq][Nx/2+1], double pPhaseY[Nq][Nky])
{

  double   pSpec[Nq][Nx]; pSpec[:][:] = 0.;
  CComplex pFreq[Nq][Nx]; pFreq[:][:] = 0.;

  // Note :  take care that the FFT output is of for e.g. an 8 number sequence [0, 1, 2, 3, 4,-3,-2,-1]
  // Note : atan2 is not a linear function, thus phase has to be calculated at last
  if(parallel->Coord[DIR_VMS] == 0) {

    // Note : We have domain decomposition in X but not in Y
    fft->solve(FFT_Type::X_FIELDS, FFT_Sign::Forward, &Field0[1][NzLlD][NkyLlD][NxLlD]);
         
    for(int q = 1; q <= Nq; q++) {
                
      // Mode power & Phase shifts for X (domain decomposed) |\phi(k_x)|
      omp_for(int x_k=fft->K1xLlD; x_k <= fft->K1xLuD; x_k++) {

        pSpec[q-1][x_k] = __sec_reduce_add(cabs(kXOut[q][NzLlD:NzLD][NkyLlD:NkyLD][x_k]))/fft->Norm_X; 
        pFreq[q-1][x_k] = __sec_reduce_add(     kXOut[q][NzLlD:NzLD][NkyLlD:NkyLD][x_k] )/fft->Norm_X;
      }
            
      
      // Mode power & Phase shifts for Y (not decomposed) |\phi(k_y)|
      omp_for(int y_k = NkyLlD; y_k <=  NkyLuD ; y_k++) { 
        
        pSpecY [q-1][y_k] =      __sec_reduce_add(cabs(Field0[q][NzLlD:NzLD][y_k][NxLlD:NxLD])); 
        pPhaseY[q-1][y_k] = carg(__sec_reduce_add(     Field0[q][NzLlD:NzLD][y_k][NxLlD:NxLD]));
      }
      
    }
             
    // Finalized calculations for x-domain (decomposed and negative frequency modes)

    // Sum up with X-values from other processes
    parallel->reduce(&pSpec[0][0], Op::SUM, DIR_XYZ, Nq * Nx);        
    parallel->reduce(&pFreq[0][0], Op::SUM, DIR_XYZ, Nq * Nx);        
    
    // Sum up with Y-values from other processes
    parallel->reduce(&pSpecY [0][0], Op::SUM, DIR_XYZ, Nq * Nky);        
    parallel->reduce(&pPhaseY[0][0], Op::SUM, DIR_XYZ, Nq * Nky);        
         
         
    // map back from fftw [k0, k1, k2, ..., k_Ny, -k_(N/y-2), ... -k_1]
    for(int x_k = 1; x_k < Nx/2; x_k++) pFreq [:][x_k] += pFreq[:][Nx-x_k];
    for(int x_k = 1; x_k < Nx/2; x_k++) pSpec [:][x_k] += pSpec[:][Nx-x_k];

    pPhaseX[:][0:Nx/2+1] = carg(pFreq[:][0:Nx/2+1]);
    pSpecX [:][0:Nx/2+1] =      pSpec[:][0:Nx/2+1];
         
  } // if DIR_XYZ

  return;
};


//////////////////////// Calculate scalar values ///////////////////////////


void Analysis::calculateScalarValues(const CComplex f [NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB], 
                                     const CComplex f0[NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB], 
                                     const double V[NvGB], const double M[NmGB], const int s,
                                     ScalarValues &scalarValues) 

{
  
    for(int s = NsLlD; s <= NsLuD; s++) {
      
    double number = 0.e0;
    double kineticEnergy=0.e0;
    double entropy = 0.;
    double particle = 0.; 
    double heat = 0.;
    
    for(int m=NmLlD; m<= NmLuD;m++) {
      
    ////////////////////////////// Calculate Particle Number ////////////////////////
    const double pn_d6Z = M_PI  * dv * grid->dm[m] * grid->dXYZV;
   
    //#pragma omp parallel for reduction(+:number)
    for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) 
    {
               number +=  __sec_reduce_add(f[s][m][NzLlD:NzLD][y_k][NxLlD:NxLD][NvLlD:NvLD]) * pn_d6Z;
    } 

    //////////// Calculate Kinetic Energy  //////////////////////////
    
    const double v2_d6Z = 0.5 * plasma->species[s].m * M_PI * plasma->species[s].n0 * dv * grid->dm[m] * grid->dXYZ * pow2(plasma->species[s].scale_v) ;
              
    //#pragma omp parallel for reduction(+:kineticEnergy) collapse (2)
    for(int z=NzLlD; z<= NzLuD;z++) {  for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int x=NxLlD; x <= NxLuD; x++) {
             
                 kineticEnergy += cabs(__sec_reduce_add(f[s][m][z][y_k][x][NvLlD:NvLD] * (pow2(V[NvLlD:NvLD])+M[m]) ) * v2_d6Z);
    } } }
    // substract initial kinetic energy or not ?
    // return  (parallel->reduce(kineticEnergy, OP_SUM, DIR_ALL) - initialEkin(sp))/((initialEkin(sp) == 0.) ? 1. : initialEkin(sp));

    ////////////////////////////// Calculate Entropy ////////////////////////////////////
    
    //#pragma omp parallel for reduction(+:kineticEnergy) collapse (2)
    for(int z=NzLlD; z<= NzLuD;z++) {  for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) {
    for(int x=NxLlD; x<= NxLuD;x++) { 

      entropy += cabs(pow2(__sec_reduce_add(f [s][m][z][y_k][x][NvLlD:NvLD])))/
                           __sec_reduce_add(f0[s][m][z][y_k][x][NvLlD:NvLD]);
    } } }
    /////////////////////////// Calculate Total Heat & Particle Flux ////////////////////////
    CComplex ParticleFlux[NkyLD][NxLD], HeatFlux[NkyLD][NxLD];
    ParticleFlux[:][:] = 0.;
    HeatFlux[:][:] = 0.;

    getParticleHeatFlux(m, s, ParticleFlux, HeatFlux, f, (A6zz) fields->Field, V, M);

    particle += __sec_reduce_add(cabs(ParticleFlux[:][:]));
    heat     += __sec_reduce_add(cabs(    HeatFlux[:][:]));
    // so bad .... (looks to ugly, to be the right way ...)
    } // m 

    // Communicate with other groups (make reduce whole structre)
    heat          = parallel->reduce(heat         , Op::SUM, DIR_M);
    particle      = parallel->reduce(particle     , Op::SUM, DIR_M);
    kineticEnergy = parallel->reduce(kineticEnergy, Op::SUM, DIR_M);
    number        = parallel->reduce(number       , Op::SUM, DIR_M);

    #pragma omp single
    {
       scalarValues.particle_number[s-1]  = number        ;
       scalarValues.entropy        [s-1]  = entropy        ;
       scalarValues.kinetic_energy [s-1]  = kineticEnergy;
       scalarValues.particle_flux  [s-1]  = particle;
       scalarValues.heat_flux      [s-1]  = heat ;
    }  

    } // s

    return;
};


///////////////////////////////// Calculate Moments  /////////////////////////////////////////////////

void Analysis::getNumberDensity(const  CComplex f[NsLB][NmLB][NzLB][NkyLD][NxLB][NvLB],
                                       CComplex R[NsLD][NzLD][NkyLD][NxLD], 
                                const double V[NvGB], const double M[NmGB])
{
  R[:][:][:][:] = 0.;

  for(int s=NsLlD; s<= NsLuD; s++) { for(int m = NmLlD; m <= NmLuD; m++) { 

  const double pn_d6Z = M_PI  * dv * grid->dm[m] * grid->dXYZV;

  for(int z=NzLlD; z<= NzLuD;z++) { for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { 
  for(int x=NxLlD; x<= NxLuD;x++) { 

          R[s-NsLlD][z-NzLlD][y_k][x-NxLlD] +=  __sec_reduce_add(f[s][m][z][y_k][x][NvLlD:NvLD]) * pn_d6Z;

  } } }
  
  } } 

  parallel->reduce( (CComplex *) R, Op::SUM, DIR_VM, NsLD * NzLD * NkyLD * NxLD);

  return;

};


void Analysis::getTemperature(const  CComplex f[NsLB][NmLB][NzLB][NkyLD][NxLB][NvLB],
                                     CComplex R[NsLD][NzLD][NkyLD][NxLD], 
                              const double V[NvGB], const double M[NmGB])
{
  // We cannot rely that R is set to zero (and we use += later)
  R[:][:][:][:] = 0.;


  for(int s = NsLlD; s <= NsLuD; s++) { omp_for(int m=NmLlD; m<=NmLuD; m++) { 

     const double d6Z = M_PI * plasma->species[s].n0 * plasma->species[s].T0 * plasma->B0 * dv * grid->dm[m] * grid->dXYZ;

     omp_C2_for(int z=NzLlD; z <= NzLuD; z++) { for(int y_k=NkyLlD; y_k<= NkyLuD; y_k++) { 
            
        for(int x=NxLlD; x<=NxLuD ; x++) { 

            // calculate temperature
            R[s-NsLlD][z-NzLlD][y_k][x-NxLlD] += __sec_reduce_add((pow2(V[NvLlD:NvLD]) + M[m] * plasma->B0) * f[s][m][z][y_k][x][NvLlD:NvLD]) * d6Z;
         
        } 
      
     } } // y_k, z
  
    } } // m, s


    parallel->reduce((CComplex *) R, Op::SUM, NsLD * NzLD * NkyLD * NxLD, DIR_VM);

    return;
};





void Analysis::getParticleHeatFlux(const int m, const int s, 
                                   CComplex ParticleFlux[NkyLD][NxLD], CComplex HeatFlux[NkyLD][NxLD],
                                   const CComplex     f    [NsLB][NmLB][NzLB][NkyLD][NxLB  ][NvLB],
                                   const CComplex Field[Nq][NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                                   const double V[NvGB], const double M[NmGB])
{

    const double d6Z = M_PI * plasma->species[s].n0 * plasma->species[s].T0 * plasma->B0 * dv * grid->dm[m] * grid->dXYZ;

    ParticleFlux[:][:] = 0.;
        HeatFlux[:][:] = 0.;

    omp_for(int z=NzLlD; z<= NzLuD;z++) { 
  
      CComplex   ky_phi[NkyLD][NxLD], 
               Particle[NkyLD][NxLD], 
                 Energy[NkyLD][NxLD],
                      R[NkyLD][NxLD];

    omp_for(int y_k=NkyLlD; y_k<= NkyLuD; y_k++) { 
       
          // Geometry ?!
          const CComplex ky = ((CComplex) (0. + 1.j))  * fft->ky(y_k);
        
          for(int x=NxLlD; x <= NxLuD ; x++) { 

                ky_phi[y_k-NkyLlD][x-NxLlD] = ky * Field[Field::phi][s][m][z][y_k][x];

              // integrate over velocity space (calculate particle number and temperature)
                Energy[y_k-NkyLlD][x-NxLlD] = __sec_reduce_add((pow2(V[NvLlD:NvLD]) +  M[m] * plasma->B0) * f[s][m][z][y_k][x][NvLlD:NvLD]) * d6Z;
              Particle[y_k-NkyLlD][x-NxLlD] = __sec_reduce_add(                                             f[s][m][z][y_k][x][NvLlD:NvLD]) * d6Z;

      } } // x, y_k
              
         // Multiply electric field with temperature in real space 
         fft->multiply(ky_phi, Particle, R);
        
         // is atomic valid here (how about support, look gcc libatomic)?
         #pragma omp atomic
         ParticleFlux[:][:] += R[:][:];
         
         // Multiply electric field with temperature in real space 
         fft->multiply(ky_phi,   Energy, R);
         #pragma omp atomic
         HeatFlux[:][:]     += R[:][:];
      
    } // add for z,m
    
    return;


};

        
////////////////////////     Calculate x-dependent values     /////////////////////

void Analysis::initData(Setup *setup, FileIO *fileIO) {
        
     analysisGroup = fileIO->newGroup("Analysis");
        
     hsize_t offset0[] = { 0, 0, 0, 0, 0, 0, 0 };
    
     //---------------  Analysis - Heat fluxes ------------
     
     // Heat Flux ky and Particle FluxKy ( per species) 
     hid_t fluxGroup = fileIO->newGroup("Flux", analysisGroup);
     
     hsize_t FSky_dsdim[] = { Nq, Ns     , Nky   , Nx     , 1 }; 
     hsize_t FSky_dmdim[] = { Nq, Ns     , Nky   , Nx     , H5S_UNLIMITED} ;
     hsize_t FSky_cDdim[] = { Nq, NsLD   , Nky   , NxLD   , 1 };
     hsize_t FSky_cBdim[] = { Nq, NsLD   , Nky   , NxLD   , 1 };
     hsize_t FSky_cdoff[] = { 0,  NsLlB-1, NkyLlD, NxLlD-3, 0};

     FA_heatKy      = new FileAttr("Heat"    , fluxGroup, fileIO->file, 5, FSky_dsdim, FSky_dmdim, FSky_cDdim, offset0,  FSky_cBdim, FSky_cdoff, parallel->Coord[DIR_XYZ] == 0);
     FA_particleKy  = new FileAttr("Particle", fluxGroup, fileIO->file, 5, FSky_dsdim, FSky_dmdim, FSky_cDdim, offset0,  FSky_cBdim, FSky_cdoff, parallel->Coord[DIR_XYZ] == 0);
    
     H5Gclose(fluxGroup);
     
     //----------------  Moments - Heat fluxes ------------------
     
     hid_t momGroup = check(H5Gcreate(fileIO->getFileID(), "/Moments", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), DMESG("Error creating group file for Phi : H5Gcreate"));

     hsize_t mom_dsdim[] =  { grid->NzGD , NkyLD      , grid->NxGD , NsLD , 1            };
     hsize_t mom_dmdim[] =  { grid->NzGD , NkyLD      , grid->NxGD , Ns   , H5S_UNLIMITED};
     hsize_t mom_cBdim[] =  { NzLD       , NkyLD      , NxLD       , NsLD , 1            };
     hsize_t mom_cDdim[] =  { NzLD       , NkyLD      , NxLD       , NsLD , 1            };
     hsize_t mom_cmoff[] =  { 0, 0, 0, 0, 0                                              };
     hsize_t mom_cdoff[] =  { NzLlD-3    , 0          , NxLlD-3    , 0     ,  0          };
     
     bool momWrite = (parallel->Coord[DIR_VM] == 0);
     
     FA_Mom_Tp        = new FileAttr("Temperature_v", momGroup, fileIO->file, 5, mom_dsdim, mom_dmdim, mom_cDdim, mom_cmoff, mom_cBdim, mom_cdoff, momWrite, fileIO->complex_tid);
     FA_Mom_HeatFlux  = new FileAttr("HeatFlux"     , momGroup, fileIO->file, 5, mom_dsdim, mom_dmdim, mom_cDdim, mom_cmoff, mom_cBdim, mom_cdoff, momWrite, fileIO->complex_tid);
     FA_Mom_Density   = new FileAttr("Density"      , momGroup, fileIO->file, 5, mom_dsdim, mom_dmdim, mom_cDdim, mom_cmoff, mom_cBdim, mom_cdoff, momWrite, fileIO->complex_tid);
     FA_Mom_Time  = fileIO->newTiming(momGroup);
        
     H5Gclose(momGroup);
      
     dataOutputMoments   = Timing(setup->get("DataOutput.Moments.Step", -1)       , setup->get("DataOutput.Moments.Time", -1.));


     ////////////////////////////////////// X-Dependent data /////////////////////
     
     hid_t XDepGroup = fileIO->newGroup("XDep", analysisGroup);

     hsize_t XDep_dsdim[] =  { Ns     , Nx      , 1            };
     hsize_t XDep_dmdim[] =  { Ns     , Nx      , H5S_UNLIMITED};
     hsize_t XDep_cBdim[] =  { NsLD   , NxLD    , 1            };
     hsize_t XDep_cDdim[] =  { NsLD   , NxLD    , 1            };
     hsize_t XDep_cmoff[] =  { 0      , 0       , 0            };
     hsize_t XDep_cdoff[] =  { NsLlB-1, NxLlD-3 , 0            };
     
     bool XDepWrite = ( (parallel->Coord[DIR_VM] == 0) && (parallel->Coord[DIR_Z] == 0));
     
     FA_XDep_Tp   = new FileAttr("Temperature", XDepGroup, fileIO->file, 3, XDep_dsdim, XDep_dmdim, XDep_cDdim, XDep_cmoff, XDep_cBdim, XDep_cdoff, XDepWrite);
     FA_XDep_n    = new FileAttr("Density"    , XDepGroup, fileIO->file, 3, XDep_dsdim, XDep_dmdim, XDep_cDdim, XDep_cmoff, XDep_cBdim, XDep_cdoff, XDepWrite);
     FA_XDep_Time = fileIO->newTiming(XDepGroup);
        
     H5Gclose(XDepGroup);
     
     dataOutputXDep = Timing(setup->get("DataOutput.XDep.Step", -1), setup->get("DataOutput.XDep.Time", -1.));


     //---------------------------   Power Spectrum  -------------------------------
     
     // X-scalarValue
     hid_t growGroup = fileIO->newGroup("PowerSpectrum", analysisGroup);

     hsize_t grow_x_dsdim[] = { Nq, Nx, 1 }; 
     hsize_t grow_x_dmdim[] = { Nq, Nx, H5S_UNLIMITED} ;
     hsize_t grow_x_cDdim[] = { Nq, Nx, 1 };
     hsize_t grow_x_cBdim[] = { Nq, Nx, 1 };
     FA_grow_x  = new FileAttr("X", growGroup, fileIO->file, 3, grow_x_dsdim, grow_x_dmdim, grow_x_cDdim, offset0,  grow_x_cBdim, offset0, parallel->myRank == 0);

     // Y-scalarValue
     hsize_t grow_y_dsdim[] = { Nq, Nky, 1 };
     hsize_t grow_y_dmdim[] = { Nq, Nky, H5S_UNLIMITED };
     hsize_t grow_y_cDdim[] = { Nq, Nky, 1 };
     hsize_t grow_y_cBdim[] = { Nq, Nky, 1 };
   
     FA_grow_y  = new FileAttr("Y", growGroup, fileIO->file, 3, grow_y_dsdim, grow_y_dmdim, grow_y_cDdim, offset0,  grow_y_cBdim, offset0, parallel->myRank == 0);

     FA_grow_t  = fileIO->newTiming(growGroup);

     H5Gclose(growGroup);
    

     hid_t freqGroup = fileIO->newGroup("PhaseShift", analysisGroup);
     FA_freq_x  = new FileAttr("X", freqGroup, fileIO->file, 3, grow_x_dsdim, grow_x_dmdim, grow_x_cDdim, offset0,  grow_x_cBdim, offset0, parallel->myRank == 0);
     FA_freq_y  = new FileAttr("Y", growGroup, fileIO->file, 3, grow_y_dsdim, grow_y_dmdim, grow_y_cDdim, offset0,  grow_y_cBdim, offset0, parallel->myRank == 0);
     FA_freq_t  = fileIO->newTiming(freqGroup);
     H5Gclose(freqGroup);

     //////////////////////////////////////////////////////////////// Setup Table for scalar data ////////////////////////////////////////////////////////
              
     ScalarValues_t scalarValues;
     
     size_t SV_cdoff[] = { HOFFSET( ScalarValues_t, timestep       ), HOFFSET( ScalarValues_t, time     ), HOFFSET( ScalarValues_t, phiEnergy       ),
                           HOFFSET( ScalarValues_t, ApEnergy       ), HOFFSET( ScalarValues_t, BpEnergy ), HOFFSET( ScalarValues_t, particle_number ),
                           HOFFSET( ScalarValues_t, kinetic_energy ), HOFFSET( ScalarValues_t, entropy  ), HOFFSET( ScalarValues_t, heat_flux       ),
                           HOFFSET( ScalarValues_t, particle_flux  ) };

     size_t SV_sizes[] = { sizeof(scalarValues.timestep), sizeof(scalarValues.time    ), sizeof(scalarValues.phiEnergy), 
                           sizeof(scalarValues.ApEnergy), sizeof(scalarValues.BpEnergy), Ns * sizeof(scalarValues.particle_number[0]), Ns * sizeof(scalarValues.kinetic_energy[0]), 
                           Ns * sizeof(scalarValues.entropy[0]), Ns * sizeof(scalarValues.heat_flux[0]), Ns * sizeof(scalarValues.particle_flux[0])};

     hid_t SV_types[] = { H5T_NATIVE_INT, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, 
                          fileIO->species_tid, fileIO->species_tid, fileIO->species_tid, fileIO->species_tid, fileIO->species_tid } ;
  
     const char *SV_names[] = { "Timestep", "Time", "phiEnergy", "ApEnergy", "BpEnergy", "ParticleNumber", "KineticEnergy", "Entropy", "HeatFlux", "ParticleFlux" };

     SVTable = new TableAttr(analysisGroup, "scalarValues", 10, SV_names, SV_cdoff, SV_types, SV_sizes, &scalarValues); 

     dataOutputStatistics  = Timing(setup->get("DataOutput.Statistics.Step", -1), setup->get("DataOutput.Statistics.Time", -1.));

}


void Analysis::getFieldEnergy(double& phiEnergy, double& ApEnergy, double& BpEnergy)
{

   fields->getFieldEnergy(phiEnergy, ApEnergy, BpEnergy);

};
  
void Analysis::writeData(const Timing &timing, const double dt)

{
  
  ////////////////// Get Moments /////////////////////

  if (timing.check(dataOutputMoments, dt)       )   {

               CComplex 
               A4_Tp[NsLD][NzLD][NkyLD][NxLD], 
               A4_n [NsLD][NzLD][NkyLD][NxLD], 
               A4_Xi[NsLD][NzLD][NkyLD][NxLD], 
               A4_Q [NsLD][NzLD][NkyLD][NxLD]; 


             getTemperature  ((A6zz) vlasov->f, (A4zz) A4_Tp, V, M);
             getNumberDensity((A6zz) vlasov->f, (A4zz) A4_n , V, M);
  
             FA_Mom_Density->write((CComplex *) A4_n);
             FA_Mom_Tp     ->write((CComplex *) A4_Tp);
             FA_Mom_Time->write(&timing);

    parallel->print("Data I/O : Moments output");
    
  }


  ////////////////// Store X-dependent data /////////////
  if (timing.check(dataOutputXDep, dt)       )   {

    {   
      CComplex A4_T[NsLD][NzLD][NkyLD][NxLD], 
               A4_p[NsLD][NzLD][NkyLD][NxLD], 
               A4_n[NsLD][NzLD][NkyLD][NxLD]; 
      double   A1_T[NsLD][NxLD], 
               A1_n[NsLD][NxLD], 
               A1_p[NsLD][NxLD];

      getNumberDensity( (A6zz) vlasov->f, (A4zz) A4_n, V, M);
      getTemperature  ( (A6zz) vlasov->f, (A4zz) A4_T, V, M);

      // Reduce over y_k and z
      for(int s = NsLlD; s <= NsLuD; s++) {  for(int x = NxLlD;  x <= NxLuD; x++) {

         A1_n[s-NsLlD][x-NxLlD] = __sec_reduce_add(cabs(A4_n[s-NsLlD][:][:][x-NxLlD]));
         A1_T[s-NsLlD][x-NxLlD] = __sec_reduce_add(cabs(A4_T[s-NsLlD][:][:][x-NxLlD]));
      } } 
   
      FA_XDep_Tp  ->write((double *) A1_T); 
      FA_XDep_n   ->write((double *) A1_n); 
      FA_XDep_Time->write(&timing);

      parallel->print("Data I/O : X-Dep output");
    }  

   ///////////////////Get Store /Particle Flux /////////////////////////////
   {
       Complex A4_q[Nq][NsLD][NkyLD][NxLD], A2_q[NkyLD][NxLD], // Heat flux      (take care stack variables have no offset)
               A4_x[Nq][NsLD][NkyLD][NxLD], A2_x[NkyLD][NxLD]; // Particle flux

        A4_q[:][:][:][:] = 0.;
        A4_x[:][:][:][:] = 0.;

      for(int s = NsLlD; s <= NsLuD; s++) {  for(int m = NmLlD;  m <= NmLuD; m++) {

           // support only phi for now
           getParticleHeatFlux(m, s, (A2zz) A2_x, (A2zz) A2_q, (A6zz) vlasov->f, (A6zz) fields->Field, V, M);
           A4_q[0][s-NsLlD][:][:] +=  A2_q[:][:];
           A4_x[0][s-NsLlD][:][:] +=  A2_x[:][:];

      } };

      parallel->reduce((CComplex *) A4_q, Op::SUM, DIR_M, Nq*NsLD*NkyLD*NxLD);
      parallel->reduce((CComplex *) A4_x, Op::SUM, DIR_M, Nq*NsLD*NkyLD*NxLD);

      FA_heatKy    ->write((CComplex *) A4_q);
      FA_particleKy->write((CComplex *) A4_x);
    }
  }

  if (timing.check(dataOutputStatistics, dt)       )   {
 

  
    ////////////// Scalar Variables /////////////////

    // Stack allocation (size ok ?) and alignement ? (well, speed not critical here...)
    double pSpecX [Nq][Nx], pSpecY [Nq][Nky],
           pPhaseX[Nq][Nx], pPhaseY[Nq][Nky];
     
    getPowerSpectrum((A4zz) fft->kXOut, (A4zz) fields->Field0, pSpecX, pSpecY, pPhaseX, pPhaseY);
    
    // Seperatly writing ? Hopefully it is buffered ... (passing stack pointer ... OK ?)
    // Memory leak !?
    FA_grow_x->write( &pSpecX [0][0]); FA_grow_y->write(&pSpecY [0][0]); FA_grow_t->write(&timing);
    FA_freq_x->write( &pPhaseX[0][0]); FA_freq_y->write(&pPhaseY[0][0]); FA_freq_t->write(&timing);
    
    ScalarValues scalarValues;
    
    // calculate kinetic Energy first, need for initial_e ! sum over domain
    scalarValues.timestep = timing.step;
    scalarValues.time     = timing.time;
    
    fields->getFieldEnergy(scalarValues.phiEnergy, scalarValues.ApEnergy, scalarValues.BpEnergy);
    
    //  Get scalar Values for every species ( this is bad calculate them alltogether)
    // gives bus error ?!
    calculateScalarValues((A6zz) vlasov->f, (A6zz) vlasov->f0, V, M, 1, scalarValues); 

    
    SVTable->append(&scalarValues);
    
    // write out to Terminal/File
    
    std::stringstream messageStream;
    messageStream << std::endl << std::endl << "Analysis | " << std::setprecision(3);
    messageStream << "Field Energy : (phi) " << scalarValues.phiEnergy  << "  (Ap) " << scalarValues.ApEnergy  <<  "  (Bp) " << scalarValues.BpEnergy << std::endl; 
    double charge = 0., kinetic_energy=0.;
    
    for(int s = NsGlD; s <= NsGuD; s++) {
    
      messageStream << "         | " << plasma->species[s].name 
                     << " N : "             << scalarValues.particle_number[s-1]  
                     << " Kinetic Energy: " << scalarValues.kinetic_energy[s-1] 
                     << " Particle Flux :"  << scalarValues.particle_flux[s-1]  
                     << " Heat Flux : "     << scalarValues.heat_flux[s-1] << std::endl;
        charge += plasma->species[s].q  * scalarValues.particle_number[s-1];
      kinetic_energy += scalarValues.kinetic_energy[s-1];
    }
    
    messageStream << //"------------------------------------------------------------------" <<
      "         | Total Energy " << kinetic_energy+scalarValues.phiEnergy + scalarValues.ApEnergy + scalarValues.BpEnergy << "    Total Charge = " << ((plasma->species[0].n0 != 0.) ? 0. : charge) 
      << std::endl;  
    parallel->print(messageStream.str());
    
  }
  
  return;
  
}

     //###################################### Spectrum  ################################################
/*
     // XY-Spectrum
     hsize_t spec_xy_dsdim[]      = { Nx/2+1, Ny/2+1, 1};
     hsize_t spec_xy_dmdim[]   = { Nx/2+1, Ny/2+1, H5S_UNLIMITED};
     hsize_t spec_xy_cDdim[] = { kNxL, kNyL, 1};

     FA_spec_xy  = new FileAttr("Unnamed",3, spec_xy_dsdim, spec_xy_dmdim, spec_xy_cDdim, offset0,  spec_xy_cDdim, parallel->isFFTGroup);
     
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
     
  delete FA_XDep_Tp; delete FA_XDep_n; delete FA_XDep_Time;
       
  delete FA_Mom_Tp;
  delete FA_Mom_HeatFlux;
  delete FA_Mom_Density;
  delete FA_Mom_Time;

  delete SVTable;
  H5Gclose(analysisGroup);


}
    
void Analysis::printOn(std::ostream &output)  const
{ 



};
