/*
 * =====================================================================================
 *
 *       Filename: Diagnostics.cpp
 *
 *    Description: Diagnostic functions for gkc
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include "Diagnostics.h"
#include "Tools/System.h"

Diagnostics::Diagnostics(Parallel *_parallel, Vlasov *_vlasov, Fields *_fields, Grid *_grid, Setup *_setup, FFTSolver *_fft, FileIO *fileIO, Geometry *_geo) : 
 
parallel(_parallel),setup(_setup), vlasov(_vlasov), grid(_grid), fields(_fields), geo(_geo),  fft(_fft)

{

  initData(setup, fileIO);

  moments   = new Moments(setup, vlasov, fields, grid, parallel);
  auxiliary = new Auxiliary(setup, fileIO, parallel, fields, vlasov, grid, fft);
}

Diagnostics::~Diagnostics() 
{
  delete auxiliary;
  delete moments;
  closeData();
}


void Diagnostics::getPowerSpectrum(CComplex  kXOut  [Nq][NzLD][Nky][FFTSolver::X_NkxL], 
                                CComplex  Field0 [Nq][NzLD][Nky][NxLD]      ,
                                double    pSpecX [Nq][Nx/2+1], double pSpecY [Nq][Nky],
                                double    pPhaseX[Nq][Nx/2+1], double pPhaseY[Nq][Nky])
{

  double   pSpec[Nq][Nx]; pSpec[:][:] = 0.;
  CComplex pFreq[Nq][Nx]; pFreq[:][:] = 0.;
  // Note : take care that the FFT output is of for e.g. an 8 number sequence [0, 1, 2, 3, 4,-3,-2,-1]
  // Note : atan2 is not a linear function, thus phase has to be calculated at last
  if(parallel->Coord[DIR_VMS] == 0) {

    // Note : We have domain decomposition in X but not in Y
    fft->solve(FFT_Type::X_FIELDS, FFT_Sign::Forward, &Field0[0][NzLlD][NkyLlD][NxLlD]);
         
    for(int q = 0; q < Nq; q++) {
                
      // Mode power & Phase shifts for X (domain decomposed) |\phi(k_x)|
      for(int x_k=fft->K1xLlD; x_k <= fft->K1xLuD; x_k++) {

        pSpec[q][x_k] = __sec_reduce_add(cabs(kXOut[q][NzLlD:NzLD][NkyLlD:Nky][x_k]))/fft->Norm_X * dy * dz; 
        pFreq[q][x_k] = __sec_reduce_add(     kXOut[q][NzLlD:NzLD][NkyLlD:Nky][x_k] )/fft->Norm_X * dy * dz;

      }
      
      // Mode power & Phase shifts for Y (not decomposed) |\phi(k_y)|
      for(int y_k = NkyLlD; y_k <=  NkyLuD ; y_k++) { 
        
        pSpecY [q][y_k] =      __sec_reduce_add(cabs(Field0[q][NzLlD:NzLD][y_k][NxLlD:NxLD])) * dx * dz; 
        pPhaseY[q][y_k] = carg(__sec_reduce_add(     Field0[q][NzLlD:NzLD][y_k][NxLlD:NxLD])) * dx * dz;

      }
      
    }
             
    // Finalized calculations for x-domain (decomposed and negative frequency modes)

    // Sum up with X-values from other processes
    parallel->reduce(&pSpec[0][0], Op::sum, DIR_XYZ, Nq * Nx);        
    parallel->reduce(&pFreq[0][0], Op::sum, DIR_XYZ, Nq * Nx);        
    
    // Sum up with Y-values from other processes
    parallel->reduce(&pSpecY [0][0], Op::sum, DIR_XYZ, Nq * Nky);        
    parallel->reduce(&pPhaseY[0][0], Op::sum, DIR_XYZ, Nq * Nky);        
         
         
    // map back from fftw [k0, k1, k2, ..., k_Ny, -k_(N/y-2), ... -k_1]
    for(int x_k = 1; x_k < Nx/2; x_k++) pFreq [:][x_k] += pFreq[:][Nx-x_k];
    for(int x_k = 1; x_k < Nx/2; x_k++) pSpec [:][x_k] += pSpec[:][Nx-x_k];

    pPhaseX[:][0:Nx/2+1] = carg(pFreq[:][0:Nx/2+1]);
    pSpecX [:][0:Nx/2+1] =      pSpec[:][0:Nx/2+1];
         
  } // if DIR_XYZ
}

//////////////////////// Calculate scalar values ///////////////////////////

void Diagnostics::calculateScalarValues(const CComplex f [NsLD][NmLD][NzLB][Nky][NxLB][NvLB], 
                                        const CComplex f0[NsLD][NmLD][NzLB][Nky][NxLB][NvLB],
                                        const CComplex Mom[8][NsLD][NzLD][Nky][NxLD], 
                                        const double   ParticleFlux[Nq][NsLD][Nky][NxLD], 
                                        const double   HeatFlux[Nq][NsLD][Nky][NxLD],
                                        ScalarValues &scalarValues) 

{
  
  for(int s = NsGlD; s <= NsGuD; s++) {
    
      
    double number        = 0.;
    double kineticEnergy = 0.;
    double entropy       = 0.;
    double particle[Nq]; particle[:] = 0.;
    double heat[Nq];         heat[:] = 0.;
    
    if((s >= NsLlD && s <= NsLuD)) {
      
    ////////////////////////////// Calculate Entropy /////////////////////////////////////

    // Kinetic energy is calculated in gyro-center coordinates
    // Y.Idomura et al., J.Comp.Phys 2007, New conservative gk ..., Eq.(11)  [ but we use different normalization]
    // Only y_k==0 has contributions, as other one cancel with integration over y
    for(int m = NmLlD; m <= NmLuD; m++) {

    const double d6Z = grid->dXYZ * dv * grid->dm[m] * pow2(species[s].m) * plasma->B0;
    for(int z = NzLlD; z <= NzLuD; z++) { for(int x = NxLlD; x <= NxLuD; x++) { 
      
      kineticEnergy += creal(__sec_reduce_add((species[s].m * pow2(V[NvLlD:NvLD]) + M[m] * plasma->B0) * f[s][m][z][0][x][NvLlD:NvLD] * d6Z));
    
    } } // z, y_k, x
    
    // calculation through moments, (which one is correct ?)
    //    kineticEnergy = (species[s].m * __sec_reduce_add(creal(Mom[1][s-NsLlD][0:NzLD][0][0:NxLD])) +
    //                                    __sec_reduce_add(creal(Mom[2][s-NsLlD][0:NzLD][0][0:NxLD])) ) * grid->dXYZ;
    //kineticEnergy = (0.5 * species[s].m * __sec_reduce_add(creal(Mom[1][s-NsLlD][0:NzLD][0][0:NxLD])) +
    //                                      __sec_reduce_add(creal(Mom[2][s-NsLlD][0:NzLD][0][0:NxLD])) ) * grid->dXYZ;
     
    //#pragma omp parallel for reduction(+:entropy) collapse (2)
    for(int z = NzLlD; z <= NzLuD; z++) { for(int x = NxLlD; x <= NxLuD; x++) { 

      entropy += creal(pow2(__sec_reduce_add(f [s][m][z][0][x][NvLlD:NvLD])))/
                            __sec_reduce_add(f0[s][m][z][0][x][NvLlD:NvLD]);
    } } // z, y_k, x

    } // m 
    
    //// Calculate Total Heat & Particle Flux (reduction over moments variables (no v_par, no mu))
    if (parallel->Coord[DIR_VM] == 0) { 
   
      ////////////////////////////// Calculate Particle Number /////////////////////////////
      number =  creal(__sec_reduce_add(Mom[0][s-NsLlD][0:NzLD][0][0:NxLD]));
      
      for(int q = 0; q < Nq; q++) {
        particle[q] = __sec_reduce_add(ParticleFlux[q][s-NsLlD][:][:]);
        heat    [q] = __sec_reduce_add(    HeatFlux[q][s-NsLlD][:][:]);
      } 
    }
    
    } // if(local s) 

    // Communicate with other groups (reduce whole structure)
    parallel->reduce(particle     , Op::sum, DIR_ALL, Nq);
    parallel->reduce(heat         , Op::sum, DIR_ALL, Nq);
    kineticEnergy = parallel->reduce(kineticEnergy, Op::sum);
    number        = parallel->reduce(number       , Op::sum);

    {
       scalarValues.particle_number[s-1]  = number        ;
       scalarValues.entropy        [s-1]  = entropy        ;
       scalarValues.kinetic_energy [s-1]  = kineticEnergy;
       scalarValues.particle_flux  [Nq*(s-1):Nq]  = particle[0:Nq];
       scalarValues.heat_flux      [Nq*(s-1):Nq]  = heat[0:Nq] ;
    }  

    } // s
}

void Diagnostics::getParticleHeatFlux( 
                                   double    ParticleFlux[Nq]   [NsLD][Nky][NxLD], 
                                   double        HeatFlux[Nq]   [NsLD][Nky][NxLD],
                                   double      CrossPhase[Nq][3][NsLD][Nky][NxLD],
                                   const CComplex  Field0[Nq][NzLD][Nky][NxLD],
                                   const CComplex Mom[8][NsLD][NzLD][Nky][NxLD] 
                                  )
{

  ParticleFlux[:][:][:][:]    = 0.;
      HeatFlux[:][:][:][:]    = 0.;
    CrossPhase[:][:][:][:][:] = 0.;

  // Triad condition. Heat/Particles are only transported by the y_k = 0, as the other y_k > 0 
  // modes cancels out. Thus multiplying y_k = 0 = A(y_k)*B(-y_k) = A(y_k)*[cc : B(y_k)]
  // where the complex conjugate values is used as physical value is a real quantity.
  
  // Note heat/particle fluxes are summed over z, cross-phases are averaged
  
  if(parallel->Coord[DIR_V] == 0) { for(int s = NsLlD; s <= NsLuD; s++) {  
  
  double norm[3] = { 1., -species[s].v_th, species[s].T0 / (species[s].q * plasma->B0) };

  norm[:] *= species[s].n0 / geo->C;
  
  // Heat/Particle fluxes are calculates for each field quantity (phi, A_par, B_par)
  for(int q = 0    ; q <  Nq   ; q++) {
 
  for(int z = NzLlD; z <= NzLuD; z++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { 
  for(int x = NxLlD; x <= NxLuD; x++) { 
      
    // Do I have to include geometry terms ?! partial_y (phi, A_par, B_par)
    const CComplex iky_field  = -_imag * fft->ky(y_k) * Field0[q][z][y_k][x];

    // take only real part as it gives the radial direction ?!
    ParticleFlux[q][s-NsLlD][y_k][x-NxLlD] += norm[q] * creal( iky_field * conj(Mom[3*q+0][s-NsLlD][z-NzLlD][y_k][x-NxLlD]));
        HeatFlux[q][s-NsLlD][y_k][x-NxLlD] += norm[q] * creal( iky_field * conj(Mom[3*q+1][s-NsLlD][z-NzLlD][y_k][x-NxLlD]
                                                                             + Mom[3*q+2][s-NsLlD][z-NzLlD][y_k][x-NxLlD]));

    // Get Cross-phases for the (phi, Ap, Bp) x ( M00, T, T) ( need to divide over Nz to get average after reduce operator later)
    CrossPhase[q][0][s-NsLlD][y_k][x-NxLlD] += carg( iky_field * conj(Mom[0][s-NsLlD][z-NzLlD][y_k][x-NxLlD])) / Nz;
    CrossPhase[q][1][s-NsLlD][y_k][x-NxLlD] += carg( iky_field * conj(Mom[1][s-NsLlD][z-NzLlD][y_k][x-NxLlD])) / Nz;
    CrossPhase[q][2][s-NsLlD][y_k][x-NxLlD] += carg( iky_field * conj(Mom[2][s-NsLlD][z-NzLlD][y_k][x-NxLlD])) / Nz;
    
  }    // x 
  } }  // y_k, z
  }    // q

  } }  // s,v

  parallel->reduce(&ParticleFlux[0][0][0][0] , Op::sum, DIR_Z, Nq * NsLD * Nky * NxLD);
  parallel->reduce(&    HeatFlux[0][0][0][0] , Op::sum, DIR_Z, Nq * NsLD * Nky * NxLD);

  parallel->reduce(&CrossPhase[0][0][0][0][0], Op::sum, DIR_Z, Nq * NsLD * Nky * NxLD * 3);
}

        
void Diagnostics::initData(Setup *setup, FileIO *fileIO) 
{
        
  analysisGroup = fileIO->newGroup("Analysis");
        
  hsize_t offset0[] = { 0, 0, 0, 0, 0, 0, 0 };
    
  //----------------  Moments - Heat fluxes ------------------
     
  hid_t momGroup = fileIO->newGroup("Moments");

  hsize_t mom_dsdim[] =  { Ns     , Nz     , Nky, Nx     , 1             };
  hsize_t mom_dmdim[] =  { Ns     , Nz     , Nky, Nx     , H5S_UNLIMITED };
  hsize_t mom_cBdim[] =  { NsLD   , NzLD   , Nky, NxLD   , 1             };
  hsize_t mom_cDdim[] =  { NsLD   , NzLD   , Nky, NxLD   , 1             };
  hsize_t mom_cmoff[] =  { 0      , 0      , 0  , 0      , 0             };
  hsize_t mom_cdoff[] =  { NsLlD-1, NzLlD-3, 0  , NxLlD-3, 0             };
     
  bool momWrite = (parallel->Coord[DIR_VM] == 0);
     
  FA_Mom_00       = new FileAttr("Mom00"   , momGroup, fileIO->file, 5, mom_dsdim, mom_dmdim, mom_cDdim, mom_cmoff, mom_cBdim, mom_cdoff, momWrite, fileIO->complex_tid);
  FA_Mom_20       = new FileAttr("Mom20"   , momGroup, fileIO->file, 5, mom_dsdim, mom_dmdim, mom_cDdim, mom_cmoff, mom_cBdim, mom_cdoff, momWrite, fileIO->complex_tid);
  FA_Mom_02       = new FileAttr("Mom02"   , momGroup, fileIO->file, 5, mom_dsdim, mom_dmdim, mom_cDdim, mom_cmoff, mom_cBdim, mom_cdoff, momWrite, fileIO->complex_tid);
  
  FA_Mom_10       = new FileAttr("Mom10"   , momGroup, fileIO->file, 5, mom_dsdim, mom_dmdim, mom_cDdim, mom_cmoff, mom_cBdim, mom_cdoff, momWrite, fileIO->complex_tid);
  FA_Mom_30       = new FileAttr("Mom30"   , momGroup, fileIO->file, 5, mom_dsdim, mom_dmdim, mom_cDdim, mom_cmoff, mom_cBdim, mom_cdoff, momWrite, fileIO->complex_tid);
  FA_Mom_12       = new FileAttr("Mom12"   , momGroup, fileIO->file, 5, mom_dsdim, mom_dmdim, mom_cDdim, mom_cmoff, mom_cBdim, mom_cdoff, momWrite, fileIO->complex_tid);

  FA_Mom_HeatFlux = new FileAttr("HeatFlux", momGroup, fileIO->file, 5, mom_dsdim, mom_dmdim, mom_cDdim, mom_cmoff, mom_cBdim, mom_cdoff, momWrite);
  FA_Mom_PartFlux = new FileAttr("PartFlux", momGroup, fileIO->file, 5, mom_dsdim, mom_dmdim, mom_cDdim, mom_cmoff, mom_cBdim, mom_cdoff, momWrite);
                                   

  FA_Mom_Time     = fileIO->newTiming(momGroup);
        
  H5Gclose(momGroup);
      
  dataOutputMoments   = Timing(setup->get("DataOutput.Moments.Step", -1), setup->get("DataOutput.Moments.Time", -1.));

  ////////////////////////////////////// X-Dependent Moments data /////////////////////
  { 
    hid_t XDepGroup = fileIO->newGroup("XDep", analysisGroup);

    bool XDepWrite  = (parallel->Coord[DIR_VM] == 0) && (parallel->Coord[DIR_Z] == 0);

    hsize_t XDep_dsdim[] =  { 8, Ns     , Nx      , 1            };
    hsize_t XDep_dmdim[] =  { 8, Ns     , Nx      , H5S_UNLIMITED};
    hsize_t XDep_cBdim[] =  { 8, NsLD   , NxLD    , 1            };
    hsize_t XDep_cDdim[] =  { 8, NsLD   , NxLD    , 1            };
    hsize_t XDep_cdoff[] =  { 0, NsLlB-1, NxLlD-3 , 0            };
     
     
    FA_XDep_Mom  = new FileAttr("Mom", XDepGroup, fileIO->file, 4, XDep_dsdim, XDep_dmdim, XDep_cDdim, offset0, XDep_cBdim, XDep_cdoff, XDepWrite);
    FA_XDep_Time = fileIO->newTiming(XDepGroup);
        
    H5Gclose(XDepGroup);
     
    dataOutputXDep = Timing(setup->get("DataOutput.XDep.Step", -1), setup->get("DataOutput.XDep.Time", -1.));
  }

  //////////////////////////////////////  Analysis - Heat fluxes     /////////////////////////////////////////
  {
    hid_t fluxGroup = fileIO->newGroup("Flux", analysisGroup);
    bool isXYZ  = (parallel->Coord[DIR_VM] == 0) && (parallel->Coord[DIR_Z] == 0);  // we average over z
   
    // Heat Flux ky and Particle FluxKy ( per species) 
    hsize_t FSky_dsdim[] = { Nq, Ns     , Nky   , Nx     , 1            }; 
    hsize_t FSky_dmdim[] = { Nq, Ns     , Nky   , Nx     , H5S_UNLIMITED};
    hsize_t FSky_cDdim[] = { Nq, NsLD   , Nky   , NxLD   , 1            };
    hsize_t FSky_cBdim[] = { Nq, NsLD   , Nky   , NxLD   , 1            };
    hsize_t FSky_cdoff[] = { 0 , NsLlB-1, 0     , NxLlD-3, 0            };

    FA_HeatFluxKy = new FileAttr("Heat"    , fluxGroup, fileIO->file, 5, FSky_dsdim, FSky_dmdim, FSky_cDdim, offset0,  FSky_cBdim, FSky_cdoff, isXYZ);
    FA_PartFluxKy = new FileAttr("Particle", fluxGroup, fileIO->file, 5, FSky_dsdim, FSky_dmdim, FSky_cDdim, offset0,  FSky_cBdim, FSky_cdoff, isXYZ);
   

    // Data-Output for cross phase calculations for (phi, Ap, Bp) x (n, Tp, To) & (Ky, Species, X)
    hsize_t cph_dsdim[] = { Nq, 3, Ns     , Nky   , Nx     , 1            }; 
    hsize_t cph_dmdim[] = { Nq, 3, Ns     , Nky   , Nx     , H5S_UNLIMITED};
    hsize_t cph_cDdim[] = { Nq, 3, NsLD   , Nky   , NxLD   , 1            };
    hsize_t cph_cBdim[] = { Nq, 3, NsLD   , Nky   , NxLD   , 1            };
    hsize_t cph_cdoff[] = { 0 , 0, NsLlB-1, 0     , NxLlD-3, 0            };
  
    FA_CrossPhase = new FileAttr("CrossPhases", fluxGroup, fileIO->file, 6, cph_dsdim, cph_dmdim, cph_cDdim, offset0, cph_cBdim, cph_cdoff, isXYZ);
    
    H5Gclose(fluxGroup);
  }

  //---------------------------   Power Spectrum  -------------------------------
  { 
    hid_t growGroup = fileIO->newGroup("PowerSpectrum", analysisGroup);
    
    bool isMaster = parallel->myRank == 0;

    // X-scalarValue
    hsize_t grow_x_dsdim[] = { Nq, Nx, 1 }; 
    hsize_t grow_x_dmdim[] = { Nq, Nx, H5S_UNLIMITED} ;
    hsize_t grow_x_cDdim[] = { Nq, Nx, 1 };
    hsize_t grow_x_cBdim[] = { Nq, Nx, 1 };
    FA_grow_x  = new FileAttr("X", growGroup, fileIO->file, 3, grow_x_dsdim, grow_x_dmdim, grow_x_cDdim, offset0,  grow_x_cBdim, offset0, isMaster);

    // Y-scalarValue
    hsize_t grow_y_dsdim[] = { Nq, Nky, 1 };
    hsize_t grow_y_dmdim[] = { Nq, Nky, H5S_UNLIMITED };
    hsize_t grow_y_cDdim[] = { Nq, Nky, 1 };
    hsize_t grow_y_cBdim[] = { Nq, Nky, 1 };
   
    FA_grow_y  = new FileAttr("Y", growGroup, fileIO->file, 3, grow_y_dsdim, grow_y_dmdim, grow_y_cDdim, offset0,  grow_y_cBdim, offset0, isMaster);

    FA_grow_t  = fileIO->newTiming(growGroup);
  
    H5Gclose(growGroup);
    
    hid_t freqGroup = fileIO->newGroup("PhaseShift", analysisGroup);
    FA_freq_x  = new FileAttr("X", freqGroup, fileIO->file, 3, grow_x_dsdim, grow_x_dmdim, grow_x_cDdim, offset0,  grow_x_cBdim, offset0, isMaster);
    FA_freq_y  = new FileAttr("Y", growGroup, fileIO->file, 3, grow_y_dsdim, grow_y_dmdim, grow_y_cDdim, offset0,  grow_y_cBdim, offset0, isMaster);
    FA_freq_t  = fileIO->newTiming(freqGroup);
    H5Gclose(freqGroup);
  }
  //////////////////////////////////////////////////////////////// Setup Table for scalar data ////////////////////////////////////////////////////////
              
  ScalarValues_t scalarValues;
     
  size_t SV_cdoff[] = { HOFFSET( ScalarValues_t, walltime       ),
                        HOFFSET( ScalarValues_t, timestep       ), HOFFSET( ScalarValues_t, time     ), HOFFSET( ScalarValues_t, dt),
                        HOFFSET( ScalarValues_t, phiEnergy       ),
                        HOFFSET( ScalarValues_t, ApEnergy       ), HOFFSET( ScalarValues_t, BpEnergy ), HOFFSET( ScalarValues_t, particle_number ),
                        HOFFSET( ScalarValues_t, kinetic_energy ), HOFFSET( ScalarValues_t, entropy  ), HOFFSET( ScalarValues_t, heat_flux       ),
                        HOFFSET( ScalarValues_t, particle_flux  ) };

  size_t SV_sizes[] = { sizeof(scalarValues.walltime),
                        sizeof(scalarValues.timestep), sizeof(scalarValues.time    ), sizeof(scalarValues.dt),
                        sizeof(scalarValues.phiEnergy), 
                        sizeof(scalarValues.ApEnergy), sizeof(scalarValues.BpEnergy), Ns * sizeof(scalarValues.particle_number[0]), 
                        Ns * sizeof(scalarValues.kinetic_energy[0]), Ns * sizeof(scalarValues.entropy[0]), 
                        Ns * Nq * sizeof(scalarValues.heat_flux[0]), Ns * Nq * sizeof(scalarValues.particle_flux[0])};

  hid_t SV_types[] = { H5T_NATIVE_DOUBLE, H5T_NATIVE_INT, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, 
                       fileIO->species_tid, fileIO->species_tid, fileIO->species_tid, fileIO->specfield_tid, fileIO->specfield_tid } ;
  
  const char *SV_names[] = { "WallTime", "Timestep", "Time", "dt", "phiEnergy", "ApEnergy", "BpEnergy", "ParticleNumber", 
                             "KineticEnergy", "Entropy", "HeatFlux", "ParticleFlux" };

  SVTable = new TableAttr(analysisGroup, "scalarValues", 12, SV_names, SV_cdoff, SV_types, SV_sizes, &scalarValues); 

  dataOutputStatistics  = Timing(setup->get("DataOutput.Statistics.Step", -1), setup->get("DataOutput.Statistics.Time", -1.));
}

void Diagnostics::getFieldEnergy(double& phiEnergy, double& ApEnergy, double& BpEnergy)
{

  fields->getFieldEnergy(phiEnergy, ApEnergy, BpEnergy);

}
  
void Diagnostics::writeData(const Timing &timing, const double dt)

{

  if (timing.check(dataOutputMoments   , dt) || timing.check(dataOutputXDep, dt) ||
      timing.check(dataOutputStatistics, dt)) 
  {

    CComplex  Mom[8][NsLD][NzLD][Nky][NxLD];
        
    double  HeatFlux[Nq]   [NsLD][Nky][NxLD],  // Heat     flux   
            PartFlux[Nq]   [NsLD][Nky][NxLD],  // Particle flux
          CrossPhase[Nq][3][NsLD][Nky][NxLD];  // Cross phases of fluxes

    // Get Moments of Vlasov equation
    moments->getMoments((A6zz) vlasov->f, (A4zz) fields->Field0, Mom);
        
    getParticleHeatFlux(PartFlux, HeatFlux, CrossPhase, (A4zz) fields->Field0, Mom);

    ////////////////// Output Moments /////////////////////

    if (timing.check(dataOutputMoments, dt) )   {
      
      FA_Mom_HeatFlux->write((double *) &HeatFlux[0][0][0][0]);
      FA_Mom_PartFlux->write((double *) &PartFlux[0][0][0][0]);

      FA_Mom_00->write((CComplex *) &Mom[0][0][0][0][0]);
      FA_Mom_20->write((CComplex *) &Mom[1][0][0][0][0]);
      FA_Mom_02->write((CComplex *) &Mom[2][0][0][0][0]);
      
      FA_Mom_10->write((CComplex *) &Mom[3][0][0][0][0]);
      FA_Mom_30->write((CComplex *) &Mom[4][0][0][0][0]);
      FA_Mom_12->write((CComplex *) &Mom[5][0][0][0][0]);
      
      FA_Mom_Time->write(&timing);

      parallel->print("Data I/O : Moments output");
    }

    ////////////////// Store X-dependent data /////////////
    if (timing.check(dataOutputXDep, dt) )   {

      FA_PartFluxKy->write((double *) PartFlux);
      FA_HeatFluxKy->write((double *) HeatFlux);
      FA_CrossPhase->write((double *) CrossPhase);
      
      double Mom_XDep[8][NsLD][NxLD]; 

      // Reduce moments over y_k and z
      for(int n = 0    ; n <      8; n++) {
      for(int s = NsLlD; s <= NsLuD; s++) {  for(int x = NxLlD; x <= NxLuD; x++) {

        Mom_XDep[n][s-NsLlD][x-NxLlD] = __sec_reduce_add(creal(Mom[n][s-NsLlD][:][0][x-NxLlD]));

      } } }
   
      FA_XDep_Mom ->write((double *) &Mom_XDep); 
      FA_XDep_Time->write(&timing);

      parallel->print("Data I/O : X-Dep output");

    }  

    ////////////// Scalar Variables /////////////////
    if (timing.check(dataOutputStatistics, dt) )   {
    { 
      // calculate mode spectrum of fields (phi, Ap, Bp)

      double pSpecX [Nq][Nx], pSpecY [Nq][Nky],
             pPhaseX[Nq][Nx], pPhaseY[Nq][Nky];
     
      getPowerSpectrum((A4zz) fft->kXOut, (A4zz) fields->Field0, pSpecX, pSpecY, pPhaseX, pPhaseY);
    
      // Separately writing ? Hopefully it is buffered ... (passing stack pointer ... OK ?)
      FA_grow_x->write( &pSpecX [0][0]); FA_grow_y->write(&pSpecY [0][0]); FA_grow_t->write(&timing);
      FA_freq_x->write( &pPhaseX[0][0]); FA_freq_y->write(&pPhaseY[0][0]); FA_freq_t->write(&timing);
    }

    ScalarValues scalarValues;
   

    static struct timeval walltime = System::getTimeOfDay();  

    // Need to communicate wall clock time to other processes (as HDF-5 table requires same value)
    scalarValues.walltime = parallel->bcast(System::getTimeDifference(walltime), parallel->myRank == 0) ;
    scalarValues.timestep = timing.step;
    scalarValues.time     = timing.time;
    scalarValues.dt       = dt;
    
    fields->getFieldEnergy(scalarValues.phiEnergy, scalarValues.ApEnergy, scalarValues.BpEnergy);
    
    //  Get scalar values for every species ( this is bad, should calculate them all together)
    calculateScalarValues((A6zz) vlasov->f, (A6zz) vlasov->f0,
                           Mom, PartFlux, HeatFlux, scalarValues); 

    SVTable->append(&scalarValues);
    
    ////////////////// print out some statistics /////////////////////////////
    
    std::stringstream messageStream;

    messageStream << std::endl << std::endl << "Analysis | " << std::setprecision(3) 
                  << "Time : " << timing.time << " Step : " << timing.step << "   ";

    messageStream << std::setprecision(2) << std::scientific << 
                     "Field Energy : (φ) " << scalarValues.phiEnergy  << 
                     "  (A∥) " << ((Nq >= 2) ? Setup::num2str(scalarValues.ApEnergy) : "off") << 
                     "  (B∥) " << ((Nq >= 3) ? Setup::num2str(scalarValues.BpEnergy) : "off") << std::endl; 

    double total_charge = 0., kinetic_energy = 0.;
    
    for(int s = NsGlD; s <= NsGuD; s++) {
    
      messageStream << "         | "   << std::setw(10) << species[s].name << " " 
                    << std::showpos 
                    << " N : "             << scalarValues.particle_number[s-1]  
                    << " Kinetic Energy: " << scalarValues.kinetic_energy[s-1] 
                    << " Particle Flux : " << scalarValues.particle_flux[s-1]  
                    << " Heat Flux : "     << scalarValues.heat_flux[s-1] << std::endl;

      // calculate charge/kin energy for verification of the run
      total_charge   += species[s].q  * scalarValues.particle_number[s-1];
      kinetic_energy += scalarValues.kinetic_energy[s-1];

    }
      
    /////////////   Output some non-mandatory values  //////////////////////////
   
    messageStream << std::setprecision(4);
    if(vlasov->doNonLinearParallel) {

      double field_energy     = scalarValues.phiEnergy ;// skip for first+ scalarValues.ApEnergy + scalarValues.BpEnergy;
      double rel_energy_consv = 100*std::abs((kinetic_energy+field_energy)/(field_energy == 0. ? 1.e-99 : field_energy));
      messageStream << "         | Energy conservation (es) : " <<  rel_energy_consv << "%";
      //messageStream << "         | Total Energy Ratio : " <<  std::noshowpos << std::fixed << rel_energy_consv << "%";
    }

    if(Ns > 1) messageStream << "    Total Charge = " << ((species[0].n0 != 0.) ? 0. : total_charge);
    messageStream << std::endl; 

    parallel->print(messageStream.str());
    
   }
  
  }

  auxiliary->writeData(timing, dt);

}

// separately we set MPI struct
void Diagnostics::setMPIStruct()
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

}

void Diagnostics::closeData() 
{

  delete FA_HeatFluxKy; 
  delete FA_PartFluxKy;
  delete FA_CrossPhase;
  
  delete FA_grow_x; delete FA_grow_y; delete FA_grow_t;     
  delete FA_freq_x; delete FA_freq_y; delete FA_freq_t;    
     
  delete FA_XDep_Mom; delete FA_XDep_Time;
       
  delete FA_Mom_00; delete FA_Mom_20; delete FA_Mom_02;
  delete FA_Mom_10; delete FA_Mom_30; delete FA_Mom_12;
  
  delete FA_Mom_HeatFlux; delete FA_Mom_PartFlux;
  delete FA_Mom_Time;

  delete SVTable;

  H5Gclose(analysisGroup);

}
    
void Diagnostics::printOn(std::ostream &output) const
{ 
  // Something useful ? e.g. Output
}
