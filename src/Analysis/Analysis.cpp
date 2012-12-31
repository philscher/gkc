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

#include <iomanip>
#include <complex.h>

Analysis::Analysis(Parallel *_parallel, Vlasov *_vlasov, Fields *_fields, Grid *_grid, Setup *_setup, FFTSolver *_fft, FileIO *fileIO, Geometry *_geo) : 
 
parallel(_parallel),setup(_setup), vlasov(_vlasov), grid(_grid), fields(_fields), geo(_geo),  fft(_fft)

{

  initData(setup, fileIO);

  moments = new Moments(setup, vlasov, fields, grid, parallel);
}



Analysis::~Analysis() 
{
  delete moments;
  closeData();
}


void Analysis::getPowerSpectrum(CComplex  kXOut  [Nq][NzLD][Nky][FFTSolver::X_NkxL], 
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

  return;
}


//////////////////////// Calculate scalar values ///////////////////////////


void Analysis::calculateScalarValues(const CComplex f [NsLD][NmLD][NzLB][Nky][NxLB][NvLB], 
                                     const CComplex f0[NsLD][NmLD][NzLB][Nky][NxLB][NvLB],
                                     const CComplex Mom00[NsLD][NzLD][Nky][NxLD], 
                                     const CComplex Mom20[NsLD][NzLD][Nky][NxLD], 
                                     const CComplex Mom02[NsLD][NzLD][Nky][NxLD], 
                                     double ParticleFlux[Nq][NsLD][Nky][NxLD], 
                                     double HeatFlux[Nq][NsLD][Nky][NxLD],
                                     ScalarValues &scalarValues) 

{
  
  for(int s = NsGlD; s <= NsGuD; s++) {
    
      
    double number = 0.e0;
    double kineticEnergy=0.e0;
    double entropy = 0.;
    double particle = 0.; 
    double heat = 0.;
    
    if((s >= NsLlD && s <= NsLuD) && (parallel->Coord[DIR_VM] == 0))  {
    
      
    ////////////////////////////// Calculate Particle Number /////////////////////////////
    
    number =  creal(__sec_reduce_add(Mom00[s-NsLlD][0:NzLD][0][0:NxLD]));
    
    //////////// Calculate Kinetic Energy  //////////////////////////
    
    kineticEnergy = (__sec_reduce_add(creal(Mom20[s-NsLlD][0:NzLD][0][0:NxLD])) +
                     __sec_reduce_add(creal(Mom02[s-NsLlD][0:NzLD][0][0:NxLD])) ) * grid->dXYZ;
    std::cout <<   __sec_reduce_add(creal(Mom02[s-NsLlD][0:NzLD][0][0:NxLD]))  * grid->dXYZ << std::endl;
    
    kineticEnergy = 0.5 * species[s].m * (__sec_reduce_add(creal(Mom20[s-NsLlD][0:NzLD][0][0:NxLD]))) * grid->dXYZ;

    ////////////////////////////// Calculate Entropy /////////////////////////////////////
    
    //#pragma omp parallel for reduction(+:kineticEnergy) collapse (2)
    for(int m = NmLlD; m <= NmLuD; m++) {
    for(int z = NzLlD; z <= NzLuD; z++) {  for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) {
    for(int x = NxLlD; x <= NxLuD; x++) { 

      entropy += creal(pow2(__sec_reduce_add(f [s][m][z][0][x][NvLlD:NvLD])))/
                            __sec_reduce_add(f0[s][m][z][0][x][NvLlD:NvLD]);
    } } }

    } // m 
    
    /////////////////////////// Calculate Total Heat & Particle Flux ////////////////////////
       
    particle = __sec_reduce_add(ParticleFlux[Field::phi][s-NsLlD][:][:]);
    heat     = __sec_reduce_add(    HeatFlux[Field::phi][s-NsLlD][:][:]);
    
    
    } // if(local s) 

    // Communicate with other groups (make reduce whole structre)
    heat          = parallel->reduce(heat         , Op::sum);
    particle      = parallel->reduce(particle     , Op::sum);
    kineticEnergy = parallel->reduce(kineticEnergy, Op::sum);
    number        = parallel->reduce(number       , Op::sum);

    {
       scalarValues.particle_number[s-1]  = number        ;
       scalarValues.entropy        [s-1]  = entropy        ;
       scalarValues.kinetic_energy [s-1]  = kineticEnergy;
       scalarValues.particle_flux  [s-1]  = particle;
       scalarValues.heat_flux      [s-1]  = heat ;
    }  

    } // s

    return;
}


void Analysis::getParticleHeatFlux( 
                                   double ParticleFlux[Nq][NsLD][Nky][NxLD], 
                                   double HeatFlux[Nq][NsLD][Nky][NxLD],
                                   const CComplex  Field0[Nq][NzLD][Nky][NxLD],
                                   const CComplex Mom00[NsLD][NzLD][Nky][NxLD], 
                                   const CComplex Mom20[NsLD][NzLD][Nky][NxLD], 
                                   const CComplex Mom02[NsLD][NzLD][Nky][NxLD] 
                                  )
{
  // Triad condition. Heat/Particles are only transported by the y_k = 0, as the other y_k > 0 
  // modes cancels out. Thus multiplying y_k = 0 = A(y_k)*B(-y_k) = A(y_k)*[cc : B(y_k)]
  // where the complex conjugate values is used as physcial value is a real quantity.
  
  // We average over z here
  
  // Get electro-static heat/particle flux
  for(int s = NsLlD; s <= NsLuD; s++) { 
  for(int z = NzLlD; z <= NzLuD; z++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { 
  for(int x = NxLlD; x <= NxLuD; x++) { 
      
    // Do I have to include geometry terms ?!
    const CComplex iky_phi  = _imag * fft->ky(y_k) * Field0[Field::phi][z][y_k][x];

    // take only real part as it gives the radial direction ?!
    ParticleFlux[Field::phi][s-NsLlD][y_k][x-NxLlD] = - creal( iky_phi * conj(Mom00[s-NsLlD][z-NzLlD][y_k][x-NxLlD]));
        HeatFlux[Field::phi][s-NsLlD][y_k][x-NxLlD] = - creal( iky_phi * conj(Mom20[s-NsLlD][z-NzLlD][y_k][x-NxLlD]
                                                                             +Mom02[s-NsLlD][z-NzLlD][y_k][x-NxLlD]));

    } // x 
  } } // y_k, z
  }   // s

  // Get heat/particle flux from magnetic flutter

  return;
}

        
void Analysis::initData(Setup *setup, FileIO *fileIO) 
{
        
  analysisGroup = fileIO->newGroup("Analysis");
        
  hsize_t offset0[] = { 0, 0, 0, 0, 0, 0, 0 };
    
  //----------------  Moments - Heat fluxes ------------------
     
  hid_t momGroup = check(H5Gcreate(fileIO->getFileID(), "/Moments", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), DMESG("Error creating group file for Phi : H5Gcreate"));

  hsize_t mom_dsdim[] =  { Nz  , Nky, Nx  , Ns  , 1             };
  hsize_t mom_dmdim[] =  { Nz  , Nky, Nx  , Ns  , H5S_UNLIMITED };
  hsize_t mom_cBdim[] =  { NzLD, Nky, NxLD, NsLD, 1             };
  hsize_t mom_cDdim[] =  { NzLD, Nky, NxLD, NsLD, 1             };
  hsize_t mom_cmoff[] =  { 0   , 0  , 0   , 0   , 0             };
  hsize_t mom_cdoff[] =  { NzLlD-3, 0 , NxLlD-3, NsLlD-1,  0    };
     
  bool momWrite = (parallel->Coord[DIR_VM] == 0);
     
  FA_Mom_Tp        = new FileAttr("Temperature_v", momGroup, fileIO->file, 5, mom_dsdim, mom_dmdim, mom_cDdim, mom_cmoff, mom_cBdim, mom_cdoff, momWrite, fileIO->complex_tid);
  FA_Mom_To        = new FileAttr("Temperature_o", momGroup, fileIO->file, 5, mom_dsdim, mom_dmdim, mom_cDdim, mom_cmoff, mom_cBdim, mom_cdoff, momWrite, fileIO->complex_tid);
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
     
  FA_XDep_Tp   = new FileAttr("Tp", XDepGroup, fileIO->file, 3, XDep_dsdim, XDep_dmdim, XDep_cDdim, XDep_cmoff, XDep_cBdim, XDep_cdoff, XDepWrite);
  FA_XDep_To   = new FileAttr("To", XDepGroup, fileIO->file, 3, XDep_dsdim, XDep_dmdim, XDep_cDdim, XDep_cmoff, XDep_cBdim, XDep_cdoff, XDepWrite);
  FA_XDep_n    = new FileAttr("n" , XDepGroup, fileIO->file, 3, XDep_dsdim, XDep_dmdim, XDep_cDdim, XDep_cmoff, XDep_cBdim, XDep_cdoff, XDepWrite);
  FA_XDep_Time = fileIO->newTiming(XDepGroup);
        
  H5Gclose(XDepGroup);
     
  dataOutputXDep = Timing(setup->get("DataOutput.XDep.Step", -1), setup->get("DataOutput.XDep.Time", -1.));
  
  //---------------  Analysis - Heat fluxes ------------
     
  // Heat Flux ky and Particle FluxKy ( per species) 
  hid_t fluxGroup = fileIO->newGroup("Flux", analysisGroup);
     
  hsize_t FSky_dsdim[] = { Nq, Ns     , Nky   , Nx     , 1 }; 
  hsize_t FSky_dmdim[] = { Nq, Ns     , Nky   , Nx     , H5S_UNLIMITED} ;
  hsize_t FSky_cDdim[] = { Nq, NsLD   , Nky   , NxLD   , 1 };
  hsize_t FSky_cBdim[] = { Nq, NsLD   , Nky   , NxLD   , 1 };
  hsize_t FSky_cdoff[] = { 0,  NsLlB-1, NkyLlD, NxLlD-3, 0};

  bool isXYZ = parallel->Coord[DIR_XYZ] == 0;

  FA_heatKy      = new FileAttr("Heat"    , fluxGroup, fileIO->file, 5, FSky_dsdim, FSky_dmdim, FSky_cDdim, offset0,  FSky_cBdim, FSky_cdoff, isXYZ);
  FA_particleKy  = new FileAttr("Particle", fluxGroup, fileIO->file, 5, FSky_dsdim, FSky_dmdim, FSky_cDdim, offset0,  FSky_cBdim, FSky_cdoff, isXYZ);
    
  H5Gclose(fluxGroup);
     

  //---------------------------   Power Spectrum  -------------------------------
   
  bool isMaster = parallel->myRank == 0;
  // X-scalarValue
  hid_t growGroup = fileIO->newGroup("PowerSpectrum", analysisGroup);

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

  return;

}


void Analysis::getFieldEnergy(double& phiEnergy, double& ApEnergy, double& BpEnergy)
{

  fields->getFieldEnergy(phiEnergy, ApEnergy, BpEnergy);

  return;
}
  
void Analysis::writeData(const Timing &timing, const double dt)

{


  if (timing.check(dataOutputMoments   , dt) || timing.check(dataOutputXDep, dt) ||
      timing.check(dataOutputStatistics, dt)) 
  {

      CComplex  Mom_Tp[NsLD][NzLD][Nky][NxLD], 
                Mom_To[NsLD][NzLD][Nky][NxLD], 
                Mom_n [NsLD][NzLD][Nky][NxLD]; 
        
      double     HeatFlux[Nq][NsLD][Nky][NxLD],  // Heat flux   
             ParticleFlux[Nq][NsLD][Nky][NxLD];  // Particle flux


      // Get Moments of Vlasov equation
      moments->getMoments((A6zz) vlasov->f, (A4zz) fields->Field0, Mom_n, Mom_Tp, Mom_To);
        
      // only electro-static flux is calculated
      getParticleHeatFlux(ParticleFlux, HeatFlux, (A4zz) fields->Field0, Mom_n, Mom_Tp, Mom_To);

      ////////////////// Output Moments /////////////////////

      if (timing.check(dataOutputMoments, dt)       )   {

        FA_Mom_Density->write((CComplex *) Mom_n );
        FA_Mom_Tp     ->write((CComplex *) Mom_Tp);
        FA_Mom_To     ->write((CComplex *) Mom_To);
        FA_Mom_Time->write(&timing);

        parallel->print("Data I/O : Moments output");
      }

      ////////////////// Store X-dependent data /////////////
      if (timing.check(dataOutputXDep, dt)       )   {
      
        double   A1_Tp[NsLD][NxLD], A1_To[NsLD][NxLD], 
                 A1_n [NsLD][NxLD],  A1_p[NsLD][NxLD];

        // Reduce over y_k and z
        for(int s = NsLlD; s <= NsLuD; s++) {  for(int x = NxLlD;  x <= NxLuD; x++) {

         A1_n [s-NsLlD][x-NxLlD] = __sec_reduce_add(creal(Mom_n [s-NsLlD][:][0][x-NxLlD]));
         A1_Tp[s-NsLlD][x-NxLlD] = __sec_reduce_add(creal(Mom_Tp[s-NsLlD][:][0][x-NxLlD]));
         A1_To[s-NsLlD][x-NxLlD] = __sec_reduce_add(creal(Mom_To[s-NsLlD][:][0][x-NxLlD]));

        } } 
   
        FA_XDep_Tp  ->write((double *) A1_Tp); 
        FA_XDep_To  ->write((double *) A1_To); 
        FA_XDep_n   ->write((double *) A1_n ); 
        FA_XDep_Time->write(&timing);

        parallel->print("Data I/O : X-Dep output");


        FA_particleKy->write((double *) ParticleFlux);
        FA_heatKy    ->write((double *)     HeatFlux);
      }  

    if (timing.check(dataOutputStatistics, dt)       )   {
 
      ////////////// Scalar Variables /////////////////

      { // calculate mode spectrum of fields (phi, Ap, Bp)

        double pSpecX [Nq][Nx], pSpecY [Nq][Nky],
               pPhaseX[Nq][Nx], pPhaseY[Nq][Nky];
     
        getPowerSpectrum((A4zz) fft->kXOut, (A4zz) fields->Field0, pSpecX, pSpecY, pPhaseX, pPhaseY);
    
        // Seperatly writing ? Hopefully it is buffered ... (passing stack pointer ... OK ?)
        FA_grow_x->write( &pSpecX [0][0]); FA_grow_y->write(&pSpecY [0][0]); FA_grow_t->write(&timing);
        FA_freq_x->write( &pPhaseX[0][0]); FA_freq_y->write(&pPhaseY[0][0]); FA_freq_t->write(&timing);
      }

      ScalarValues scalarValues;
    
      // calculate kinetic Energy first, need for initial_e ! sum over domain
      scalarValues.timestep = timing.step;
      scalarValues.time     = timing.time;
    
      fields->getFieldEnergy(scalarValues.phiEnergy, scalarValues.ApEnergy, scalarValues.BpEnergy);
    
      //  Get scalar Values for every species ( this is bad calculate them alltogether)
      calculateScalarValues((A6zz) vlasov->f, (A6zz) vlasov->f0,
                             Mom_n, Mom_Tp, Mom_To, ParticleFlux, HeatFlux,
                             scalarValues); 

      SVTable->append(&scalarValues);
    
      ////////////////// print out some statistics /////////////////////////////
    
      std::stringstream messageStream;
      messageStream << std::endl << std::endl << "Analysis | " << std::setprecision(3) << "Time : " << timing.time << " Step : " << timing.step << "   ";
      messageStream << std::setprecision(2) << std::scientific << 
                       "Field Energy : (φ) " << scalarValues.phiEnergy  << 
                       "  (A∥) " << ((Nq >= 2) ? Setup::num2str(scalarValues.ApEnergy) : "off") << 
                       "  (B∥) " << ((Nq >= 3) ? Setup::num2str(scalarValues.BpEnergy) : "off") << std::endl; 
      double charge = 0., kinetic_energy=0.;
    
      for(int s = NsGlD; s <= NsGuD; s++) {
    
        messageStream << "         | "   << std::setw(10) << species[s].name << " " 
//                      << std::setprecision(2) << std::scientific << std::showpos 
                      << std::showpos 
                      << " N : "             << scalarValues.particle_number[s-1]  
                      << " Kinetic Energy: " << scalarValues.kinetic_energy[s-1] 
                      << " Particle Flux : " << scalarValues.particle_flux[s-1]  
                      << " Heat Flux : "     << scalarValues.heat_flux[s-1] << std::endl;
        charge += species[s].q  * scalarValues.particle_number[s-1];
       kinetic_energy += scalarValues.kinetic_energy[s-1];

      }
      
      /////////////   Output some non-mandatory values  //////////////////////////
      const double field_energy = scalarValues.phiEnergy + scalarValues.ApEnergy + scalarValues.BpEnergy;
   
      messageStream << std::setprecision(4);
      if(vlasov->doNonLinearParallel) messageStream << "         | Total Energy Ratio : " <<  std::noshowpos << std::fixed << 100*abs((kinetic_energy+field_energy)/(field_energy == 0. ? 1.e-99 : field_energy)) << "%";
      if(Ns > 1                     ) messageStream << "    Total Charge = " << ((species[0].n0 != 0.) ? 0. : charge);
      messageStream << std::endl; 
      messageStream << std::setprecision(5) << field_energy/(kinetic_energy == 0. ? 1.e-99 : kinetic_energy) << " " << kinetic_energy/(field_energy == 0. ? 1.e-99 : field_energy) << std::endl;

      parallel->print(messageStream.str());
    
     }
  
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

}



void Analysis::closeData() 
{

  delete FA_heatKy; 
  delete FA_particleKy;
  delete FA_grow_x; delete FA_grow_y; delete FA_grow_t;     
  delete FA_freq_x; delete FA_freq_y; delete FA_freq_t;    
     
  delete FA_XDep_Tp; delete FA_XDep_To; 
  delete FA_XDep_n ; delete FA_XDep_Time;
       
  delete FA_Mom_Tp;
  delete FA_Mom_To;
  delete FA_Mom_HeatFlux;
  delete FA_Mom_Density;
  delete FA_Mom_Time;

  delete SVTable;

  H5Gclose(analysisGroup);

  return;
}
    
void Analysis::printOn(std::ostream &output) const
{ 


  return;
}
