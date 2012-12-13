/*
 * =====================================================================================
 *
 *       Filename: Analysis.h
 *
 *    Description: Analysis functions for gkc
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef ANALYSIS_H__
#define ANALYSIS_H__

#include "Global.h"

#include "Setup.h"
#include "Parallel/Parallel.h"
#include "Grid.h"
#include "Vlasov/Vlasov.h"
#include "FFTSolver/FFTSolver.h"
#include "Geometry/Geometry.h"
#include "Timing.h"
#include "Plasma.h"

/**
*    @brief Data analysis and output
*
*
*    @note This kitchen is a fucking mess [(c) Ramsay] !!! Clean up & polish
*
*
*    Parallelization, 
*
*    @todo Writting operations are done step-by-step, is HDF-5 buffering them and
*          perform all at once, do we have to set some options, or does it not
*          buffer at all ?
*
*    @todo Write MPI reduce operations
*
**/
class Analysis : public IfaceGKC {

  hid_t analysisGroup; ///< HDF-5 id for group label

  /**
  *    @brief Structure for storing scalar values diagnostics
  *
  *    Structure represents scalar values of quantities calculated over
  *    the whole simulation domain.
  *
  **/ 
  typedef struct ScalarValues_t
  { 
    int    timestep;                        ///< Current time step 
    double time;                            ///< Current time
    double phiEnergy;                       ///< Electric field energy
    double ApEnergy;                        ///< Magnetic field energy (from \f$ A_\parallel \f$)
    double BpEnergy;                        ///< Magnetic field energy (from \f$ B_\parallel \f$)
    double particle_number[SPECIES_MAX];    ///< Total particle number (per species)
    double kinetic_energy [SPECIES_MAX];    ///< Total kinetic energy (per species)
    double entropy        [SPECIES_MAX];    ///< Total entropy (per species)
    double heat_flux      [SPECIES_MAX];    ///< Total heat flux (per species)
    double particle_flux  [SPECIES_MAX];    ///< Total particle flux (per species)

  } ScalarValues;
   
  //////////////  Data Output Stuff //////////////////
   
  Timing dataOutputStatistics, ///< Timing to define output of scalarValues
         dataOutputMoments   , ///< Timing to define output of moments
         dataOutputXDep      ; ///< Timing to define output of X-dependent variables

  ///@{
  ///@inb group HDF-5 Attributes 
  FileAttr *FA_Mom_Density ,  ///< Density  
           *FA_Mom_Tp      ,  ///< Temperature 
           *FA_Mom_HeatFlux,  ///< Heat flux
           *FA_Mom_Time    ;  ///< Time 

  FileAttr  *FA_heatKy,     ///< Heat Flux as  Q(s,x,ky)
            *FA_particleKy; ///< Particle Flux as X(s,x,ky)

  FileAttr  *FA_grow_x,     ///< Mode-Power (kx)
            *FA_grow_y,     ///< Mode-Power (ky) 
            *FA_grow_t;     ///< Mode-Power Time

  FileAttr  *FA_freq_x,     ///< Phase (kx)
            *FA_freq_y,     ///< Phase (ky)
            *FA_freq_t;     ///< Time

     
  FileAttr *FA_XDep_Tp,     ///< X-Dependence (Temperature)
           *FA_XDep_n,      ///< X-Dependence (Density)
           *FA_XDep_Time;   ///< X-Depdndence (Time)

  TableAttr *SVTable;       ///< Table for scalar Values

  ///@}

  //////////////////////////////////////////////////////////////
  Parallel *parallel;
  Setup *setup;
  Vlasov *vlasov;
  Grid   *grid;
  Fields *fields;
  FFTSolver *fft;   ///< required for heat flux calculations
  Geometry *geo;


  // MPI structures
  void setMPIStruct();

 public:

  Analysis(Parallel *_parallel, Vlasov *vlasov, Fields *_fields, Grid *_grid, Setup *_setup, FFTSolver *fft, FileIO *fileIO, Geometry *geometry);

 ~Analysis() ;


  /** \brief Calculate the Power Spectrum of the Phi/B potential 
  *
  *  NOTE : Only r2c format is supported. Also we require that the FFT library
  *         has r2c in x-y and only the half-format spectrum in z. If other
  *         decomposition are required very minor modification need to be done
  *         for this function.
  *
  *         Calculates the power of each mode and corresponding phase (better do in
  *         post processing ?)
  *
  *         Power \f[ abs( F(n, x,k_y) )^2 \f]
  *         Phase \f[ arctan2( F(n,x,k_y)) \f]
  *
  *         where $n=\phi,A_{1\parallel}, B_{1\parallel}$.
  *
  *
  *  @return double Array of Power density spectra in direction dir
  *
  **/
  void getPowerSpectrum(CComplex  kXOut[Nq][NzLD][NkyLD][FFTSolver::X_NkxL], 
                        CComplex Field0[Nq][NzLD][NkyLD][NxLD]             , 
                        double   pSepcX [plasma->nfields][Nx/2]            , 
                        double   pSpecY [plasma->nfields][Nky ]            ,
                        double   pPhaseX[plasma->nfields][Nx/2]            , 
                        double   pPhaseY[plasma->nfields][Nky ]           );


    
  /** 
  *
  *  @brief Calculate the kinetic Energy for particles
  *  
  *  If TOTAL or nothing is provided, the total kinetic energy of
  *  all particles with be given. It calculates according to
  *
  *  Physical :
  *
  *  \f[
  *      E_{\textrm{kin}} = int_{v_parallel=-infty}^infty \int_{\mu=0}i^infty
  *                         left( frac{m_sigma}{2} v_parallel^2 + \mu B_0 / hat{T} right) f_sigma d^6Z
  *
  *  \f]
  *
  *  Normalized :
  *
  *  \f[
  *      E_{\textrm{kin}} = int_{v_parallel=-infty}^infty \int_{\mu=0}i^infty
  *                         left( hat{v}^2 + \mu B_0 / hat{T} right) f_sigma d^6Z
  *  \f]
  *
  *  Normalized : Note that \f$ d6Z = m\hat{B} d\mu dv_\parallel \f$
  *
  *
  *  Calculate entropy using Imadera et al.
  *
  *  @param species particle species \f$ \sigma \f$
  *  @return   the total energy of species
  *
  **/
  void calculateScalarValues(const CComplex f [NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB], 
                              const CComplex f0[NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB], 
                              const double V[NvGB], const double M[NmGB], 
                              ScalarValues &SV); 


  /**
  *  @brief calculates and returns the heat flux across magnetic flux surfaces
  *         for species \f$ \sigma \f$
  *
  *  Obtained from Lapillione, PhD Thesis (2010).
  *
  *  For the particle flux Eq. (2.118)
  *
  *  \f[
  *     \Gamma_\sigma = \left< \int F_{\sigma j} v_{Ex} d^3v  \right> 
  *  \f]
  *
  *
  *  and the heat flux Eq.(2.123)
  *  
  *  \f[
  *     Q_\sigma = -\frac{1}{2C} \left< 
  *                \int \left(m_\sigma v_\parallel^2 + 2 \mu B \right) F_{1\sigma}
  *                \frac{\partial \bar{\phi}}{\partial y} 
  *                \right> 
  *  \f]
  *
  *  @todo include FLR correction terms
  *        include Geometry factors
  * 
  *
  **/
  void getParticleHeatFlux(const int m, const int s, 
                           double ParticleFlux[NkyLD][NxLD], double HeatFlux[NkyLD][NxLD],
                           const CComplex         f[NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const CComplex Field[Nq][NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           const double V[NvGB], const double M[NmGB]);




  //////////////// Calculate Moments (4-dimensional s,x,y_k,z) ////////////////////////
  
  /** 
  *
  *  @brief Calculates the temperature, note the 
  * 
  *  Calculates the temperature by
  *
  *  \f[ T_\sigma = \frac{ \int_{-infty}{infty} 2 v^2 f}{n_\sigma} \f]
  *
  *
  *  \f[ \int_{-\infty}{\infty} v^2 \exp^{-\frac{v^2}{T}} = \frac{1}{2} \sqrt{\pi}T^{3/2} \f]
  *  so we need to multiply with factor 1./2. Also divide by n.
  *
  *  @param f    the perturbed phase space function \f$ f(s,m,z,y_k, x,v) \f$
  *  @param A4_z the number density in \f$ n4(\sigma,z,k_y,x) \f$
  *
  **/
  void getTemperature(const CComplex f[NsLB][NmLB][NzLB][NkyLD][NxLB][NvLB],
                            CComplex T[NsLD][NzLD][NkyLD][NxLD]); 


  void getNumberDensity(const CComplex f[NsLB][NmLB][NzLB][NkyLD][NxLB][NvLB],
                              CComplex n[NsLD][NzLD][NkyLD][NxLD]) ;


  //////////////////////////////////  Data-I/0 stuff ////////////////////////////////////////

  void initData(Setup *setup, FileIO *fileIO) ;
  void writeData(const Timing &timing, const double dt);
  void closeData();
      
  void printOn(std::ostream &output) const;

};



#endif // ANALYSIS_H__
