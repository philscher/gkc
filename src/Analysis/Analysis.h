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

   hid_t analysisGroup;

   /// scalarValue - Structure to store scalar data.
   typedef struct ScalarValues_t
   { 
      int    timestep;
      double time;
      double phiEnergy;
      double ApEnergy;
      double BpEnergy;
      double particle_number[SPECIES_MAX];
      double kinetic_energy [SPECIES_MAX];
      double entropy        [SPECIES_MAX];
      double heat_flux      [SPECIES_MAX];
      double particle_flux  [SPECIES_MAX];

      // Below are variables not included in HDF-5 file
      // however also included in MPI reduce operation


   } ScalarValues;

   /// Moments Variables 
   FileAttr *FA_Mom_Density, *FA_mom_up, *FA_Mom_Tp, *FA_mom_To, *FA_mom_qp, *FA_mom_qo, * FA_mom_QES, *FA_Mom_Time, *FA_Mom_HeatFlux;
   /// Spectrum Variables   
   FileAttr *FA_spec_yz, *FA_spec_xy, *FA_spec_xz, *FA_spec_time;
   /// Data Output Stuff
   FileAttr  *FA_heatKy, *FA_particleKy;
   FileAttr  *FA_grow_x, *FA_grow_y, *FA_grow_z, *FA_grow_t;
   FileAttr  *FA_freq_x, *FA_freq_y, *FA_freq_z, *FA_freq_t;
   TableAttr *SVTable;

   Parallel *parallel;
   Setup *setup;
   Vlasov *vlasov;
   Grid   *grid;
   Fields *fields;
   FFTSolver *fft;   ///< required for heat flux calculations
   Geometry *geo;


   CComplex *A4;
   nct::allocate ArrayA4;
   
   Timing dataOutputStatistics, ///< Timing to define output of scalarValues
          dataOutputMoments;    ///< Timing to define output of moments

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


    
   /** \brief Calculate the kinetic Energy for particles
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
   *   Calculate entropy using Imadera et al.
   *
   *
   *  @param species particle species \f$ \sigma \f$
   *  @return   the total energy of species
   *
   *
   **/
   void calculateScalarValues(const CComplex f [NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB], 
                              const CComplex f0[NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB], 
                              const double V[NvGB], const double M[NmGB], const int s,
                              ScalarValues &SV); 


   /** Calculates the temperature, note the 
   * 
   * Calculates the temperature by
   *
   *  \f[ T_\sigma = \frac{ \int_{-infty}{infty} 2 v^2 f}{n_\sigma} \f]
   *
   *
   *  \f[ \int_{-\infty}{\infty} v^2 \exp^{-\frac{v^2}{T}} = \frac{1}{2} \sqrt{\pi}T^{3/2} \f]
   *  so we need to multiply with factor 1./2. Also divide by n.
   *
   *
   **/

   /**
   *  calculates and returns the heat flux across magnetic flux surfaces
   *  for species s by \f[ E \times B_0 \f]
   *  
   *  \f[
   *     Q_s = \left< \int \left(v_\parallel^2 + \mu B  \right) F_sj k_y \phi(x,k_y, z, \mu, s) d^3v \right> 
   *  \f]
   *
   *  calculates and returns the heat flux across magnetic flux surfaces
   *  for species s by
   *  
   *  \f[
   *     \Gamma_\sigma = \left< \int F_{\sigma j} v_{Ex} d^3v  \right> 
   *  \f]
   * 
   *  Ref : Dannert PhD Thesis
   *
   **/
   void getParticleHeatFlux(const int m, const int s, 
                            CComplex ParticleFlux[NkyLD][NxLD], CComplex HeatFlux[NkyLD][NxLD],
                            const CComplex f[NsLB][NmLB][NzLB][NkyLD][NxLB][NvLB],
                            const CComplex phi[NsLD][NmLD][NzLD][NkyLD][NxLB+4],
                            const double V[NvGB], const double M[NmGB]);


   int updateSpectrum      (unsigned int dir);
   
   /**
   *    Get the field energy
   *
   *    @note internally it is just deliagting it to the
   *    Fields class.
   *
   *    @todo find out if this is a bad habit for programm design.
   *          I would guess so.
   *
   **/ 
   void getFieldEnergy(double& phiEnergy, double& ApEnergy, double& BpEnergy);

   //////////////// Calculate Moments (4-dimensional s,x,y_k,z) ////////////////////////

   /**
   *   What for is this exactly ?
   *
   **/
   
   void getTemperatureParallel(const  CComplex f[NsLB][NmLB][NzLB][NkyLD][NxLB][NvLB],
                                      CComplex A4_z[NsLD][NzLD][NkyLD][NxLD], 
                               const double V[NvGB], const double M[NmGB]);
   // depreciated
   void getTemperatureOthogonal();
   void getHeatFluxOrthogonal();
   void getHeatFluxParallel ();
   
   void getNumberDensity(const  CComplex f[NsLB][NmLB][NzLB][NkyLD][NxLB][NvLB],
                                CComplex D[NsLD][NzLD][NkyLD][NxLD], 
                         const double V[NvGB], const double M[NmGB]);
   void getMomentumParallel();


   //////////// Calculate X-Dependent variables ( 2-dimensional species & X) /////////////





   //////////////////////////////////  Data-I/0 stuff ////////////////////////////////////////

   void initDataOutput(Setup *setup, FileIO *fileIO) ;
   void writeData(Timing timing, double dt);
   void closeData();
      
   void printOn(std::ostream &output) const;
};



#endif // ANALYSIS_H__
