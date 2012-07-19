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
#include "config.h"
#include "Parallel.h"
#include "Grid.h"
#include "Vlasov.h"
#include "FFTSolver.h"
#include "Geometry.h"
#include "GeometryShear.h"
#include "GeometrySlab.h"
#include "Geometry2D.h"
#include "Timing.h"
#include "Plasma.h"

#define TOTAL 0


/**
*    @brief Data analysis and output
*
*
*
*
**/
class Analysis : public IfaceHelios {

   hid_t analysisGroup;
   /// scalarValue - Structure to store scalar data.
   typedef struct ScalarValues
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
   } _ScalarValues;

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
   FFTSolver *fft;
   Geometry<HELIOS_GEOMETRY> *geo;
   Array1d  initialEkin, dT;
   Array2c spectrumXZ, spectrumYZ, spectrumXY;
   Array3d pSpec, pPhase;
   Array3z pFreq;
   Array2d heatFluxKy;
   Array3z A_xyz;
   Array3z phi_k;
   Array4d A4;
   Array4z A4_z;

   template<typename tV, typename tM> cmplxd  F (int  x, int y, int z, tV v, tM m, int s) {

      if (plasma->global == false)               return sum(vlasov->f(x,y,z,v,m,s) + vlasov->f0(x,y,z,v,m,s));
      else                                       return sum(vlasov->f(x,y,z,v,m,s));

   }

   template<typename tV, typename tM> cmplxd  dF(int  x, int y, int z, tV v, tM m, int s) {

    if (plasma->global == false) return sum(vlasov->f(x,y,z,v,m,s));
    else                         return sum(vlasov->f(x,y,z,v,m,s) - vlasov->f0(x,y,z,v,m,s));

   }
 
   Timing dataOutputStatistics, dataOutputMoments;

  public:

   Analysis(Parallel *_parallel, Vlasov *vlasov, Fields *_fields, Grid *_grid, Setup *_setup, FFTSolver *fft, FileIO *fileIO, Geometry<HELIOS_GEOMETRY> *geometry);

   ~Analysis() ;


   /** \brief Calculate the Power Spectrum of the Phi/B potential 
   *
   *  NOTE : Only r2c format is supported. Also we require that the FFT library
   *         has r2c in x-y and only the half-format spectrum in z. If other
   *         decomposition are required very minor modification need to be done
   *         for this function.
   *
   *  @param dir Direction either DIR_X, DIR_Y or DIR_Z
   *
   *  @return double Array of Power density spectra in direction dir
   *
   **/
   Array3d getPowerSpectrum();
    
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
   * @param species particle species \f$ \sigma \f$
   *
   * @return   the total energy of species
   **/
   double getKineticEnergy (int species=TOTAL);
   
   double getEntropy       (int species=TOTAL);
   
   double getParticelNumber(int species=TOTAL);



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
   **/
   Array4z getHeatFlux      (int species=TOTAL);

   double  getTotalHeatFlux      (int species=TOTAL);
    
   Array3d getHeatFluxKy      (int species=TOTAL);
    
   /**
   *  calculates and returns the heat flux across magnetic flux surfaces
   *  for species s by
   *  
   *  \f[
   *     \Gamma_\sigma = \left< \int F_{\sigma j} v_{Ex} d^3v  \right> 
   *  \f]
   * 
   *  Ref : Dannert PhD Thesis
   **/
   double getTotalParticleFlux   ( int species=TOTAL);
    
   Array3d getParticleFluxKy ( int species=TOTAL);

   /**
   *
   *  Get electrostatic field energy, note that the field energy consist of various contributions
   *  depending on what kind of simulations are performed. This is
   *
   *  \f[
   *         \Gamma_{es} = - \left< F_{j1} k_y \phi(x,y_k,z,m,s) d^3 \right> 
   *  \f]
   *
   *  Ref : Dannert PhD Thesis
   *
   **/
   void getFieldEnergy(double& phiEnergy, double& ApEnergy, double& BpEnergy);
    
   int updateSpectrum      (unsigned int dir);
   Array2c getSpectrum(unsigned int dir);

   double getMaxTimeStep(Timing timing, int dir=DIR_ALL, const double maxCFL=0.4);

   Array4d getNumberDensity(const bool total=true);
   Array4d getMomentumParallel();
   Array4z getTemperatureParallel();
   Array4d getTemperatureOthogonal();
   Array4d getHeatFluxOrthogonal();
   Array4d getHeatFluxParallel ();

   void initDataOutput(Setup *setup, FileIO *fileIO) ;
   int  writeData(Timing timing, double dt);
   void closeData();
      
   void printOn(ostream &output) const;
};



#endif // ANALYSIS_H__
