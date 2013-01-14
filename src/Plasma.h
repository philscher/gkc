/*
 * =====================================================================================
 *
 *       Filename: Plasma.h
 *
 *    Description: Properties of plasma
 *
 *         Author: Paul P. Hilscher (2009-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __GKC_PLASMA_H__
#define __GKC_PLASMA_H__

#include "Global.h"
#include "Setup.h"
#include "FileIO.h"
#include "Grid.h"
#include "Geometry.h"



// Maximum number of species (excluding adiabatic species) to be 
// included without recompilation
#define SPECIES_MAX 4

/**
*   @brief Information about Species
*
*   Definition of species included.
*
**/
struct Species 
{

  Species() : q(0.), m(0.), doGyro(true), T0(0.), n0(0.) 
  { 
    alpha     = 0.;
    sigma     = 0.;
    gyroModel = "";
      
//    ArrayProfiles = nct::allocate(nct::Range(NxGlB-2, NxGB+4))(&T, &n, &w_T, &w_n, &src, &krook);
  };
   
  double q;          ///< Charge 
  double m;          ///< Mass 
  double T0;         ///< Temperature normalization
  double n0;         ///< Density normalization

  bool   doGyro;     ///< Set if gyro-averaging is performed
   
  double v_th;       ///< Velocity scale / Thermal velocity
  double sigma;      ///< sigma 
  double alpha;      ///< alpha

  char name[64];     ///< name of species
  char n_name[64];   ///< Density function string    , e.g  n(x) = "1"
  char T_name[64];   ///< Temperature function string, e.g. T(x) = "1/x"
  
  double T[1024]    , ///< Temperature profile
         n[1024]    , ///< Density profile
         w_T[1024]  , ///< Temperature scale length
         w_n[1024]  , ///< Density scale length
         src[1024]  , ///< (Energy source term)
         krook[1024]; ///< Krook operator used for damping
  
  //std::string spec_name, ///< name of species
  //            func_dens, ///< Density function string    , e.g  n(x) = "1"
  //            func_Temp; ///< Temperature function string, e.g. T(x) = "1/x"
  
  /////////// Below non-POD types (and not saved [yet] in HDF-5 file) ////////////
  std::string gyroModel;
  std::string f0_str;
  nct::allocate ArrayProfiles;
      
  // stupid fix, but we have to otherwise all stuff is private
  void update(Geometry *geo, double cs) { 

    v_th = sqrt(2.*T0/m); 
    alpha   = v_th*  1./(cs*sqrt(geo->eps_hat));
    sigma   = q / T0;
   
  };

  /// Calculate debye legnth
  double debye2(const int x) const { return T[x]/(4.*M_PI*n[x]*q*q); };
   
};

/**
*   @brief Hold information about the plasma with species and normalizations
*
*   General normalization and intialization of species, including 
*   initial temperature profiles / density profiles etc.
*
*
**/
class Plasma : public IfaceGKC {

 public:

  double n_ref, ///< Reference 
         L_ref, ///< Reference scale length
       rho_ref, ///< Gyro-radius reference length
         c_ref, ///< Sound speed reference
         T_ref; ///< Reference temperature
         
  double rho_star; ///< Sound speed reference
   
  double cs;  ///<  speed of sound of ions \f$ c_s = \sqrt{\frac{T_{e0}}{m_i}} \f$ 

    
  /**  The normalized Debye length (from e.g. G\"orler PhD Eq. (2.82))
  *
  * \f[
  *      \lambda_D = \lambda_D / \rho_{ref} = \sqrt{T_{ref}}{4\pi\rho_{ref}^2 n_{ref} e^2}
  * \f]
  *  @note needs to be 1-d over X
  **/
  double debye2;

  /**
  *   @todo is it necessary ?
  *
  *
  **/
  bool global;

  double B0,    ///< Scale of Equilibrium Magnetic field \f$ B_0 \f$
         beta,  ///< Plasma pressure  \f$ \beta \f$  
         w_p;   ///< Plasma pressure scale length \f$ \omega_p \f$

  /**
  *   @todo is it necessary ? It is defined in Fields.cpp as Nq anyway, ..
  *
  **/
  int nfields;
   
  /**
  *    @brief intializes species 
  * 
  **/
  Plasma(Setup *setup, FileIO *fileIO, Geometry *geo, const int nfields=1);
   
  /**
  *
  *  @brief Free up resources
  *
  **/
 ~Plasma();


 protected:

  virtual void printOn(std::ostream &output) const;

  void initData(FileIO *fileIO) ;
  virtual void writeData(const Timing &timing, const double dt) {};
  virtual void closeData() {};

};


#endif // __GKC_PLASMA_H__

