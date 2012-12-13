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

#ifndef PLASMA_H__
#define PLASMA_H__

#include "Global.h"
#include "Setup.h"
#include "FileIO.h"
#include "Grid.h"
#include "Geometry.h"


// BUG : if SPECIES_MAX >= 8 HDF-5 gives bus error
// why ?1

// Excluding adiabatic species
#define SPECIES_MAX 4
/**
*   @brief Hold information about the plasma with species and normalizations
*
**/
class Plasma : public IfaceGKC {

  public:

   /**
   *   @brief Information about Species
   *
   **/
   typedef struct Species 
   {

     Species() : q(0.), m(0.), collision(0.), w_n(0.), w_T(0.), doGyro(true), T0(0.), n0(0.) 
     { 
       alpha   = 0.;
       sigma   = 0.;
       nct::allocate(nct::Range(NxLlB, NxLB))(&T, &n);
       gyroModel = "";
     };
   
     double q;          ///< Charge 
     double m;          ///< Mass 
     double collision;  ///< Collisional Frequency 
     double w_n;        ///< Density gradient
     double w_T;        ///< Temperature gradient
     double T0;         ///< Temperature normalization
     double n0;         ///< Density normalization
     bool   doGyro;     ///< Set if gyro-averaging is performed
   
     double scale_v;    ///< Velocity scale / Thermal velocity
     double scale_n;    ///< Density scale 
     double sigma;      ///< sigma 
     double alpha;      ///< alpha
   
   
     char name[64];     ///< name of species
     char n_name[64];   ///< dont know
     char T_name[64];   ///< dont know
  
     /////////// Below non-POD types (and not saved [yet] in HDF-5 file) ////////////
     std::string gyroModel;
     std::string f0_str;
      
     // stupid fix, but we have to otherwise all stuff is private
   
   void update(Geometry *geo, double cs) { 
        scale_v = sqrt(2.*T0/m); 
        scale_n = n0;
        alpha = scale_v*  1./(cs*sqrt(geo->eps_hat));
        sigma = q / T0;
   
     };

     /// Calculate debye legnth
     double debye2(const int x) const { return T[x]/(4.*M_PI*n[x]*q*q); };
   
     double *T, ///< Temperature profile
            *n; ///< Density profile
   
   } _Species;

   // Check what is really necessary also in plasma

   double n_ref, ///< Reference 
          L_ref, ///< Reference scale length
          T_ref; ///< Reference temperature
   
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
   **/
   bool global;

   double B0,    ///< Scale of Equilibrium Magnetic field \f$ B_0 \f$
          beta,  ///< Plasma pressure  \f$ \beta \f$  
          w_p;   ///< Plasma pressure scale length \f$ \omega_p \f$

   /**
   *   @todo is it necessary ?
   **/
   int nfields;
   
   /**
   *   Note Species goes from 0 ... SPECIES_MAX, where 0 is an
   *   adiabatic species.
   *   @todo check impact on speed
   **/
   Species species[SPECIES_MAX+1];
   
   /**
   *    @brief intializes species 
   *
   *
   * 
   **/
   Plasma(Setup *setup, FileIO *fileIO, Geometry *geo, const int nfields=1);
   
   /**
    *
    *  @brief Free up resources
    *
   **/
   virtual ~Plasma() {};


protected:

    virtual void printOn(std::ostream &output) const;

    void initData(FileIO *fileIO) ;
    virtual void writeData(const Timing &timing, const double dt) {};
    virtual void closeData() {};
};


#endif // CONFIG_H__

