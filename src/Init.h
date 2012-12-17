/*
 * =====================================================================================
 *
 *       Filename: Init.h
 *
 *    Description: Definition of Initial conditions 
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#ifndef __INIT_H_
#define __INIT_H_

#include "Global.h"



#include "Grid.h"
#include "Vlasov/Vlasov.h"
#include "Geometry/Geometry.h"
#include "Parallel/Parallel.h"


/**
*   @brief Initialized the background field and intial perturbation
*
*
*
**/ 
class Init   : public IfaceGKC 
{

   Geometry *geo;
   
   double epsilon_0, ///< Strengh of perturbation
          sigma;     ///< Additional value e.g. FWHM of exp

   /**
   *
   *   @brief Defines seed to be used for random number generator (RNG)
   *  
   *   If set to zero, a combination of the system's time and 
   *   process id will be used, with the effect that due to
   *   the initial conditions, simulations are not reproduction
   *   unless seed is same. 
   *
   **/
   int random_seed; 
   

   /**
   *   @brief inits the Maxwellian 
   **/
   void initBackground(Setup *setup, Grid *grid, 
                       CComplex f0[NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB],
                       CComplex f [NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB]);

  public :

   /**
   *   @brief document me please
   **/
   std::string PerturbationMethod;
   /**
   *   @brief document me please
   **/
   Init(Parallel *parallel, Grid *grid, Setup *setup, FileIO *fileIO, Vlasov *vlasov, Fields *fields, Geometry *geo);
   /**
   *   @brief document me please
   **/
   ~Init() {};

   /**
    *  @brief Function to set Initital conditions
    *
    *  This function is used to set the initial conditions for the
    *  gyrokinetic simulation. For the electrostatic simulation
    *  one only need to set the temperature profile (T), its
    *  derivative Tx and the inital Phasespace function f0
    *  and the perturbed one f
    *
    *  T  = ...
    *  f0 = ....
    *  f  = f0
    *
    *  The phase space function takes the function of
    *  \f[ 
    *      \hat{F}_0(v_\parallel, \mu) = \pi^{3/2} e^{-(v_\parallel^2 + \mu \hat{B}}
    *  \f]
    * where a scaling is applied according to
    * \f[
    *      x rho_i = \hat{x} \\ 
    *      \hat{\phi} = \frac{e}{\hat{T}_{ref}}\phi
    *
    * \f]
    *
    *  thus for the dependent variables
    * \f[
    *  \hat{F}_0 = \frac{\hat{n}_{0\sigma}}{\hat{v}_{T\sigma}^3} F_{0\sigma} \\
    *  
    * \f]
    **/ 
   void PerturbationPSFMode(const CComplex f0[NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB],
                                  CComplex f [NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB]);
   /**
   *   @brief document me please
   **/
   void PerturbationPSFExp(const CComplex f0[NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB],
                                 CComplex f [NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB]);
   /**
   *   @brief document me please
   **/
   void PerturbationPSFNoise(const CComplex f0[NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB],
                                   CComplex f [NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB]);
   /**
   *   @brief document me please
   **/
   void PerturbationHermitePolynomial(const CComplex f0[NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB],
                                            CComplex f [NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB]);

  protected:

    virtual void printOn(std::ostream &output) const ;
     
    virtual void initData(FileIO *fileIO) {};
    virtual void writeData(const Timing &timing, const double dt) {};
    virtual void closeData() {};



};


#endif // __INIT_H_
