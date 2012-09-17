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

#include<random/uniform.h>


#include "Grid.h"
#include "Vlasov.h"
#include "Geometry.h"

#include "Parallel.h"


/**
*   @brief Initialized the background field and intial perturbation
*
*
*
**/ 
class Init   : public IfaceGKC {

   /**
   *   @brief document me please
   **/
   double epsilon_0, sigma;  
   
   Geometry *geo;

   /**
   *   @brief document me please
   **/
   void setFieldFromDataFile(Setup *setup, Array4C field, int n, std::string path);
   /**
   *   @brief document me please
   **/
   void setFieldFromFunction(Setup *setup, Array4C field, int n , std::string func);
/**
 *   @brief inits the Maxwellian 
 **/
   void initMaxwellian(Setup *setup, CComplex f0[NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB],
                                     CComplex f [NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB],
                                    const double V[NvGB], const double M[NmGB]);

  public :

   /**
   *   @brief document me please
   **/
   std::string PerturbationMethod;
   /**
   *   @brief document me please
   **/
   Init(Parallel *parallel, Grid *grid, Setup *setup, Vlasov *vlasov, Fields *fields, Geometry *geo);
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
   int PerturbationHermitePolynomial(const CComplex f0[NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB],
                                           CComplex f [NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB]);

  protected:

    virtual void printOn(ostream &output) const {


         output << " Init      | " << PerturbationMethod << std::endl;

    };
     virtual void initDataOutput(FileIO *fileIO) {};
     virtual void writeData(Timing *timing) {};
     virtual void closeData() {};



};


#endif // __INIT_H_
