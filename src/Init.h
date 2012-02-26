/*
 * =====================================================================================
 *
 *       Filename:  Init.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/19/2009 06:20:20 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */
#ifndef __INIT_H_
#define __INIT_H_

#include "Global.h"

#include<blitz/array.h>
#include<random/uniform.h>


#include "Grid.h"
#include "Vlasov.h"
#include "Geometry.h"
#include "GeometrySlab.h"
#include "GeometryShear.h"
#include "Geometry2D.h"

#include "Special/HermitePoly.h"


class Init   : public IfaceHelios {
double epsilon_0, sigma;  
  Geometry<HELIOS_GEOMETRY> *geo;
private:
  static double Phase(const int q, const int N) { return 2.*M_PI*((double) (q-1)/N); };
  void    setFieldFromDataFile(Setup *setup, Array4d field, int n, std::string path);
    void setFieldFromFunction(Setup *setup, Array4d field, int n , std::string func);
  public :
  std::string PerturbationMethod;
  Init(Grid *grid, Setup *setup, Vlasov *vlasov, Fields *fields, Geometry<HELIOS_GEOMETRY> *geo);
   ~Init() {};



//*  Function to set Initital conditions
/*!
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
 */ 

int PerturbationPSFMode(Vlasov *vlasov, int s = 1, double pre=0.);
int PerturbationPSFExp(Vlasov *vlasov, int s, double pre);
int PerturbationPSFNoise(Vlasov *vlasov, int s, double pre);
int PerturbationHermitePolynomial(Vlasov *vlasov, int s, double pert, int l);




/* 
double inline Perturbation(int x,int y,int z, const double epsilon_0, const double sigma);
static int setPerturbation3D(Array3z A); 
static int setPerturbation3D(Array3d A);
static int setPerturbation1D(Array1z A);
 * */

    virtual void printOn(ostream &output) const {


         output << " Init      | " << PerturbationMethod << std::endl;

    };
     virtual void initDataOutput(FileIO *fileIO) {};
     virtual void writeData(Timing *timing) {};
     virtual void closeData() {};

    double inline Perturbation(int x,int y,int z, const double epsilon_0, const double sigma);


};


#endif // __INIT_H_
