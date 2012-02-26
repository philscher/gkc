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

#include "Geometry.h"




#define SPECIES_MAX 8

  typedef struct Species {
    Species() : q(0.), m(0.), collision(0.), w_n(0.), w_T(0.), doGyro(true), T0(0.), n0(0.) 
  { 
    alpha   = 0.;
    sigma   = 0.;
    T.resize(RxLB); T = 0.;
    n.resize(RxLB); n = 0.;
  };
    double q;
    double m;
    double collision;
    double w_n;
    double w_T;
    bool   doGyro;
    double T0;
    double n0;
    /** perform gyro average for this species / for adiab. species use flux average */
    double scale_v;
    double scale_n;
    double sigma;
    double alpha;
    //std::string name;
    char name[64];
    char n_name[64];
    char T_name[64];
    // stupid fix, but we have to otherwise all stuff is private
    void update(Geometry<HELIOS_GEOMETRY> *geo, double cs) { 
        scale_v = sqrt(2.*T0/m); 
        scale_n = n0;
        alpha = scale_v*  1./(cs*sqrt(geo->eps_hat));
        sigma = q / T0;
    };
    double debye2(const int x) { return T(x)/(4.*M_PI*n(x)*q*q); };
    Array1d T;
    Array1d n;
} _Species;

class Plasma : public IfaceHelios {



public:
double n_ref, L_ref, T_ref;

  /**  The normalized Debye length (from e.g. G\"orler PhD Eq. (2.82))
     *
     * \f[
     *      \lambda_D = \lambda_D / \rho_{ref} = \sqrt{T_{ref}}{4\pi\rho_{ref}^2 n_{ref} e^2}
     * \f]
     *
     * */
   double debye2;
   bool global;
   double B0, beta, w_p;
   int nfields;
  // const int nfields;
    /**
     *   Note Species goes from 0 ... SPECIES_MAX, where 0 is an
     *   adiabatic species.
     *
     *
     *
     *
     */
    Array<Species, 1>  species;
    bool adiab_species;

    Plasma(Setup *setup, FileIO *fileIO, Geometry<HELIOS_GEOMETRY> *geo, const int nfields=1);
    virtual ~Plasma() {};

    /** Sound speed of ions $c_s = \sqrt{\frac{T_{e0}}{m_i}}$ */
    double cs;
protected:
    virtual void printOn(ostream &output) const;
public:
    void initDataOutput(FileIO *fileIO) ;
     virtual void writeData(Timing *timing) {};
     virtual void closeData() {};
};


#endif // CONFIG_H__

