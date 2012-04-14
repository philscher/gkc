/*
 * =====================================================================================
 *
 *       Filename: Vlasov_Cilk.h
 *
 *    Description: Implementation of GK Vlasov's equation using 
 *                 Intel Cilk (Array Notation)
 *
 *         Author: Paul P. Hilscher (2011-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#include "Vlasov.h"

#ifndef __VLASOV_CILK_H
#define __VLASOV_CILK_H


#include<iostream>
#include<fstream>

#include "Global.h"
#include "Vlasov_LenardBernstein.h"

class Event;

//class VlasovCilk : public Vlasov_Lorentz {
class VlasovCilk : public Vlasov_LenardBernstein {
  friend class Event;
        
        void setBoundaryXY(Array3d A, int dir=DIR_XY) ;
        
	    Array3d transform2Real(Array5z A, const int m, const int s);
        Array3d transform2Real(Array6z A, const int v, const int m, const int s);
        // Some temporary arrays (not all are necesserally intialized !)
        Array3z dphi_dy, dphi_dx, dAp_dx, dAp_dy;
        Array4z k2p_phi;
        Array4z calculatePhiNonLinearity(Array5z phi, Array6z fs, const int m, const int s);

        // needed for nonLinear transforms
        Array4z nonLinearTerms; 
        Array3d xy_dphi_dx, xy_dphi_dy, xy_df1_dy, xy_df1_dx, xy_phi, xy_f1;

        Array3d  SendXuR, SendXlR, RecvXuR, RecvXlR;

        
        
        void    Vlasov_2D(
                           cmplxd fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd vf0[NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd phi[NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           //cmplxd k2_phi[NzLD][NkyLD][NxLD],
                           cmplxd k2_phi[NzLD][NkyLD][NxLD],
                           cmplxd nonLinear[NzLD][NkyLD][NxLD][NvLD],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           Fields *fields,
                           const double dt, const int rk_step, Array6z _fs);
       void    Vlasov_2D_Global(
                           cmplxd fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd vf0[NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd phi[NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           //cmplxd k2_phi[NzLD][NkyLD][NxLD],
                           cmplxd k2_phi[NzLD][NkyLD][NxLD],
                           Fields *fields,
                           const double dt, const int rk_step, Array6z _fs);


void  Vlasov_2D_Island(
                           cmplxd fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd vf0[NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd phi[NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           //cmplxd k2_phi[NzLD][NkyLD][NxLD],
                           cmplxd k2_phi[NzLD][NkyLD][NxLD],
                           cmplxd nonLinear[NzLD][NkyLD][NxLD][NvLD],
                           cmplxd dphi_dx[NzLB][NkyLD][NxLB],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           Fields *fields,
                           const double dt, const int rk_step, Array6z _fs);

void  Landau_Damping(
                           cmplxd fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd vf0[NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd phi[NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           //cmplxd k2_phi[NzLD][NkyLD][NxLD],
                           cmplxd k2_phi[NzLD][NkyLD][NxLD],
                           cmplxd dphi_dx[NzLB][NkyLD][NxLB],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           Fields *fields,
                           const double dt, const int rk_step);

void setupXiAndG(
                           const cmplxd g       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd f0       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd phi       [NsLD][NmLD][NzLB][NkyLD][NxLB+4  ],
                           const cmplxd Ap       [NsLD][NmLD][NzLB][NkyLD][NxLB+4  ],
                           const cmplxd Bp       [NsLD][NmLD][NzLB][NkyLD][NxLB+4  ],
                           cmplxd Xi       [NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd G        [NzLB][NkyLD][NxLB  ][NvLB],
                           const double V[NvGB], const double M[NmGB],
                           const int m, const int s);


void Vlasov_EM(
                           cmplxd fs       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd fss      [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd vf0[NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd f1 [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           cmplxd ft       [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                           const cmplxd phi[NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           const cmplxd Ap [NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           const cmplxd Bp [NsLD][NmLD][NzLB][NkyLD][NxLB+4],
                           //cmplxd k2_phi[NzLD][NkyLD][NxLD],
                           cmplxd k2_phi[NzLD][NkyLD][NxLD],
                           cmplxd dphi_dx[NzLB][NkyLD][NxLB],
                           cmplxd Xi       [NzLD][NkyLD][NxLD  ][NvLD],
                           cmplxd G        [NzLD][NkyLD][NxLD  ][NvLD],
                           const double X[NxGB], const double V[NvGB], const double M[NmGB],
                           Fields *fields,
                           const double dt, const int rk_step);



        //void Vlasov_2D(A6z fs, A6z fss, A6z ft);

       /* 
        void Vlasov_2D(cmplxd fs [NsLlB,NsLB][NmLlB,NmLB][NvLlB,NvLB][NzLlB,NzLB][NkyLlD,NkyLB][NxLlB,NxLB],
                       cmplxd fss[NsLlB,NsLB][NmLlB,NmLB][NvLlB,NvLB][NzLlB,NzLB][NkyLlD,NkyLB][NxLlB,NxLB],
                       cmplxd vf0[NsLlB,NsLB][NmLlB,NmLB][NvLlB,NvLB][NzLlB,NzLB][NkyLlD,NkyLB][NxLlB,NxLB],
                       cmplxd f [NsLlB,NsLB][NmLlB,NmLB][NvLlB,NvLB][NzLlB,NzLB][NkyLlD,NkyLB][NxLlB,NxLB],
                       cmplxd ft [NsLlB,NsLB][NmLlB,NmLB][NvLlB,NvLB][NzLlB,NzLB][NkyLlD,NkyLB][NxLlB,NxLB], Fields *fields,
                       const double dt, const int rk_step);

        * */ 
        
//        void Vlasov_2D_Linear_FixedAp(Array6z fs, Array6z fss, Array6z ft, Fields *fields, double dt, int rk_step);
//        void Vlasov_2D_NonLinear_FixedAp(Array6z fs, Array6z fss, Array6z ft, Fields *fields, double dt, int rk_step);
        /**
         *   Note from Goerles PhD Thesis:
         *
         *   Calculated the gyro-average modified potential (Eq. 2.32)
         *
         *   \f[ \bar{\Xi} = \bar{\phi}_1 - \frac{v_\parallel}{c}\bar{A}_{1\parallel} 
         *                          + \frac{\mu}{q_\sigma} \bar{B}_{1\parallel}
         *   \f]
         *
         *  as well as the Gamma abbrevation functions (Eq. 2.50 below)
         *
         *  \f[ \Gamma_{\sigma,\nu} = \partial_\nu F_{1\sigma} + 
         *          \frac{F_{0\sigma}}{T_{0\sigma}} \partial_\nu \left( q_\sigma \bar{\phi}_1 + \mu \bar{B}_{1\parallel}
         *  \f]
         *
         *  and F1 which needs to be reconstructed from g and is needed for the magnetic mirror term (Eq. 2.50)
         *
         *  \f[ F_1 = g_{1\sigma} - \frac{q_sigma}{c}\bar{A}_{1\parallel}\frac{v_\parallel}{T_{0\sigma}} F_{0\sigma} \f]
         *
         *
         */
  public:
        VlasovCilk(Grid *_grid, Parallel *_parallel, Setup *_setup, FileIO * fileIO, Geometry<HELIOS_GEOMETRY> *_geo, FFTSolver *fft); 
        
        int solve(std::string equation_tyoe, Fields *fields, Array6z fs, Array6z fss, double dt, int rk_step, int user_boundary_type=BOUNDARY_CLEAN);
 
        double L_RF, RF0;
  protected :

        void initDataOutput(hid_t vlasovGroup, FileIO *fileIO) {
             
                check(H5LTset_attribute_string(vlasovGroup, ".", "CollisionModel","LennardBernstein"), DMESG("H5LTset_attribute"));
               check(H5LTset_attribute_double(vlasovGroup, ".", "CollisionBeta"   ,  &collisionBeta, 1), DMESG("H5LTset_attribute"));
        
        };
        void printOn(ostream &output) const;

};


#endif //  __VLASOV_BLITZ_H

