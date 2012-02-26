/*
 * =====================================================================================
 *
 *       Filename:  Vlasov.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/19/2009 07:17:13 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */


#ifndef __VLASOV_BLITZ_H
#define __VLASOV_BLITZ_H


#include<iostream>
#include<fstream>

#include "Global.h"
#include "Vlasov_LenardBernstein.h"

//class VlasovBlitz : public Vlasov_Lorentz {
class VlasovBlitz : public Vlasov_LenardBernstein {
        
        // Some temporary arrays (not all are necesserally intialized !)
        Array3z k2p_phi, dphi_dy, dphi_dx, dAp_dx, dAp_dy;

        // needed for nonLinear transforms
        Array4z nonLinearTerms; 
        Array3z xy_dphi_dx, xy_dphi_dy, xy_df1_dy, xy_df1_dx, xy_phi, xy_f1; 
//        void Vlasov_EM_Linear   (Array6z fs, Array6z fss, Array6z ft, Fields *fields, double dt, int rk_step);
//        void Vlasov_EM_NonLinear(Array6z fs, Array6z fss, Array6z ft, Fields *fields, double dt, int rk_step);
        void Vlasov_2D(Array6z fs, Array6z fss, Array6z ft, Fields *fields, double dt, int rk_step);
        void Vlasov_2D_Island(Array6z fs, Array6z fss, Array6z ft, Fields *fields, double dt, int rk_step);
        
        Array4z calculatePhiNonLinearity(Array5z phi, Array6z fs, const int m, const int s);
        Array4z calculatePhiNonLinearityAA(Array5z phi, Array6z fs, const int m, const int s);
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
        void setupXiAndG(Array6z fs, Fields *fields);
  public:
        VlasovBlitz(Grid *_grid, Parallel *_parallel, Setup *_setup, FileIO *fileIO, Geometry<HELIOS_GEOMETRY> *_geo, FFTSolver *fft); 
        
        int solve(std::string equation_tyoe, Fields *fields, Array6z fs, Array6z fss, double dt, int rk_step, int user_boundary_type=BOUNDARY_CLEAN);
 
        double L_RF, RF0;
  protected :

        void initDataOutput(hid_t groupID, FileIO *fileIO) {};
        void printOn(ostream &output) const;

};


#endif //  __VLASOV_BLITZ_H

