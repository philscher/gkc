/*
 * =====================================================================================
 *
 *       Filename: Vlasov.h
 *
 *    Description: Vlasov Solver Interface
 *
 *         Author: Paul P. Hilscher (2009-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __VLASOV_H
#define __VLASOV_H

#include "Global.h"
#include "Parallel.h"
#include "Setup.h"
#include "Fields.h"
#include "Grid.h"
#include "FFTSolver.h"

#include "Collisions.h"


enum Boundary  { BOUNDARY_DIRTY=0,  BOUNDARY_CLEAN=1};
class Event;
//! class for solving the Vlasov equation
/*! We solve the 6 dimensional Vlasov equation (X, Y, Z, Vp), m, s 
 *
 */
class Vlasov : public IfaceHelios {

  friend class Event;
  protected:
        double collisionBeta;
  private:
		
        Array6z f_boundary;
        
        // Array for boundary
        Array6z  SendYu, SendXu, SendYl, SendXl, SendVl, SendVu, SendZl, SendZu; 
        Array6z  RecvYu, RecvXu, RecvYl, RecvXl, RecvVl, RecvVu, RecvZl, RecvZu;
      
        bool boundary_isclean;
        int cleanBoundary(Array6z A);
  protected:
        bool useAntiAliasing;
        FFTSolver *fft;
        Parallel *parallel;
        Grid *grid;
        Collisions     *collisions;
        Setup *setup;
        Geometry<HELIOS_GEOMETRY> *geo;
        double hyper_visc;        
        std::string equation_type;
        bool calculate_nonLinear;

        /**
         *   
         *   Solves gyro-kinetic equation in 2-d plane in sheared geometry. For the derivation, the gyro-kinetic equation
         *   in slab geometry 
         *    
         *    \f[ 
         *          \frac{g_1}{\partial t} = 
         *          - \frac{B_0}{B_0^\star} \left( \omega_n + \omega_T \left(v^2+\mu B_t \right) \right) 
         *            \frac{\partial \Xi}{\partial y} F_{0\sigma} - \alpha v_\parallel \frac{\partial G}{\partial z} 
         *    \f]
         *
         *   Sheared geometry is expressed as \f[ B_0 = \left( 0, -x/L_s, 1 \right) \f] with L_s the shearing length. The parallel component
         *   along the magnetic field line is thus calculated according to 
         *   \f[ k_\parallel = B_0 * k = \left( 0, -x/L_s, 1 \right) \cdot (k_x, k_y, k_z) = k_z - x \hat{s} k_y \f]
         *   where in our normalization \\f[ \hat{s} \f]  is defined as \f[ \hat{s} = 1/L_s \f] . 
         *   As $k_z \ll k_y $ we assumtion $k_\parallel = x \hat{s} k_y$ is valid. So in the Vlasov equation the z-derivative
         *   is replaced by \f[ \partial_z = \hat{s} x \partial_y \f].
         *
         *   Note : Thus is not field lign average, thus k_z corresponds to z-direction which is NOT along the magnetic field line
         *
         *  References :
         *
         *          Wang, Z. X.; Li, J. Q.; Kishimoto, Y.; Dong, J. Q.;  Magnetic-island-induced ion temperature gradient mode ; PoP (2009)
         *          Dong, J. Q.; Guzdar, P. N.; Lee, Y. C.            ;  Finite beta effects on ion temperature gradient driven modes ; Phys.of Fluid (1987)
         *
         *
         * */
  public:
	    Array1d  Xi_max;
        Array6z f0, f, fs, fss, ft, f1;
        Array4z G, Xi;

        Vlasov(Grid *grid, Parallel *parallel, Setup *setup, FileIO *fileIO, Geometry<HELIOS_GEOMETRY> *geo, FFTSolver *fft);
        virtual ~Vlasov() {};

        // ! Solve Vlasov equation for current timestep
        int solve(Fields *fields, Array6z fs, Array6z fss, double dt, int rk_step, int user_boundary_type=BOUNDARY_CLEAN);
        virtual int solve(std::string equation_type, Fields *fields, Array6z fs, Array6z fss, double dt, int rk_step, int user_boundary_type=BOUNDARY_CLEAN) = 0;
        int setBoundary(Array6z  A, int boundary_type=BOUNDARY_CLEAN);


        std::string getEquationType() { return equation_type;};

        // for global f simulation we update f
        int updateMaxwellian() { return HELIOS_SUCCESS;};

        /**
         *   Calculates the timestep according to a defined CFL number. For the Vlasov equation
         *   several terms needs to be calculated for the highest possible timestep. These are
         *
         *
         *   \[f \frac{\D \phi} {\phi x} {\frac{\D \phi}{\D y} \f] 
         *   
         *
         *
         *
         *
         * */
        double getMaxTimeStep(int dir, const double maxCFL);



        inline void updateCFL(const cmplxd dphi_dx, const cmplxd dphi_dy, const cmplxd dphi_dz) {

             Xi_max(DIR_X) = max(Xi_max(DIR_X), abs(dphi_dx));
             Xi_max(DIR_Y) = max(Xi_max(DIR_Y), abs(dphi_dy));
             Xi_max(DIR_Z) = max(Xi_max(DIR_Z), abs(dphi_dz));

        };
  protected :

        void printOn(ostream &output) const;


     virtual void initDataOutput(hid_t groupID, FileIO *fileIO) = 0;
     virtual void initDataOutput(FileIO *fileIO) {

        // Data Output
        //
     //##########################     // For the phase space function //###########################
     
          //////////////////////////////////////////////////////////////// Phasespace Group ////////////////////////////////////////////////////////
          hid_t psfGroup = check(H5Gcreate(fileIO->getFileID(), "/Vlasov",H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), DMESG("Error creating group file for Phasespace : H5Gcreate"));
          initDataOutput(psfGroup, fileIO);
     /*
          hsize_t psf_offset[7] =  { NsLlB-1, NmLlB-1, NvLlB-1, NzLlB-1, NyLlB-1, NxLlB-1, 0 };
        
         // FA_psf->create(psfGroup, "Data", psf_offset);
         // FA_psfTime->create(psfGroup, "Timing", offset0);
        
        
        
        hsize_t psf_dim[]       = { grid->NsGD, grid->NmGD, grid->NvGD, grid->NzGD, grid->NyGD, grid->NxGD,  1 };
     hsize_t psf_maxdim[]    = { grid->NsGD, grid->NmGD, grid->NvGD, grid->NzGD, grid->NyGD, grid->NxGD, H5S_UNLIMITED};
     hsize_t psf_moffset[]   = { 0, 0, 2, 2, 2, 2, 0 };
     hsize_t psf_chunkBdim[] = { grid->NsGD, grid->NmGD, grid->NvLB, grid->NzLB,  grid->NyLB, grid->NxLB, 1};
     hsize_t psf_chunkdim[]  = {NsLD, NmLD, NvLD, NzLD, NyLD, NxLD, 1};
     
    // FA_psf      = new FileAttr("Unnamed",7, psf_dim, psf_maxdim, psf_chunkdim, psf_moffset,  psf_chunkBdim, true);
    // FA_psfTime  = new FileAttr("Unnamed",1, time_dim, timing_maxdim, timing_chunkdim, offset0,  timing_chunkdim, parallel->myRank == 0, timing_tid);
     */
          
          H5Gclose(psfGroup);
     };
     virtual void writeData(Timing *timing) {};
     virtual void closeData() {};

};


#endif //  __VLASOV_H

