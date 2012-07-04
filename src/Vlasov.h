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
class PETScMatrixVector;
//! class for solving the Vlasov equation
/*! We solve the 6 dimensional Vlasov equation (X, Y, Z, Vp), m, s 
 *
 */
class Vlasov : public IfaceHelios {

  /**
   *    Please Document Me !
   *
   **/
  friend class Event;

  /**
   *    Please Document Me !
   *
   **/
  friend class PETScMatrixVector;

  protected:
  /**
   *    Please Document Me !
   *
   **/
        double collisionBeta;
  private:
  /**
   *    Please Document Me !
   *
   **/
        std::string useBoundary;

        Array6z f_boundary;
        
        /**
         *
         *    Temperary buffer for boundary exchange using message passing for
         *    X,Y,Z,V
         *
         *    Note : Add also M, (later S) 
         *
         *    which will be needed for neo-classical collisions
         *
         */ 
        Array6z  SendYu, SendXu, SendYl, SendXl, SendVl, SendVu, SendZl, SendZu; 
        Array6z  RecvYu, RecvXu, RecvYl, RecvXl, RecvVl, RecvVu, RecvZl, RecvZu;
      
  /**
   *    Please Document Me !
   *
   **/
        bool boundary_isclean;
  /**
   *    Please Document Me !
   *
   **/
        int cleanBoundary(Array6z A);
  protected:
  /**
   *    Please Document Me !
   *
   **/
        bool useAntiAliasing;
        FFTSolver *fft;
        Parallel *parallel;
        Grid *grid;
        Collisions     *collisions;
        Setup *setup;
        Geometry<HELIOS_GEOMETRY> *geo;

        /**
         *   Stabilize simulation by adding a small amount of hyper-viscosity.
         *   Especially needed in x-direction to avoid odd-even decoupling.
         *
         *   Is this similar effect as mentioned in 
         *    JJ. Quirk , A contribution to the great Riemann solver debate, 1994 ?
         *
         *   Values between 1.e-5 - 1.e-3 are OK, otherwise it will have impact on
         *   physical results.
         *
         **/
        double hyper_visc[DIR_ALL];        
  /**
   *    Please Document Me !
   *
   **/
        std::string equation_type;
  /**
   *    Please Document Me !
   *
   **/
        bool calculate_nonLinear;
    
        /**
         *
         *  Set the Krook operator 
         *
         *  \frac{\partial g_{1sigma}}{\partial t} = \dots - \nu(x) g_{1\sigma}
         *
         * Is used to damp oscillations close to the simulation boundary.
         *
         *  Note : 
         *          * Is this the no-slip boundary condition ?
         *          * Violates conservation of particles, energy and momentum and
         *            needs to be fixed by modifing the fields. See Lapillone.
         *
         *
         * */
        double krook_nu;

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
        
        /**
         *
         *
         *
         *  Solve Vlasov equation for the current integration step. This is just a wrapper to 
         *  Unfortunatelty right now this is hard-wired. 
         *
         *
         *
         **/
        virtual int solve(std::string equation_type, Fields *fields, Array6z fs, Array6z fss, double dt, int rk_step, int user_boundary_type=BOUNDARY_CLEAN) = 0;
  public:
  /**
   *    Please Document Me !
   *
   **/
	    Array1d  Xi_max;
  /**
   *    Please Document Me !
   *
   **/
        Array6z f0, f, fs, fss, ft, f1;
  /**
   *    Please Document Me !
   *
   **/
        Array4z G, Xi;

  /**
   *    Please Document Me !
   *
   **/
        Vlasov(Grid *grid, Parallel *parallel, Setup *setup, FileIO *fileIO, Geometry<HELIOS_GEOMETRY> *geo, FFTSolver *fft);
        /**
         *  Destructor
         *  Free Arrays
         *
         **/
        virtual ~Vlasov();

        /**
         *
         *  Calls Vlasov::solve(equation type ...)
         *  Handles boundary conditions
         *
         **/
        int solve(Fields *fields, Array6z fs, Array6z fss, double dt, int rk_step, int user_boundary_type=BOUNDARY_CLEAN);

  /**
   *    Please Document Me !
   *
   **/
        int setBoundary(Array6z  A, int boundary_type=BOUNDARY_CLEAN);


  /**
   *    Please Document Me !
   *
   **/
        std::string getEquationType() const { return equation_type;};

        /**
         *    for global f simulation we update f
         *
         *
         *   ToDo : Skeleton, update new Naxwellian background
         *
         */ 
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


        /**
         *
         *
         *   Updated the CFL (Courant-Friedrich-Levy number). For explicit time-stepening
         *   the CFL value has to be always < 1 to ensure stability of the system
         *   (practically < 0.4).
         *   
         *   Note : Stability is still not guranteed. As the system is unstable. Thus the
         *          time-steppening scheme needs to allows imaginary values e.g.
         *          (RK-3, RK-4, Heun method).
         *
         *   Calculated using ....
         *
         * */
        inline void updateCFL(const cmplxd dphi_dx, const cmplxd dphi_dy, const cmplxd dphi_dz) {

             Xi_max(DIR_X) = max(Xi_max(DIR_X), abs(dphi_dx));
             Xi_max(DIR_Y) = max(Xi_max(DIR_Y), abs(dphi_dy));
             Xi_max(DIR_Z) = max(Xi_max(DIR_Z), abs(dphi_dz));

        };
  protected :

  /**
   *    Please Document Me !
   *
   **/
     void printOn(ostream &output) const;

  /**
   *    Please Document Me !
   *
   **/
     virtual void initDataOutput(FileIO *fileIO);
  /**
   *    Please Document Me !
   *
   **/
     virtual void writeData(Timing *timing) {};
  /**
   *    Please Document Me !
   *
   **/
     virtual void closeData() {};

};


#endif //  __VLASOV_H

