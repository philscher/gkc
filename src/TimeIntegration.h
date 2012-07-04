/*
 * =====================================================================================
 *
 *       Filename: TimeIntegration.h
 *
 *    Description:Time Integration Interface for gkc. 
 *
 *         Author: Paul P. Hilscher (2011), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */



#include "config.h"
#include "Global.h"

#ifndef TIMEINTEGRATION_H__
#define TIMEINTEGRATION_H__

#include "TestParticle.h"
#include "Setup.h"
#include "Grid.h"
#include "Parallel.h"
#include "Init.h"
#include "Eigenvalue.h"
#include "Fields.h"

class TimeIntegration {
  protected:
       std::string timeIntegrationScheme;
       Parallel *parallel;
    public:
         bool useCFL;
         double maxCFLNumber;
         double linearSafetyFactor;
         double dt;
    time_t start_time;
    Timing maxTiming;
    double maxLinearTimeStep ;
 /** 
   *
   * The maximum TimeStep is defined from istability consideration of the Time Integration
   * method. Following  Karniadakis(2009), Parallel Scientific Computing in C++, it is defined
   * by the Eigenvalue problem
   * \f[
		\frac{dU}{dt} = \lambda U
   * for Runge-Kutta Time Step we get
   * \f[
   *    U^{n+1} = U^n + \frac{1}{2} \delta t \left[ X_1 + X_2 + 2 X_3 + X_4  \right] 
   * \f]
   * with growth factor \f[ G = [ 1 + ... + \frac{\lambda^4 dt^4}{24} ] \f] 
   * thus we get the stability aarea by solving following equation
   * \f[ 
   *    1 + \mu + \frac{\mu^2}{2}{ + \frac{\mu^3}{6} + \frac{\mu^4}{24} = e^{i \theta}
   * \f]
   *
   * we can apply it also also to other  
   * with maximum mu
   *
   *       Integration Scheme  max(|\mu|) 
   *             RK-1   | 2.0
   *             RK-2   | 2.197
   *             RK-3   | 2.537
   *             RK-4   | 2.96
   *            SRK-4   |
   *

 
   */

  /* In eigenvalue
   *  Out Maxumum timestep
   */ 
    double getMaxTimeStepFromEigenvalue(cmplxd max_abs_eigv);
    virtual double solveTimeStep(Vlasov *vlasov, Fields *fields, TestParticles *particles, Timing &timing);
    
    TimeIntegration(Setup *setup, Grid *grid, Parallel *parallel, Vlasov *vlasov, Fields *fields, Eigenvalue *eigenvalue);
    int writeTimeStep(Timing timing, Timing maxTiming, double dt);
   
    
  private:


      /**
       *    Solve Gyro-kinetic equation using explicit Runge-Kutta fourth order (RK4) Integration 
       *
       *    y_{n+1} = y_n + \frac{1}{6} dt \left[ k_1 + 2k_2 + 2 k_3 +k_4 \right]
       * 
       *    \f[
       *        k_1  = f(t_n, y_n)
       *        k_2  = f(t_n+1/2dt, y_n + 1/2 h k_1)
       *        k_3  = f(t_n + 1/2dt, y_n + 1/2 h k_2)
       *        k_4 = f(t_n+dt, y_n+hk_3)
       *    \f]
       *
       *    To save computational time, calculation of k_1 and the update to y_{n+1} is performed
       *    at the last timestep.
       * */
        void solveTimeStepRK4(Fields *fields, Vlasov *vlasov,TestParticles *particles, Timing timing, double  dt);

      /**
       *    Solve Gyro-kinetic equation using explicit Runge-Kutta fourth order (RK4) Integration 
       *
       *    y_{n+1} = y_n + \frac{1}{6} dt \left[ k_1 + 2k_2 + 2 k_3 +k_4 \right]
       * 
       *    \f[
       *        k_1  = f(t_n, y_n)
       *        k_2  = f(t_n+1/2dt, y_n + 1/2 dt k_1)
       *        k_3  = f(t_n + dt, y_n - dt * k_1+ 2 dt k_2)
       *    \f]
       *
       *    To save computational time, calculation of k_1 and the update to y_{n+1} is performed
       *    at the last timestep.
       * */
        void solveTimeStepRK3(Fields *fields, Vlasov *vlasov,TestParticles *particles, Timing timing, double  dt);

      /**
       *    Solve Gyro-kinetic equation using explicit Runge-Kutta second order (RK2) Integration 
       *
       * */
        void solveTimeStepRK2(Fields *fields, Vlasov *vlasov,TestParticles *particles, Timing timing, double  dt);


      /**
       *    Solve Gyro-kinetic equation using explicit Heun's method  (second order) Integration 
       *    (unstable scheme)
       * */
        void solveTimeStepHeun(Fields *fields, Vlasov *vlasov,TestParticles *particles, Timing timing, double  dt);


      /**
       *    Eigenvalues caluclation do not require an timestep integration
       *
       **/
        void solveTimeStepEigen(Fields *fields, Vlasov *vlasov,Timing timing, double  dt);

            
        
            virtual void printOn(ostream &output) const {


            output << "CFL        |  maxCFL Number : " << maxCFLNumber   << std::endl;
            output << "Max Timing |  " << maxTiming << std::endl;
        
         }

   
};


        

#endif // TIMEINTEGRATION_H__
