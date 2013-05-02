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

#ifndef TIMEINTEGRATION_H__
#define TIMEINTEGRATION_H__

#include "Global.h"

#include "Setup.h"
#include "Grid.h"
#include "Parallel/Parallel.h"
#include "Init.h"
#include "Eigenvalue/Eigenvalue.h"
#include "Benchmark/Benchmark_PAPI.h"
#include "Analysis/TestParticle.h"
#include "Fields.h"

/**
*    @brief Class which handles explicit time integration
*
*
*    The maximum time step is defined from stability consideration of the Time Integration
*    method. Following  Karniadakis(2009), Parallel Scientific Computing in C++, it is defined
*    by the Eigenvalue problem
*    \f[ 
*        \frac{dU}{dt} = \lambda U 
*    \f]
*    for Runge-Kutta Time Step we get
*    \f[
*       U^{n+1} = U^n + \frac{1}{2} \delta t \left[ X_1 + X_2 + 2 X_3 + X_4  \right] 
*    \f]
*    with growth factor \f[ G = [ 1 + ... + \frac{\lambda^4 dt^4}{24} ] \f] 
*    thus we get the stability aarea by solving following equation
*    \f[ 
*       1 + \mu + \frac{\mu^2}{2} + \frac{\mu^3}{6} + \frac{\mu^4}{24} = e^{i \theta}
*    \f]
*
*    we can apply it also also to other  
*    with maximum mu
*
*       Integration Scheme  max(|\f$ \mu \f$|) 
*             RK-1   | 2.0
*             RK-2   | 2.197
*             RK-3   | 2.537
*             RK-4   | 2.96
*            SRK-4   |
*
*
**/
class TimeIntegration  : public IfaceGKC 
{
 protected:
   
  Parallel *parallel;
  Fields   *fields;
  Vlasov   *vlasov;
  TestParticles *particles;
  Benchmark *bench;
   
  int outputRatio;            ///< Frequency (time step) of step output
  bool useCFL;                ///< Set if CFL number is used for time-step calculations
  double maxCFLNumber,        ///< Restrict time step to fulfill CFL number
         linearSafetyFactor,  ///< Safety factor for time-step got from eigenvalue calculations
         maxLinearTimeStep ;  ///< Maximum linear time-step

  std::string linearTimeStep;        ///< Type of linear time step
  std::string timeIntegrationScheme; ///< Type of timeIntegrationScheme
   
 public:

  Timing maxTiming;

  /**
  *   Constructor
  *
  **/ 
  TimeIntegration(Setup *setup, Grid *grid, Parallel *parallel, Vlasov *vlasov, Fields *fields, TestParticles *particles, 
                  Eigenvalue *eigenvalue, Benchmark *bench);

  /** 
  *  @brief calculated the maximum allowed linear time step from max. absolute eigenvalue
  *
  *  Calculated maximum allowed time step from the maximum absolute eigenvalue
  *  which is guaranteed to be stable for the current integration scheme.
  *
  *  Note : This is only an approximation, because the positive real eigenvalues
  *         (growing structure) is not necessarily captured.
  *
  *   @param   max_abs_eigv maximum absolute eigenvalue
  *   @return               maximum linear time step
  **/ 
  double getMaxTimeStepFromEigenvalue(Vlasov *vlasov, Fields *fields, Eigenvalue *eigenvalue);

  /**
  * 
  * @brief solves the current time step using predefined integration scheme
  *
  **/
  virtual double solveTimeStep(Vlasov *vlasov, Fields *fields, TestParticles *particles, Timing &timing);
   
  /**
  *
  *
  **/
  void writeTimeStep(Timing timing, Timing maxTiming, double dt);

  /**
  *    @brief ?
  *
  **/
  void setMaxLinearTimeStep(Eigenvalue *eigenvalue, Vlasov *vlasov, Fields *fields);
    
 private:

  /**
  *    Solve Gyro-kinetic equation using explicit Runge-Kutta fourth order (RK4) Integration 
  *
  *    \f[
  *        y_{n+1} = y_n + \frac{1}{6} dt \left[ k_1 + 2k_2 + 2 k_3 +k_4 \right]
  *    \f]
  * 
  *    \f[
  *        k_1  = f(t_n, y_n)
  *        k_2  = f(t_n+1/2dt, y_n + 1/2 h k_1)
  *        k_3  = f(t_n + 1/2dt, y_n + 1/2 h k_2)
  *        k_4 = f(t_n+dt, y_n+hk_3)
  *    \f]
  *
  *    To save computational time, calculation of k_1 and the update to y_{n+1} is performed
  *    at the last time step.
  *
  *    RK-4 integration is probably the best choice, between trade-off of computational cost 
  *    and maximum allowed time step for linear runs. e.g.  RK-3 scheme has a 33% 
  *    lower computational cost, maximum linear time-step is also reduced by 30%. 
  *    However, non-linear time step is restricted by the ExB velocity to fulfill the
  *    CFL condition. Thus in case time integration error can be neglected, RK-3 is the
  *    better choice.
  *
  *    @image html TimeIntegration_RK4_Stability.png
  *
  **/
  void solveTimeStepRK4(Timing timing, const double dt);

  /**
  *    Solve Gyro-kinetic equation using explicit Runge-Kutta fourth order (RK3) Integration 
  *    \f[
  *        y_{n+1} = y_n + \frac{1}{6} dt \left[ k_1 + 2k_2 + 2 k_3 +k_4 \right]
  *    \f]
  * 
  *    \f[
  *        k_1  = f(t_n, y_n)
  *        k_2  = f(t_n+1/2dt, y_n + 1/2 dt k_1)
  *        k_3  = f(t_n + dt, y_n - dt * k_1+ 2 dt k_2)
  *    \f]
  *
  *    To save computational time, calculation of k_1 and the update to y_{n+1} is performed
  *    at the last time step.
  *
  *    @image html TimeIntegration_RK3_Stability.png
  **/
  void solveTimeStepRK3(Timing timing, const double dt);

  /**
  *    Solve Gyro-kinetic equation using explicit Runge-Kutta second order (RK2) Integration 
  *    Note : This is an unstable scheme for growing setups
  *
  *    @image html TimeIntegration_RK2_Stability.png
  *
  *
  **/
  void solveTimeStepRK2(Timing timing, const double dt);

  /**
  *    Solve Gyro-kinetic equation using explicit Heun's method  (second order) Integration 
  *    (unstable scheme)
  **/
  void solveTimeStepHeun(Timing timing, const double dt);

  /**
  *    Eigenvalues calculation do not require an time step integration
  *
  **/
  void solveTimeStepEigen(Timing timing, const double dt);
        
  /**
  *
  **/
  virtual void printOn(std::ostream &output) const ;
};

#endif // TIMEINTEGRATION_H__
