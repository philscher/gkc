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

#include "TimeIntegration.h"
#include "Tools/TermColor.h"

TimeIntegration::TimeIntegration(Setup *setup, Grid *grid, Parallel *_parallel, Vlasov *_vlasov,
                                 Fields *_fields, TestParticles *_particles, Eigenvalue *eigenvalue, Benchmark *_bench) 
  : parallel(_parallel), bench(_bench), vlasov(_vlasov), fields(_fields), particles(_particles) 
{

  timeIntegrationScheme = setup->get("TimeIntegration.Scheme"            , "Explicit_RK4");
  linearSafetyFactor    = setup->get("TimeIntegration.LinearSafetyFactor", 0.7           );
  linearTimeStep        = setup->get("TimeIntegration.LinearTimeStep"    , "Eigenvalue"  );
  useCFL                = setup->get("TimeIntegration.useCFL"            , 1             );
  maxCFLNumber          = setup->get("TimeIntegration.maxCFLNumber"      , 0.4           );
  outputRatio           = setup->get("TimeIntegration.StepOutputRatio"   , 100           );
  maxTiming.time        = setup->get("TimeIntegration.MaxTime"           , -1.           );
  maxTiming.step        = setup->get("TimeIntegration.MaxSteps"          , -1            );

}

void TimeIntegration::setMaxLinearTimeStep(Eigenvalue *eigenvalue, Vlasov *vlasov, Fields *fields)
{

  if     (linearTimeStep == "Estimate"  ) maxLinearTimeStep = 1.e-99; // not implemented, use estimate of max(kp)
  else if(linearTimeStep == "Eigenvalue") maxLinearTimeStep = getMaxTimeStepFromEigenvalue(vlasov, fields, eigenvalue);
  else                                    maxLinearTimeStep = std::stod(linearTimeStep);

}

double TimeIntegration::getMaxTimeStepFromEigenvalue(Vlasov *vlasov, Fields *fields, Eigenvalue *eigenvalue)
{
  parallel->print("Using maximum absolute eigenvalue to estimate linear time step limit.");
  
  Complex max_abs_eigv = eigenvalue->getMaxAbsEigenvalue(vlasov, fields);
  
  double max_scheme_eigv = 0.;
  
  if      (timeIntegrationScheme == "Explicit_RK4") max_scheme_eigv = 2.96;
  else if (timeIntegrationScheme == "Explicit_RK3") max_scheme_eigv = 2.30;
  else    check(-1, DMESG("Config File Error : TimeIntegration.Scheme"));
       
  return (max_scheme_eigv * linearSafetyFactor / abs(max_abs_eigv));

}

double TimeIntegration::solveTimeStep(Vlasov *vlasov, Fields *fields, TestParticles *particles, Timing &timing) 
{
  
  // set time-step as minimum between (constant) linear time step and (if enabled) from non-linear dt
  double dt = std::min(maxLinearTimeStep, useCFL ? vlasov->getMaxNLTimeStep(maxCFLNumber) : 1.e99);

  if     (timeIntegrationScheme == "Explicit_RK4" ) solveTimeStepRK4 (timing, dt);
  else if(timeIntegrationScheme == "Explicit_RK3" ) solveTimeStepRK3 (timing, dt);
  else if(timeIntegrationScheme == "Explicit_Heun") solveTimeStepHeun(timing, dt);
  else   check(-1, DMESG("No such Integration Scheme"));

  #pragma omp barrier
  #pragma omp single
  {
     timing.time += dt;
     timing.step++;
     writeTimeStep(timing, maxTiming, dt);
  }

  return dt;

}

// use const reference instead
void TimeIntegration::solveTimeStepRK4(Timing timing, const double dt)
{

  // Runge-Kutta step 1
  const double rk_1[] = { 0., 1., 0.};
  
  // if(threadID == 0) fields->solve(vlasov->f0,vlasov->f, timing) 
  //                   vlasov->solve(fields, vlasov->f  , vlasov->fs ,  0.5*dt , 1, rk_1 );

  fields->solve(vlasov->f0,vlasov->f, timing);
  vlasov->solve(fields, vlasov->f  , vlasov->fs ,  0.5*dt , 1, rk_1 );
  particles->integrate(vlasov, fields, 4);
        
  
  // Runge-Kutta step 2
  const double rk_2[] = { 1., 2., 0.};
  fields->solve(vlasov->f0,vlasov->fs, timing);
  vlasov->solve(fields, vlasov->fs , vlasov->fss,  0.5*dt, 2, rk_2);
  particles->integrate(vlasov, fields, 3);

  // Runge-Kutta step 3
  const double rk_3[] = { 1., 2., 0.};
  fields->solve(vlasov->f0,vlasov->fss, timing);
  vlasov->solve(fields, vlasov->fss, vlasov->fs,  dt      , 3, rk_3);
  particles->integrate(vlasov, fields, 2);

  // Runge-Kutta step 4
  const double rk_4[] = { 1., 0., 1.};
  fields->solve(vlasov->f0,vlasov->fs, timing);
  vlasov ->solve(fields, vlasov->fs , vlasov->f ,  dt/6., 4, rk_4);
  particles->integrate(vlasov, fields, 1);

  return;
}

void TimeIntegration::solveTimeStepRK3(Timing timing, const double dt) 
{

  // Runge-Kutta step 1
  const double rk_1[] = { 0., 1., 0.};
  fields->solve(vlasov->f0,vlasov->f, timing);
  vlasov->solve(fields, vlasov->f  , vlasov->fs,  1./3. * dt , 3, rk_1);

  // Runge-Kutta step 2
  const double rk_2[] = { 1./3., 0., 0.};
  fields->solve(vlasov->f0,vlasov->fs, timing);
  vlasov->solve(fields, vlasov->fs , vlasov->fss,  2./3. *dt, 2, rk_2);

  // Runge-Kutta step 3
  const double rk_3[] = { 1., 0., 1.};
  fields->solve(vlasov->f0,vlasov->fss, timing);
  vlasov->solve(fields, vlasov->fss, vlasov->f, 3./4. * dt , 1, rk_3);
        
  return;
}
      

void TimeIntegration::solveTimeStepRK2(Timing timing, const double dt) 
{

 /* 
  
  // Runge-Kutta step 1
  fields->solve(vlasov->f0,vlasov->f, timing);
  vlasov ->solve(fields, vlasov->f  , vlasov->fs, 0.5e0*dt , 1);

  // Runge-Kutta step 2
  fields->solve(vlasov->f0,vlasov->fs, timing);
  vlasov ->solve(fields, vlasov->fs , vlasov->fss, 0.5e0*dt, 2);
 */ 
}

void TimeIntegration::solveTimeStepHeun(Timing timing, const double dt) 
{
 /*   
   // Runge-Kutta step 1
   fields->solve(vlasov->f0,vlasov->f, timing);
   vlasov ->solve(fields, vlasov->f  , vlasov->fs, 2./3.*dt , 1);

   // Runge-Kutta step 2
   fields->solve(vlasov->f0,vlasov->fs, timing);
   vlasov ->solve(fields, vlasov->fs , vlasov->f, 1./4.*dt, 2);
  * */ 
}

void TimeIntegration::solveTimeStepEigen(Timing timing, const double dt) 
{
  
  // Runge-Kutta step 1
  const double rk_0[] = { 0., 0., 0.};
  fields->solve(vlasov->f0,vlasov->f, timing);
  vlasov->solve(fields, vlasov->f  , vlasov->fs,  1. , 0, rk_0);
}
       
void TimeIntegration::writeTimeStep(Timing timing, Timing maxTiming, double dt)
{
        
   static const time_t start_time = std::time(0);

  // should I use flush ? For many CPU maybe not good.
  if(parallel->myRank == 0 && !(timing.step % outputRatio)) {
  
    std::cout   << TermColor::lyellow   << "\r"
                << "Steps  : " << timing.step  << "/" << maxTiming.step 
                << "  Time : " << timing.time  << "/" << maxTiming.time 
                << std::setprecision(3) << " Δt : "   << dt << std::flush; 
 
    std::cout << "  Wall Time : " << Timing::TimeStringFromSeconds(std::time(0) - start_time);
    std::cout << Timing::getRemainingTimeString(timing, maxTiming, start_time);
    std::cout << TermColor::cdefault;

  }
  
  return;
}

void TimeIntegration::printOn(std::ostream &output) const 
{
  output << "Time Int.  |  " << timeIntegrationScheme << "  maxCFL Number : " << maxCFLNumber   << std::endl;
  output << "Max Timing |  " << maxTiming << "  lin. TimeStep : " << linearTimeStep <<  std::endl;
}

