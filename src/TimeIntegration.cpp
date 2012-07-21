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


TimeIntegration::TimeIntegration(Setup *setup, Grid *grid, Parallel *_parallel, Vlasov *vlasov, Fields *fields, Eigenvalue *eigenvalue) : parallel(_parallel) {

    start_time = std::time(0); 
    timeIntegrationScheme = setup->get("Helios.TimeIntegration", "Explicit_RK4");


    
    linearSafetyFactor  = setup->get("Helios.LinearSafetyFactor", 0.8);

     linearTimeStep = setup->get("Helios.LinearTimeStep", "Eigenvalue");
    
    useCFL         = setup->get("Helios.useCFL", 1);
    maxCFLNumber   = setup->get("Helios.maxCFLNumber", 0.4);
    
    
    maxTiming.time = setup->get("Helios.MaxTime", -1.);
    maxTiming.step = setup->get("Helios.MaxSteps", -1);
    
    dt             = maxTiming.time / maxTiming.step;
    if (dt <= 0.) dt = maxLinearTimeStep;

};


void TimeIntegration::setMaxLinearTimeStep(Eigenvalue *eigenvalue, Vlasov *vlasov, Fields *fields, const double lin_dt)
{
  if(lin_dt <= 0.) {
    if     (linearTimeStep == "CFL"       ) maxLinearTimeStep = vlasov->getMaxTimeStep(DIR_V, maxCFLNumber);
    else if(linearTimeStep == "Eigenvalue") maxLinearTimeStep = getMaxTimeStepFromEigenvalue(eigenvalue->getMaxAbsEigenvalue(vlasov, fields));
    else if(linearTimeStep == "Fixed"     ) maxLinearTimeStep = dt;
    else                                   maxLinearTimeStep = atof(linearTimeStep.c_str());
  } 
  else maxLinearTimeStep = lin_dt;

}

double TimeIntegration::getMaxTimeStepFromEigenvalue(cmplxd max_abs_eigv)
{

        parallel->print("Using maximum absolute eigenvaue for timestep calulcations");

        // simple assume simple RK-4
	    return (2.96 * linearSafetyFactor /  abs(max_abs_eigv));

};

 double TimeIntegration::solveTimeStep(Vlasov *vlasov, Fields *fields, TestParticles *particles, Timing &timing) 
{
            // OK, we assume only posiitive timesteps ffor now
  
            // set maximum timestep. We reduce intial timsteps in case  initial condition is violant (for right side take care of overflow)
        	if(useCFL == true) dt = min(maxLinearTimeStep, vlasov->getMaxTimeStep(DIR_XY, maxCFLNumber)) ; //* ((timing.step <= 100) ? min(1., 1.e-3 * (timing.step)) : 1));

            if     (timeIntegrationScheme == "Explicit_RK4" ) solveTimeStepRK4 (fields, vlasov, particles, timing, dt);
            else if(timeIntegrationScheme == "Explicit_RK3" ) solveTimeStepRK3 (fields, vlasov, particles, timing, dt);
            else if(timeIntegrationScheme == "Explicit_Heun") solveTimeStepHeun(fields, vlasov, particles, timing, dt);
            else   check(-1, DMESG("No such Integration Scheme"));
 
            timing.time += dt;
            timing.step++;
                
            writeTimeStep(timing, maxTiming, dt);
    
            return dt;
 };


 void TimeIntegration::solveTimeStepRK4(Fields *fields, Vlasov *vlasov,TestParticles *particles, Timing timing, double  dt) {

        
        // Runge-Kutta step 1
	    fields->solve(vlasov->f0,vlasov->f, timing);
        vlasov ->solve(fields, vlasov->f  , vlasov->fs, 0.5e0*dt , 1, BOUNDARY_DIRTY);
        particles->integrate(vlasov, fields, 1);
        // Runge-Kutta step 2
        fields->solve(vlasov->f0,vlasov->fs, timing);
        vlasov ->solve(fields, vlasov->fs , vlasov->fss, 0.5e0*dt, 2, BOUNDARY_DIRTY);
        particles->integrate(vlasov, fields, 2);

    
        // Runge-Kutta step 3
	    fields->solve(vlasov->f0,vlasov->fss, timing);
        vlasov ->solve(fields, vlasov->fss, vlasov->fs,  dt      , 3, BOUNDARY_DIRTY);
        particles->integrate(vlasov, fields, 3);

        // Runge-Kutta step 4
        fields->solve(vlasov->f0,vlasov->fs, timing, 4);
        vlasov ->solve(fields, vlasov->fs , vlasov->f ,  dt/6.e0 , 4, BOUNDARY_DIRTY);
        particles->integrate(vlasov, fields, 4);


      };

   void TimeIntegration::solveTimeStepRK3(Fields *fields, Vlasov *vlasov,TestParticles *particles, Timing timing, double  dt) {

  
        // Runge-Kutta step 1
        fields->solve(vlasov->f0,vlasov->f, timing);
        vlasov ->solve(fields, vlasov->f  , vlasov->fs, 0.5e0*dt , 1, BOUNDARY_DIRTY);


        // Runge-Kutta step 2
        fields->solve(vlasov->f0,vlasov->fs, timing);
        vlasov ->solve(fields, vlasov->fs , vlasov->fss, dt, 2, BOUNDARY_DIRTY);

        // Runge-Kutta step 3
        fields->solve(vlasov->f0,vlasov->fs, timing);
        vlasov ->solve(fields, vlasov->fss, vlasov->f,  dt/6.      , 3, BOUNDARY_DIRTY);

        };
      
   void TimeIntegration::solveTimeStepRK2(Fields *fields, Vlasov *vlasov,TestParticles *particles, Timing timing, double  dt) {

  
        // Runge-Kutta step 1
        fields->solve(vlasov->f0,vlasov->f, timing);
        vlasov ->solve(fields, vlasov->f  , vlasov->fs, 0.5e0*dt , 1, BOUNDARY_DIRTY);


        // Runge-Kutta step 2
        fields->solve(vlasov->f0,vlasov->fs, timing);
        vlasov ->solve(fields, vlasov->fs , vlasov->fss, 0.5e0*dt, 2, BOUNDARY_DIRTY);

       };


      void TimeIntegration::solveTimeStepHeun(Fields *fields, Vlasov *vlasov,TestParticles *particles, Timing timing, double  dt) {

  
        // Runge-Kutta step 1
        fields->solve(vlasov->f0,vlasov->f, timing);
        vlasov ->solve(fields, vlasov->f  , vlasov->fs, 2./3.*dt , 1, BOUNDARY_DIRTY);


        // Runge-Kutta step 2
        fields->solve(vlasov->f0,vlasov->fs, timing);
        vlasov ->solve(fields, vlasov->fs , vlasov->f, 1./4.*dt, 2, BOUNDARY_DIRTY);

       };

      void TimeIntegration::solveTimeStepEigen(Fields *fields, Vlasov *vlasov,Timing timing, double  dt) {

  
        // Runge-Kutta step 1
        fields->solve(vlasov->f0,vlasov->f, timing);
        vlasov ->solve(fields, vlasov->f  , vlasov->fs, 0.5e0*dt , 1, BOUNDARY_DIRTY);

        };
       


int TimeIntegration::writeTimeStep(Timing timing, Timing maxTiming, double dt) {
        
//         TerminalIO::print(
         if(parallel->myRank == 0) {
                   std::cout << "\r"   << "Steps  : " << timing.step  << "/" << maxTiming.step 
                   << "       Time : " << timing.time << "/" << maxTiming.time << "  dt  : " << dt 
                   << std::flush;
                  if(timing.step % 20 == 0)  std::cout << Timing::getRemainingTimeString(timing, maxTiming, start_time);
         }

         return 0;
 };

 
 void TimeIntegration::printOn(ostream &output) const {


            output << "CFL        |  maxCFL Number : " << maxCFLNumber   << std::endl;
            output << "Max Timing |  " << maxTiming << std::endl;
        
         }

   
