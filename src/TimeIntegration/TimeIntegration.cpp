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


TimeIntegration::TimeIntegration(Setup *setup, Grid *grid, Parallel *_parallel, Vlasov *vlasov, Fields *fields, Eigenvalue *eigenvalue, Benchmark *_bench) 
  : parallel(_parallel), bench(_bench) 
{

    start_time = std::time(0); 
    timeIntegrationScheme = setup->get("Helios.TimeIntegration", "Explicit_RK4");

    linearSafetyFactor  = setup->get("Helios.LinearSafetyFactor", 0.7);
    linearTimeStep      = setup->get("Helios.LinearTimeStep"    , "Eigenvalue");
    useCFL              = setup->get("Helios.useCFL"            , 1);
    maxCFLNumber        = setup->get("Helios.maxCFLNumber"      , 0.4);
    outputRatio         = setup->get("Helios.OutputRatio"       , 1);
    
    
    maxTiming.time = setup->get("Helios.MaxTime" , -1.);
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
    else                                    maxLinearTimeStep = Setup::str2num(linearTimeStep);
  } 
  else maxLinearTimeStep = lin_dt;

}

double TimeIntegration::getMaxTimeStepFromEigenvalue(Complex max_abs_eigv)
{

        parallel->print("Using maximum absolute eigenvaue for timestep calulcations");

        // simple assume simple RK-4
        double max_scheme_eigv = 0.;
        if      (timeIntegrationScheme == "Explicit_RK4") max_scheme_eigv = 2.96;
        else if (timeIntegrationScheme == "Explicit_RK4") max_scheme_eigv = 2.30;
        else    check(-1, DMESG("Max Stable Eigenvalue not defined"));
       
        return (max_scheme_eigv * linearSafetyFactor /  abs(max_abs_eigv));

};

 double TimeIntegration::solveTimeStep(Vlasov *vlasov, Fields *fields, TestParticles *particles, Timing &timing) 
{
            // OK, we assume only posiitive timesteps for now
  
            // set maximum timestep. We reduce intial timsteps in case initial condition is violant (for right side take care of overflow)
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


 void TimeIntegration::solveTimeStepRK4(Fields *fields, Vlasov *vlasov,TestParticles *particles, Timing timing, double  dt)
{

        
        // Runge-Kutta step 1
        const double rk_1[] = { 0., 1., 0.};
        fields->solve(vlasov->f0,vlasov->f, timing);
        vlasov ->solve(fields, vlasov->f  , vlasov->fs, 0.5e0*dt , 1, rk_1 );
        particles->integrate(vlasov, fields, 1);
        
        // Runge-Kutta step 2
        const double rk_2[] = { 1., 2., 0.};
        fields->solve(vlasov->f0,vlasov->fs, timing);
        vlasov ->solve(fields, vlasov->fs , vlasov->fss, 0.5e0*dt, 2, rk_2);
        particles->integrate(vlasov, fields, 2);

    
        // Runge-Kutta step 3
        const double rk_3[] = { 1., 2., 0.};
        fields->solve(vlasov->f0,vlasov->fss, timing);
        vlasov ->solve(fields, vlasov->fss, vlasov->fs,  dt      , 3, rk_3);
        particles->integrate(vlasov, fields, 3);

        // Runge-Kutta step 4
        const double rk_4[] = { 1., 0., 1.};
        fields->solve(vlasov->f0,vlasov->fs, timing);
        vlasov ->solve(fields, vlasov->fs , vlasov->f ,  dt/6.e0 , 4, rk_4);
        particles->integrate(vlasov, fields, 4);



};


void TimeIntegration::solveTimeStepRK3(Fields *fields, Vlasov *vlasov,TestParticles *particles, Timing timing, double  dt) 
{

  
        // Runge-Kutta step 1
        const double rk_1[] = { 0., 1., 0.};
        fields->solve(vlasov->f0,vlasov->f, timing);
        vlasov ->solve(fields, vlasov->f  , vlasov->fs, 1./3. * dt , 1, rk_1);


        // Runge-Kutta step 2
        const double rk_2[] = { 1./3., 0., 0.};
        fields->solve(vlasov->f0,vlasov->fs, timing);
        vlasov ->solve(fields, vlasov->fs , vlasov->fss, 2./3. *dt, 2, rk_2);

        // Runge-Kutta step 3
        const double rk_3[] = { 1., 0., 1.};
        fields->solve(vlasov->f0,vlasov->fss, timing);
        vlasov ->solve(fields, vlasov->fss, vlasov->f,  3./4. * dt , 3, rk_3);
        

};
      

void TimeIntegration::solveTimeStepRK2(Fields *fields, Vlasov *vlasov,TestParticles *particles, Timing timing, double  dt) 
{

 /* 
        // Runge-Kutta step 1
        fields->solve(vlasov->f0,vlasov->f, timing);
        vlasov ->solve(fields, vlasov->f  , vlasov->fs, 0.5e0*dt , 1);


        // Runge-Kutta step 2
        fields->solve(vlasov->f0,vlasov->fs, timing);
        vlasov ->solve(fields, vlasov->fs , vlasov->fss, 0.5e0*dt, 2);
  * */ 


};



void TimeIntegration::solveTimeStepHeun(Fields *fields, Vlasov *vlasov,TestParticles *particles, Timing timing, double  dt) 
{

 /*   
        // Runge-Kutta step 1
        fields->solve(vlasov->f0,vlasov->f, timing);
        vlasov ->solve(fields, vlasov->f  , vlasov->fs, 2./3.*dt , 1);


        // Runge-Kutta step 2
        fields->solve(vlasov->f0,vlasov->fs, timing);
        vlasov ->solve(fields, vlasov->fs , vlasov->f, 1./4.*dt, 2);
  * */ 

};


void TimeIntegration::solveTimeStepEigen(Fields *fields, Vlasov *vlasov,Timing timing, double  dt) 
{
  
        // Runge-Kutta step 1
        const double rk_0[] = { 0., 0., 0.};
        fields->solve(vlasov->f0,vlasov->f, timing);
        vlasov ->solve(fields, vlasov->f  , vlasov->fs, 1. , 0, rk_0);

};
       


void TimeIntegration::writeTimeStep(Timing timing, Timing maxTiming, double dt)
{
        
  // should I use flush ? For many CPU maybe not good.
  
  if(parallel->myRank == 0 && !(timing.step % outputRatio)) {
  
    std::cout   << "\r" << "Steps  : " << timing.step  << "/" << maxTiming.step 
                << "  Time : " << timing.time  << "/" << maxTiming.time 
                << std::setprecision(3) <<   "  dt : " << dt << std::flush; 
  
    std::cout << Timing::getRemainingTimeString(timing, maxTiming, start_time);
  
  }

  
  return;

};

 

void TimeIntegration::printOn(std::ostream &output) const 
{

            output << "Time Int.  |  " << timeIntegrationScheme << "  maxCFL Number : " << maxCFLNumber   << std::endl;
            output << "Max Timing |  " << maxTiming << "  lin. TimeStep : " << linearTimeStep <<  std::endl;
        
}

  
/*  
         /////////////////////  setup Table for CFL  //////////////////////
         cfl_table = new CFLTable();
         
         cfl_offset[0] =  HOFFSET( CFLTable, timeStep );
         cfl_offset[1] =  HOFFSET( CFLTable, time );
         cfl_offset[2] =  HOFFSET( CFLTable, Fx );
         cfl_offset[3] =  HOFFSET( CFLTable, Fy );
         cfl_offset[4] =  HOFFSET( CFLTable, Fz  );
         cfl_offset[5] =  HOFFSET( CFLTable, Fv );
         cfl_offset[6] =  HOFFSET( CFLTable, total );
          

         for(int i = 1; i < 7; i++)  cfl_sizes[i] = sizeof(double); cfl_sizes[0] = sizeof(int);
         hid_t   cfl_type[7]; for(int i = 1; i < 7; i++)  cfl_type [i] = H5T_NATIVE_DOUBLE; cfl_type[0] = H5T_NATIVE_INT;

         const char *cfl_names[7];
         cfl_names[0] = "timeStep";
         cfl_names[1] = "time";
         cfl_names[2] = "Fx"; cfl_names[3] = "Fy"; cfl_names[4] = "Fz"; cfl_names[5] = "Fv"; cfl_names[6] = "Total";

          check(H5TBmake_table("cflTable", file, "cfl", (hsize_t) 7, (hsize_t) 0, sizeof(CFLTable), (const char**) cfl_names,
                               cfl_offset, cfl_type, 32, NULL, 0, cfl_table ), DMESG("H5Tmake_table : cfl"));
  int FileIO::writeCFLValues(Analysis *analysis, Fields *fields, Timing timing) {
        
    cfl_table->timeStep  = timing.step;
    cfl_table->time  = timing.time;
    cfl_table->Fx    = analysis->getMaxTimeStep(timing, DIR_X);
    cfl_table->Fy    = analysis->getMaxTimeStep(timing, DIR_Y);
    cfl_table->Fz    = analysis->getMaxTimeStep(timing, DIR_Z);
    cfl_table->Fv    = analysis->getMaxTimeStep(timing, DIR_V);
    cfl_table->total = analysis->getMaxTimeStep(timing);
    
    check(H5TBappend_records (file, "cfl", 1, sizeof(CFLTable), cfl_offset, cfl_sizes, cfl_table), DMESG("Append Table")); 

  }
   //
  //  @brief Table to store the CFL values and various contributions
   //
   // @note  Is it worth the effort to keep in the codebase ?
   //        Get rid of contributional part, just store the time-step
   // 
   typedef struct CFLTable
   {
     //  CFLTable(int _timeStep, double _time, double _Fx, double _Fy, double _Fz, double _total) : 
     //           timeStep(_timeStep), time(_time), Fx(_Fx), Fy(_Fy), Fz(_Fz), total(_total) {}; 
     int timeStep;
     double time;
     double   Fx, Fy, Fz, Fv, total;
   } _CFLTable;
   
   size_t cfl_offset[7];
   size_t cfl_sizes[7];

   CFLTable *cfl_table;
    
    */     
