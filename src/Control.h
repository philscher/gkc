/*
 * =====================================================================================
 *
 *       Filename: Control.h
 *
 *    Description: Handles signals and other triggers
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __CONTROL_H_
#define __CONTROL_H_

#include "Global.h"
#include "Setup.h"
#include "Parallel.h"
#include "Analysis.h"


#include <fenv.h>
#include <csignal>

extern int control_triggered_signal;


class Control {

  class CntrlID {
    bool ok;
    std::string error_message;
    public:
    CntrlID() { ok=true; error_message="";};

    bool check(bool isOK, std::string message) {
        // Append messeges of failure
        if(isOK == false) {
            ok = false;
            error_message += std::string("\n") + message;
        }
        return isOK && ok;
        };
    bool isOK() { return ok; };
    std::string getMessage() { return error_message;};
    };

    Parallel *parallel;
    Analysis *analysis;
    CntrlID  cntrl;
    double maxElectricEnergy, maxMagneticEnergy, maxKineticEnergy;

    void setSignalHandler();

    time_t startTime;
    int maxRunningTime;
     

public:
    std::string cntrl_file_name;

    Control(Setup *setup, Parallel *_parallel, Analysis *_analysis) : parallel(_parallel), analysis(_analysis) {

        // set our control file, this is same for all processes and set to mpi root process id 
        cntrl_file_name   = setup->get("Control.StopFileName", "helios_" + Setup::number2string(parallel->master_process_id) + ".stop");
        maxKineticEnergy  = setup->get("Control.MaxKineticEnergy", 1.e10);
        maxElectricEnergy = setup->get("Control.MaxElectricEnergy", 1.e10);
        maxMagneticEnergy = setup->get("Control.MaxMagneticEnergy", 1.e10);


        // Some Tests to check the function
        maxRunningTime = Setup::getSecondsFromTimeString(setup->get("Control.MaxRunningTime", "0s"));

        // set start Time
        startTime = time (NULL);



     setSignalHandler();
    };

    ~Control() {
        // clean up stop file
        remove(cntrl_file_name.c_str());
    }
     
    void printLoopStopReason() {

        writeMessage(std::string("\nMain Loop finished due to ") + cntrl.getMessage());


    };

    void runningException(int status, char *error_message=NULL) {
        // catch secondary exceptions
        std::cerr << error_message << std::endl;
    //    if     (status == HELIOS_FINISH) abort_run = 1;
    //    else if(status == HELIOS_EXIT  ) delete fileIO;
   //     else    check(-1, DMESG("No such status"));
    #ifdef HELIOS_PARALLEL_MPI
        parallel->barrier();
    #endif
     
    };


    bool checkOK(Timing timing, Timing maxTiming) {
        
      cntrl.check(timing <= maxTiming, "(1) : Time Limit for simulation reached");
      cntrl.check(ifstream(cntrl_file_name.c_str()) == NULL, "(1) : Manual stop bu using file.stop trigger");
      
      
      cntrl.check(((time(NULL)-startTime) < maxRunningTime) || (maxRunningTime == 0), "(1) : Running Time Limit Reached");
     
   
      // check for some physical limits exceeded
      double phiEnergy, ApEnergy, BpEnergy;
      analysis->getFieldEnergy(phiEnergy, ApEnergy, BpEnergy);
//      cntrl.check(phiEnergy <= maxElectricEnergy , "(2) Kinetic Energy is over Maximum");
      cntrl.check( phiEnergy <= maxElectricEnergy, "(2) Electric Energy is over Maximum");
      //cntrl.check( ApEnergy <= maxElectricEnergy, "(2) Electric Energy is over Maximum");
//      cntrl.check( BpEnergy <= maxMagneticEnergy, "(2) Magnetic Energy is over Maximum");
      
      // check Signals
      cntrl.check(!(control_triggered_signal & SIGINT ), "(3) Interrupted by SIGINT");
      cntrl.check(!(control_triggered_signal & SIGTERM), "(3) Interrupted by SIGTERM");
      cntrl.check(!(control_triggered_signal & SIGUSR1), "(3) Interrupted by SIGUSR1");
      cntrl.check(!(control_triggered_signal & SIGUSR2), "(3) Interrupted by SIGUSR2");

      return cntrl.isOK();
    }

        void printOn(ostream &output) const {
            output << "Control    | phi^2 " << (maxElectricEnergy > 0. ? Setup::number2string(maxElectricEnergy) : "off") << std::endl;

        }
};


#endif // __CONTROL_H_
