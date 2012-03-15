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

    Control(Setup *setup, Parallel *_parallel, Analysis *_analysis);

    ~Control();
     
    void printLoopStopReason();

    void runningException(int status, char *error_message=NULL);
    bool checkOK(Timing timing, Timing maxTiming);

        void printOn(ostream &output) const;
};


#endif // __CONTROL_H_
