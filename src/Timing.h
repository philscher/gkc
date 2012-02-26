/*
 * =====================================================================================
 *
 *       Filename: Timing.h
 *
 *    Description: Handles timing events
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __TIMING_H
#define __TIMING_H

#include "Global.h"


#include<ctime>


/** Timing class - governs timing etc.
 */
struct Timing {
    
    int    step;
    double time;
    
    Timing(int _step=0, double _time=0.) : step(_step), time(_time) {};

    bool operator<=(Timing &b) ;
    bool operator<=(int step2) ;
    bool operator!=(int step2) ;
    bool operator!=(double time2) ;
    bool operator%(Timing &b) ;
    bool check(Timing &b, double dt) ;
    double operator&(Timing h_time) ;
    //friend ostream& operator<<(ostream& output, const Timing& t);
    friend ostream& operator<<(ostream& output, const Timing& t) {
      output << "Steps : " << t.step << " Time : " << t.time << std::endl;
      return output;
    };
    /**
     *
     * Show time in hour 
     *
     *
     * */
    static std::string getRemainingTimeString(Timing timing, Timing maxTiming, time_t start_time) ;
    static std::string TimeStringFromSeconds(int secs) ;

};




#endif // __TIMING_H
