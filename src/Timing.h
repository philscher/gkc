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
    
  /**
   *    Please Document Me !
   *
   **/
    int    step;
  /**
   *    Please Document Me !
   *
   **/
    double time;
    
    //Timing(int _step=-1, double _time=-1.) = default: step(_step), time(_time) {};
  /**
   *    Please Document Me !
   *
   **/
    Timing(int _step=-1, double _time=-1.) : step(_step), time(_time) {};

  /**
   *    Please Document Me !
   *
   **/
    bool operator<=(Timing &b) ;
  /**
   *    Please Document Me !
   *
   **/
    bool operator<=(int step2) ;
  /**
   *    Please Document Me !
   *
   **/
    bool operator!=(int step2) ;
  /**
   *    Please Document Me !
   *
   **/
    bool operator!=(double time2) ;
  /**
   *    Please Document Me !
   *
   **/
    bool operator%(Timing &b) ;
  /**
   *    Please Document Me !
   *
   **/
    bool check(Timing &b, double dt) ;
  /**
   *    Please Document Me !
   *
   **/
    double operator&(Timing h_time) ;
  /**
   *    Please Document Me !
   *
   **/
    friend ostream& operator<<(ostream& output, const Timing& t) {
      output << "Steps : " << t.step << " Time : " << t.time << std::endl;
      return output;
    };
    /**
     *
     * Show time in hour 
   *    Please Document Me !
     *
     *
     * */
    static std::string getRemainingTimeString(Timing timing, Timing maxTiming, time_t start_time) ;
  /**
   *    Please Document Me !
   *
   **/
    static std::string TimeStringFromSeconds(int secs) ;

};




#endif // __TIMING_H
