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

#include<sstream>
#include<ostream>

#include "Timing.h"


bool Timing::operator<=(Timing &b) 
{
  bool check_step = (step <= b.step) ||  (step == -1 ) || (b.step == -1 );
  bool check_time = (time <= b.time) ||  (time == -1.) || (b.time == -1.);
    
  return (check_step && check_time);
}

bool Timing::operator<=(int step2) 
{
  return (step <= step2 );
}

bool Timing::operator!=(int step2) 
{
  return (step != step2 );
}

bool Timing::operator!=(double time2) 
{
  return (time != time2 );
}

bool Timing::operator%(Timing &b) 
{

  bool doWrite = false;
        
  if((step != -1 ) && (b.step != -1 )) doWrite |=  (step % b.step     == 0 );

  // Note : We can have variable time steps, thus fmod is most of the time != 0.
  if((time > 0.) && (b.time > 0.))     doWrite |=  (fmod(time,b.time) == 0. );
  return doWrite;
}


bool Timing::check(const Timing &b, const double dt) const 
{

  bool doWrite = false;
   
  if((step != -1 ) && (b.step != -1 )) doWrite |=  (step % b.step     == 0 );

  // Note : We can have variable time steps, thus fmod is most of the time != 0.
  if((time > 0.) && (b.time > 0.))     doWrite |=  (fmod(time,b.time) < dt );
  return doWrite;
}


double Timing::operator&(Timing h_time) 
{
  return (time > 0.) ? fmod(h_time.time, time) : 1.e99;
}

std::string Timing::getRemainingTimeString(Timing timing, Timing maxTiming, time_t start_time) 
{

  time_t curr_time = std::time(0); 
  int  time = (int) difftime(curr_time , start_time);

  std::string remainTimeString;
  // extrapolate time 
   
  int remainingSeconds = -1;

  if((timing.step != 0) && (timing.time != 0.)) {

    // BUG : if 1.e6 is changes to 1.e10, ifort crashes due to SIGFPE, gcc not, why ??
    int sec_step = (int) (maxTiming.step > 0 ? (int) ((1.-(double) timing.step/(double) maxTiming.step)/((double) timing.step/(double) maxTiming.step) * (double) time) : 1e6);// - time;
    int sec_time = (int) (maxTiming.time > 0 ? (int) ((1.-timing.time/maxTiming.time)/(timing.time/maxTiming.time) * (double) time) : 1.e6);// - time;

    remainingSeconds = std::min(sec_step, sec_time);
   
  }

  if(remainingSeconds > 0) remainTimeString = " ETA : ~ " + TimeStringFromSeconds(remainingSeconds);
  else                     remainTimeString = " ETA :   -d:-h:-m:-s";  

  return remainTimeString;
}

std::string Timing::TimeStringFromSeconds(int secs) 
{

  int days    = secs / (60*60*24 ); secs -= days    * (60*60*24);
  int hours   = secs / (60*60    ); secs -= hours   * (60*60)   ;
  int minutes = secs / (60       ); secs -= minutes *  60       ;
  int seconds = secs;
      
  std::stringstream ss;
            
  if(days > 0   )  ss << days    << "d ";
  if(hours > 0  )  ss << hours   << "h ";
  if(minutes > 0)  ss << minutes << "m ";
  if(seconds > 0)  ss << seconds << "s ";

  return ss.str();
}

