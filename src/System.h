/*
 * =====================================================================================
 *
 *       Filename: System.h
 *
 *    Description: Abstract layer for system calls
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __SYSTEM_H
#define __SYSTEM_H



#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <iostream>

/**
*    @brief Class which provided some system functions which may be
            not portable.
*
*
*
*
*/
class System 
{
   public:
   /**
   *    @brief returns the process id
   *    
   *    Note : For MPI parallelized version this number is not
   *           equal among different processes.
   *
   *    @return process id of this proccess     
   */
   static unsigned int getProcessID() 
   {
      return getpid();
   };


   /**
   *    @brief returns the process id
   *    @param seconds amount of seconds to wait
   *
   */
   static void doSleep(int seconds) 
   {
    
      unsigned int s = sleep((unsigned int) seconds);

   };
    
   /**
   *    @brief Get the number of seconds after 1st January 1970
   *
   *    @return number of seconds after 1st January 1970
   */
   static unsigned int getTime() {

      return time(0);

   };

   /**
   *    @brief Get current time as string
   *
   *    @return Current time as string e.g. "1st January 1996"
   */
   std::string getTimeString()
   {

      time_t start_time = std::time(0); 
      std::string currentTime(std::ctime(&start_time));
	
      return currentTime;
   }

};

#endif // __SYSTEM_H

