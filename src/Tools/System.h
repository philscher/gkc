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
#include <ctime>
#include <sys/types.h>
#include <sys/resource.h>
#include <unistd.h>
#include <execinfo.h> // stack trace

#include <iostream>
#include <fstream>

/**
*
*  @brief Class which provided some system functions which may be
*         not portable.
*
*
*
*
*/
class System 
{
 
 public:

  /**
  *   @brief Increase stack size limit
  *  
  *   Calculation of the non-linear ExB term and especially the
  *   diagonostic module uses large static arrays. 
  *   In order to avoid stack overflow, we increase the stack size
  *   limit. 
  *
  *   @note On most supercomputers using a queue system (e.g. helios),
  *         the available stack size is not limited.
  *
  *   @todo This feature is currently only supported under linux.
  *         Possible to define it in class System ?
  *
  **/
  static void set_min_stacksize(int size_Mbytes)
  {
    
    const rlim_t kStackSize = size_Mbytes * 1024 * 1024; 

    struct rlimit rl;

    int result = getrlimit(RLIMIT_STACK, &rl);
    
    if (result == 0)
    {
      if (rl.rlim_cur < kStackSize)
      {
        rl.rlim_cur = kStackSize;
        result = setrlimit(RLIMIT_STACK, &rl);
        if (result != 0) std::cerr << "setrlimit returned : " << result << std::endl;
      }
    }
    
    result = getrlimit(RLIMIT_STACK, &rl);
    //std::cout << " Stack size is : " << rl.rlim_cur/(1024*1024) << " MBytes" << std::endl;

    return;
  };

   
  /**
  *    @brief returns the process id
  *    
  *    Note : For MPI parallelized version this number is not
  *           equal among different processes.
  *
  *    @return process id of this proccess     
  */
  static int getProcessID() 
  {
      // getpid returns unsigned int (take care of casting)
      // casting as we do not want to define new datatype in MPI
      return (int) getpid();
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
  static unsigned int getTime()
  {
    return time(0);
  };

  /**
  *    @brief Get current time as string
  *
  *    @return Current time as string e.g. "1st January 1996"
  */
  static std::string getTimeString()
  {

      time_t start_time = std::time(0); 
      std::string currentTime(std::ctime(&start_time));
   
      return currentTime;
  };

  static void* align_alloc(size_t N, int align)
  {
      void *p;
           posix_memalign((void **) &p, align, N);
      return p;
  };

   
  // provided by Larry Gritz from
  // http://stackoverflow.com/questions/372484/how-do-i-programmatically-check-memory-use-in-a-fairly-portable-way-c-c
  // The amount of memory currently being used by this process, in bytes.
  /// By default, returns the full virtual arena, but if resident=true,
  /// it will report just the resident set in RAM (if supported on that OS).
  static int memoryUsage()
  {
    // Ugh, getrusage doesn't work well on Linux.  Try grabbing info
    // directly from the /proc pseudo-filesystem.  Reading from
    // /proc/self/statm gives info on your own process, as one line of
    // numbers that are: virtual mem program size, resident set size,
    // shared pages, text/code, data/stack, library, dirty pages.  The
    // mem sizes should all be multiplied by the page size.
    size_t size = 0;
    FILE *file = fopen("/proc/self/statm", "r");
    if (file) {
      unsigned long m_vsize = 0; // virtual memory
      unsigned long m_mres  = 0;  // resident memory
      unsigned long m_share = 0; // shared memory
      unsigned long m_text  = 0; // text data
      unsigned long m_lib   = 0; // library 
      unsigned long m_data  = 0; // data

      fscanf (file, "%lu", &m_vsize);  // Just need the first num: vm size
      fclose (file);
      size = (size_t) m_vsize * getpagesize();
    }
   
    return size;
   
  };

  /**
  *   @brief print the process status
  *
  **/
  static std::string getProcessStatusString()
  {
  
    std::fstream filestr ("/proc/self/status", std::fstream::in);
    std::stringstream ss;
    ss << filestr.rdbuf();
    return ss.str();
  };


  /**
  *   @brief print the stack trace of last 32 calls
  *
  *
  *
  **/
  static void printStackTrace() 
  {
    // get backtrace
    const int numStackTraces = 32;

    void *array[numStackTraces];
    size_t size;
    
    char **strings;
    size_t i;
     
    size = backtrace (array, numStackTraces);
    strings = backtrace_symbols (array, size);
    
    std::cerr << "Showing " << size << " stack frames : " << std::endl;  
    
    // Output them
    for (int i = 0; i < size; i++) std::cout << "    :  " << strings[i] << std::endl;
    free (strings);

   }
};

#endif // __SYSTEM_H

