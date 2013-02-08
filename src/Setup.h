/*
 * =====================================================================================
 *
 *       Filename: Setup.h
 *
 *    Description: Read in Config parameters
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __GKC_SETUP_H__
#define __GKC_SETUP_H__


// define to compile Fparser with complex number support
//#define FP_SUPPORT_COMPLEX_DOUBLE_TYPE

#include "Global.h"
#include "FunctionParser/fparser.hh"

#include <string>
#include <cctype>
#include <vector>
#include <map>
#include <algorithm>
#include <fstream>



/**
*   @brief Program configuration 
*
*
*   This class provides an interface to the simulation's 
*   setup, which is either provided by a file, and/or
*   command line options.
* 
*   Class also provided interface to argv arguments.
*
*   Reads option like Grid.Nx = 64
*
*
**/ 
class Setup : public IfaceGKC {

  /**
  *   @brief hash which holds all file configuration options
  *
  **/ 
  std::map   <std::string, std::string> config;
  
  /**
  *  @brief array to check access to configuration
  *
  **/ 
  std::vector<std::string> config_check;

  std::string commandLineOptions,  ///< line options specified by -x
                extraLineOptions;  ///< line options specified by -o
 public:

  /**
  *   @brief flags defined from command line input
  *
  **/ 
  enum GKCCmdLineFlags { GKC_STATISTICS=1, GKC_VERBOSE=2, GKC_OVERWRITE=4, GKC_READ_STDIN=8} ;
 
  /**
  *    @brief Converts a number to std::string
  *
  **/ 
  template<class T> static std::string num2str(T number) 
  {
      std::stringstream ss;
      ss << number;
      return ss.str();
  }

  /**
  *   @brief hold  GKCCmdLineFlags  flags
  *
  *   Hold the corresponding flags as defined from the command
  *   line options
  *
  **/ 
  int flags;

  /**
  *   @brief command line arguments after -x
  *
  *   These command line options are indended to be passed
  *   to 3rd party libraries like PETSc.
  *
  *   @todo : better use direct option string e.g.
  *           PETSc.Option = -ksp
  *
  **/ 
  std::vector<char *> ExArgv;

  /**
  *   @brief command line input parameters
  *
  **/ 
  int argc; char **argv;

  std::string setupFilename   ,  ///< filename of configuration file
              configFileString;  ///< Configuration file as string

  std::string parser_constants;  ///< Constants for function parser 

  /**
  *   @brief Function parser to parse and evaluate mathematical eqution
  *
  *   See  http://warp.povusers.org/FunctionParser/
  *
  **/ 
  FunctionParser    getFParser();

  /**
  *   @brief constructor
  *
  **/ 
  Setup(const int argc, char **argv, std::string setup_filename = "", std::string setup_decomposition= "1:1:1:1:1:1", 
         std::string setup_Xoptions="", std::string setup_ExArgv="",  int flags=0);
  
  /**
  *   @brief destructor
  *
  **/ 
 ~Setup() { };


  /**
  *   @brief erases specific characters from the string 
  **/
  static std::string eraseCharacter(std::string str, std::string chars);

  /**
  *   @brief trim all characters to lower space
  **/
  static std::string trimLower(std::string str, bool lowerCase=true);

  /**
  *   @brief split the string into 2 parst separated by delim
  *
  *   @param str   original string
  *   @param delim delim string
  *
  *   @return vector of string split by delim
  *
  **/
  static std::vector<std::string> split(std::string str, std::string delim);

   /**
   *   @brief parser configuration files
   *
   **/ 
   void parseOption(std::string line, bool fromFile = true);

   /**
   *  @brief checks if all options were red.
   *
   *  After GKC.cpp intializes all sub-modules, this 
   *  subroutine should be called to check if some
   *  options, where not requested but set in the configuration file
   *  (e.g. which may happen due to a spelling error)
   *
   **/ 
   void check_config();

   /** 
   *    @brief Access Elements from Configuration file
   *
   *    Will drop an error in case a key is duplicated.
   *
   **/
   std::string get(std::string key, const char *default_Value) { return get(key, std::string(default_Value)); };
  
   
   /** 
   *    @brief Access Elements from Configuration file
   *
   *    Will drop an error in case a key is duplicated.
   *   
   *    @note : We use explicit template intantiation
   **/
   template<class T> T get(std::string key, const T default_Value);
  
   /**
   *
   * @brief returns int of seconds from time string
   * 
   *  time_string should be of the form e.g.
   *  1y:1l:1w:1d:1h:24m:7s, 1h:24m, 1d:30m
   *
   *  month [l] = fixed to 30 days 
   *              (note as m collides with m for minutes we choose l (luna) for months
   *  one year  = fixed to 365 days
   *
   **/
   static int getSecondsFromTimeString(std::string time_string);
     
  protected:


   /** 
   *   @brief configuration summary output 
   *
   **/ 
   virtual void printOn(std::ostream &output) const;

   // we don't have FileIO object yet, so no writting is performed

};


#endif // __GKC_SETUP_H__
