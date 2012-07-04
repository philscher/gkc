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

#ifndef SETUP_H__
#define SETUP_H__


// define to compile Fparser with complex number support
#define FP_SUPPORT_COMPLEX_DOUBLE_TYPE

#include <string>
#include <cctype>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
 #include <sstream>

#include "Global.h"
#include "FunctionParser/fparser.hh"

enum SpecRange { SPEC_START = 0, SPEC_STRIDE=1, SPEC_END=2 };
enum SpecDir   {SPEC_NO=-1, SPEC_XY=0, SPEC_XZ=1, SPEC_YZ=2};
enum Decomposition { DECOMP_NO = 0, DECOMP_X=1, DECOMP_Y=2, DECOMP_XY=3, DECOMP_Z=4, DECOMP_XYZ=7,DECOMP_V=8, DECOMP_M=16, DECOMP_S=32};
enum VlasovSolverT { VL_NO=0, VL_LIN=1, VL_NONLIN=2, VL_GEO=4, VL_TRAP=8};
enum HeliosFlagsT { HELIOS_STATISTICS=1, HELIOS_VERBOSE=2, HELIOS_OVERWRITE=4, HELIOS_READ_STDIN=8};


class Setup : public IfaceHelios {

    static double sign(double T) { return ((T >= 0.) ? 1. : -1); }

    map   <std::string, std::string> config;
    vector<std::string> config_check;

    std::string commandLineOptions, extraLineOptions;
public:
 
 /** create a stringstream
        add number to the stream
        return a string with the contents of the stream
*/
template<class T> static std::string number2string(T number)
{
     stringstream ss;
     ss << number;
     return ss.str();
};
 
static double string_to_double( const std::string s )
   {
        return atof(s.c_str());
        /* 
        std::istringstream i(s);
        double x;
        if (!(i >> x)) return 0;
        return x;
         * */
   }

    int flags;

    std::vector<char *> ExArgv;
    int argc; char **argv;
    
    std::string setupFilename, configFileString;
    std::string parser_constants;

    FunctionParser    getFParser();
//    FunctionParser_cd getFParser_cd();
   // Some user function

   Setup(const int argc, char **argv, std::string setup_filename = "", std::string setup_decomposition= "1:1:1:1:1:1", std::string setup_Xoptions="", std::string setup_ExArgv="", int process_id=0, int flags=0);
  ~Setup() { };




   static std::string eraseCharacter(std::string str, std::string chars);
   static std::string trimLower(std::string str, bool lowerCase=true);
   static std::vector<std::string> split(std::string str, std::string delim);
   int parseOption(std::string line, bool fromFile = true);

   int check_config() {
       if(config_check.size() != 0) {
            vector<std::string>::const_iterator mode;
            for(mode =  config_check.begin(); mode != config_check.end(); mode++) std::cout << "Element not accessed : " << *mode << std::endl;
            check(-1, DMESG("Parsing Error, Elements not accessed"));
       }

       return HELIOS_SUCCESS;

   };

  /** 
   *    Access Elements from Configuration file
   *
   *
   *
   *
   **/
  std::string get(std::string key, const char *default_Value) { return get(key, std::string(default_Value)); };
  template<class T> T get(std::string key, const T default_Value)
    {
        if(     config.count(key) == 1) {
            std::string s = config[key];
            std::istringstream stream (s);
            T t;
            stream >> t;
           
            // delete element in check
            vector<std::string>::iterator f = find(config_check.begin(), config_check.end(), key);
            if( f != config_check.end() ) config_check.erase(f);
            
            
            return t;
        }
        else if(config.count(key) > 1) check(-1, DMESG("Parser Error elements occurse more than once"));



        return default_Value;
    };
    
     /**
     *  time_string should be of the form 1y:1m:1w:1d:1h:24m:7s, 1h:24m, 1d:30m
     *                                   while years are nonsense it 
     *  month [l] = 30 days (note as m collides with m for minutes whe choose l for luna)
     *  year  = 365 days
     *
     * */
    static int getSecondsFromTimeString(std::string time_string) {

        std::vector<std::string> const_vec = split(time_string, ":");
        int seconds = 0;

        for(unsigned int s = 0; s < const_vec.size(); s++) { 
            
            std::string token = const_vec[s];

            int c2sec = 1;
            char id = token[token.length()-1];
            
            switch(id) {
                case('s') : c2sec = 1 ; break;
                case('m') : c2sec = 60; break; 
                case('h') : c2sec = 3600; break; 
                case('d') : c2sec = 86400; break; 
                case('w') : c2sec = 604800; break; 
                case('l') : c2sec = 2592000;  break;
                case('a') : c2sec = 3153600; break;
                default : check(-1, DMESG("No such token for time string"));
            };
            seconds += c2sec*atoi(token.substr(0, token.length()-1).c_str());
        };
    
        return seconds;
    };
     
protected:
   virtual void printOn(ostream &output) const {
         output << "Parser     | Parsing Config  File" << std::endl;
         if(commandLineOptions != "")   output << "Parser   | with " << commandLineOptions << std::endl;
         if(extraLineOptions != "")   output << "Extra   | with " << commandLineOptions << std::endl;
	 if(parser_constants != "")   output <<   "         | FParser Constants : " << parser_constants << std::endl; 
    
 }

     // we don't have FileIO object yet, so no writting is performed
 //    virtual void initDataOutput(FileIO *fileIO) {};
 //    virtual void writeData(Timing *timing) {};
 //    virtual void closeData() {};

};


#endif // SETUP_H__
