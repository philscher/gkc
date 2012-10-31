/*
 * =====================================================================================
 *
 *       Filename:  Input.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/09/2011 10:59:50 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef __INPUT_H
#define __INPUT_H

#include <sstream>
#include <string>
#include <cctype>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
 #include <sstream>


#include "FunctionParser/fparser.hh"

class Input {
  std::map   <std::string, std::string> config;
  std::vector<std::string> config_check;
  public:
    template<class T> static std::string num2str(T num);
   int check_config();
   
   
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
            std::vector<std::string>::iterator f = find(config_check.begin(), config_check.end(), key);
            if( f != config_check.end() ) config_check.erase(f);
            
            
            return t;
        }
        else if(config.count(key) > 1) std::cout << "Parser Error elements occurse more than once" << std::endl;



        return default_Value;
    };
  
  // std::string trim adapted from 
  std::string static trimLower(std::string str, bool lowerCase)  { 

       // Trim Both leading and trailing spaces  
       size_t startpos = str.find_first_not_of("\t "); // Find the first character position after excluding leading blank spaces  
       size_t endpos = str.find_last_not_of("\t "); // Find the first character position from reverse af  
     
       // if all spaces or empty return an empty string  
       if(( std::string::npos == startpos ) || ( std::string::npos == endpos))  
       {  
           str = std::string("");  
       }  
      else  
           str = str.substr( startpos, endpos-startpos+1);  
    
       // Lower case
      if(lowerCase == true) for(unsigned int i=0;i<str.length();i++) str[i] = std::tolower(str[i]);
       
       return str; 
  }    
 
  int parseOption(std::string line, bool fromFile=true);
  Input(std::string setupFilename, std::string setup_Xoptions);

std::vector<std::string> split(std::string str, std::string delim);
FunctionParser getFParser();
    
};

#endif // __INPUT_H
