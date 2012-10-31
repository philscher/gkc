/*
 * =====================================================================================
 *
 *       Filename:  Input.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/09/2011 10:57:42 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include<string>
#include <iostream>
#include <fstream>

#include "Input.h"

inline int check( int status, std::string file, int line, std::string error_text, bool doAbort=false) {
        if(status == -1 ) {
            
       // check rank from MPI!
            std::stringstream ss;
            ss << std::endl; 
            ss << "\033[0;m";
            ss << "\033[1;m"  <<"!.... " << file << "(" << line << ") " << error_text; 
            ss << "\033[0m" << std::endl << std::flush;

            std::cout << ss.str();
              
            // exit through abort so we can get stack trace
       if(doAbort==true) abort();
            abort();
            exit(0);
        }
        return status;
    }


#define DMESG(mesg)  std::string(__FILE__), __LINE__, std::string(mesg)










    template<class T> static std::string num2str(T num) {
      std::stringstream ss;
        ss << num;
        return ss.str();
    };




   int Input::check_config() {

       if(config_check.size() != 0) {
            std::vector<std::string>::const_iterator mode;
            for(mode =  config_check.begin(); mode != config_check.end(); mode++) std::cout << "Element not accessed : " << *mode << std::endl;
            check(-1, DMESG("Parsing Error, Elements not accessed"));
       }

       return 1;

   };

  /** 
   *    Access Elements from Configuration file
   *
   *
   *
   *
   **/
 
  int Input::parseOption(std::string line, bool fromFile) {

               if((line[0] == '#') || (line[0] == '!')) return 1;
      // read file in as a string
               if(line.empty() == true) return 1;
              int posEqual=line.find('=');

                std::string name  = trimLower(line.substr(0,posEqual), false);
                std::string value = trimLower(line.substr(posEqual+1), false);

                 std::map<std::string, std::string> :: const_iterator el;
                 el = config.find(value);

                 if((el != config.end()) && fromFile) { 
                     std::cout << "Parsing Error : Duplicated key in File " << name << std::endl; check(-1, DMESG("Duplicate key added"));
                 }
                 if(fromFile) config_check.push_back(name) ;
                 config[name] = value;


           return 1;
}

Input::Input(std::string setupFilename, std::string setup_Xoptions) 
  {
      ///////////////////////////////////////////////////////////////////////////////////////
      // Create own, see http://stackoverflow.com/questions/1511797/convert-string-to-argv-in-c
    // now exec with &args[0], and then:
    //
       //////////////////// parse Setup file /////////////////
      if(setupFilename == "") check(-1, DMESG("No config file provided. Start helios with : -c <config_file_name> option"));
      std::string line;
      std::ifstream file;
      file.open(setupFilename.c_str(), std::ios::in);
      check(file.is_open(), DMESG("Could open find file! Wrong path ?"));

      file.seekg (0, std::ios::beg);
      // read setup file line per line
      while(std::getline(file, line)) parseOption(line);
      
      file.close();
      
      if(setup_Xoptions != "")  {
           std::vector<std::string> options_list = split(setup_Xoptions, ";");
           while(options_list.empty() == false) { 
             parseOption(options_list.back(), false);
             options_list.pop_back();
           }
      }



}
  
std::vector<std::string> Input::split(std::string str, std::string delim) {
      std::vector<std::string> items;
        std::size_t dlm_idx;
          if(str.npos == (dlm_idx = str.find_first_of(delim))) {
                items.push_back(str.substr(0, dlm_idx));
                  }
            while(str.npos != (dlm_idx = str.find_first_of(delim))) {
                  if(str.npos == str.find_first_not_of(delim)) {
                          break;
                              }
                      items.push_back(str.substr(0, dlm_idx));
                          dlm_idx++;
                              str = str.erase(0, dlm_idx);
                                  if(str.npos == str.find_first_of(delim) && "" != str) {
                                          items.push_back(str);
                                                break;
                                                    }
                                    }
              return items;
  };

FunctionParser Input::getFParser() {

    FunctionParser parser;

    parser.AddConstant("pi", M_PI);
//    parser.AddConstant("Lx", Lx); parser.AddConstant("Nx", (double) Nx);
//    parser.AddConstant("shat", geometry->shear);

    // we can define in values as Setup.Parser = { eps = 0.01; sigma = 0.1 }, we need to parse it
    // split in eps = 0.01
/* 
    std::vector<std::string> const_vec = split(parser_constants, ",");
    for(int s = 0; s < const_vec.size(); s++) { 
            std::vector<std::string> key_value = split(const_vec[s],"=");
	        //std::cout << " Key : " <<   trimLower(key_value[0], false) << " value : " <<  string_to_double(key_value[1]) << std::endl;
            parser.AddConstant(trimLower(key_value[0], false), string_to_double(key_value[1]));
    };
 * */
    
    return parser;

 
}
