/*
 * =====================================================================================
 *
 *       Filename: Setup.cpp
 *
 *    Description: Read in Config parameters
 *
 *         Author: Paul P. Hilscher (2009-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */



#include "Setup.h"

#include <iostream>
#include <string>
#include <cstdlib>
#include<vector>

#include "Global.h"
#include "Plasma.h"

#include<sstream>


template<typename T> void PrintVector(const std::vector<T>& t){
        std::cout << std::endl << "Printing vector contents" << std::endl;
        for(typename std::vector<T>::size_type i=0; i<t.size(); ++i){
                 std::cout << t[i] << '\t';
        }
        std::cout << std::endl << std::endl;
}




/** Note :
* string
* The position of the first occurrence in the string of the searched content.
*If the content is not found, the member value npos is returned.
*
* Input string original string
  @2 chars : chatacters to replace
*/

std::string eraseCharacter(std::string str, std::string chars) {

    for(int p, s = 0; s < chars.length(); s++) { 
    	while((p = str.find(chars[s])) != std::string::npos) str.erase(p,1);
    }
/*
    	while((str.find(chars[s])) != std::string::npos) {
		p = str.find(chars[s]);
		str.erase(p,1);
    } 
*/
    return str;
};




  // std::string trim adapted from 
  std::string Setup::trimLower(std::string str, bool lowerCase)  { 

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

  std::vector<std::string> Setup::split(std::string str, std::string delim) {
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
  }

  Setup::Setup(int _argc, char **_argv, std::string setup_filename, std::string setup_decomposition, std::string setup_Xoptions, std::string setup_ExArgv, int process_id, int _flags) :
   flags(_flags)
  {
    commandLineOptions = setup_Xoptions;
      setupFilename = setup_filename;
      ///////////////////////////////////////////////////////////////////////////////////////
      // Create own, see http://stackoverflow.com/questions/1511797/convert-string-to-argv-in-c
    std::istringstream iss(setup_ExArgv);
    std::string token;
    
          ExArgv.push_back("helios");
    while(iss >> token) {
          char *arg = new char[token.size() + 1];
          copy(token.begin(), token.end(), arg);
          arg[token.size()] = '\0';
          ExArgv.push_back(arg);
    }
    ExArgv.push_back(0);
    //argc = _argc;//ExArgv.size()-1;
    //argv = _argv;//&ExArgv[0];
    argc = ExArgv.size()-1;
    argv = &ExArgv[0];


    // now exec with &args[0], and then:
    //
     
      
      
      
      config["Parallel.Decomposition"] = setup_decomposition;
      config["Helios.Process_ID"]     =  number2string(process_id);
     

       //////////////////// parse Setup file /////////////////
      if(setupFilename == "") check(-1, DMESG("No config file provided. Start helios with : -c <config_file_name> option"));
      std::string line;
      std::ifstream file;
      file.open(setupFilename.c_str(), ios::in);
      check(file.is_open(), DMESG("Could open find file! Wrong path ?"));

      // read file in as a string
      std::stringstream buffer;
      buffer << file.rdbuf();
      configFileString = buffer.str();

      // now parse it
      file.seekg (0, ios::beg);
      // read setup file line per line
      while(std::getline(file, line)) parseOption(line);
      
      file.close();
      
      if(flags & HELIOS_READ_STDIN) {
       // while(std::getline(file, line)) parseOption(line);
        
      }
      // parse options from command line input
      if(setup_Xoptions != "")  {
           std::vector<std::string> options_list = split(setup_Xoptions, ";");
           while(options_list.empty() == false) { 
             parseOption(options_list.back(), false);
             options_list.pop_back();
           }
      }


     // ToDo : check some values
    
     parser_constants = eraseCharacter(get("Setup.Constants", ""), " {}");

  }



 int Setup::parseOption(std::string line, bool fromFile) {

               if((line[0] == '#') || (line[0] == '!')) return HELIOS_SUCCESS;
               if(line.empty() == true) return HELIOS_SUCCESS;
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


           return HELIOS_SUCCESS;
}


FunctionParser Setup::getFParser() {

    FunctionParser parser;

    parser.AddConstant("pi", M_PI);
    parser.AddConstant("Lx", Lx); parser.AddConstant("Nx", (double) Nx);
    parser.AddConstant("Ly", Ly); parser.AddConstant("Nky", (double) Nky);
    parser.AddConstant("Lz", Lz); parser.AddConstant("Nz", (double) Nz);
    parser.AddConstant("Lv", Lv); parser.AddConstant("Nv", (double) Nv);
    parser.AddConstant("Lm", Lm); parser.AddConstant("Nm", (double) Nm);
    parser.AddConstant("Ns", (double) Ns);
//    parser.AddConstant("shat", geometry->shear);

    // we can define in values as Setup.Parser = { eps = 0.01; sigma = 0.1 }, we need to parse it
    // split in eps = 0.01

    std::vector<std::string> const_vec = split(parser_constants, ",");
    for(int s = 0; s < const_vec.size(); s++) { 
            std::vector<std::string> key_value = split(const_vec[s],"=");
            parser.AddConstant(trimLower(key_value[0], false), string_to_double(key_value[1]));
    };
    
    return parser;

 
}
