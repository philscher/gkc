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
#include <cstdlib>

#include "Plasma.h"


// Debuging Helper function (still needed ?) depraced
template<typename T> void PrintVector(const std::vector<T>& t)
{
   
  std::cout << std::endl << "Printing vector contents" << std::endl;
  
  for(typename std::vector<T>::size_type i=0; i<t.size(); ++i){
    std::cout << t[i] << '\t';
  }
  std::cout << std::endl << std::endl;
}




std::string Setup::eraseCharacter(std::string str, std::string chars) 
{

   for(int p, s = 0; s < chars.length(); s++) { 
      while((p = str.find(chars[s])) != std::string::npos) str.erase(p,1);
   }

    return str;
};

std::string Setup::trimLower(std::string str, bool lowerCase)  
{ 

   
  // Trim Both leading and trailing spaces  
  size_t startpos = str.find_first_not_of("\t "); // Find the first character position after excluding leading blank spaces  
  size_t endpos = str.find_last_not_of("\t "); // Find the first character position from reverse af  
     
  // if all spaces or empty return an empty string  
  if(( std::string::npos == startpos ) || ( std::string::npos == endpos))   str = std::string("");  
  else  str = str.substr( startpos, endpos-startpos+1);  
    
  // Lower case
  if(lowerCase == true) for(unsigned int i=0;i<str.length();i++) str[i] = std::tolower(str[i]);
       
  return str; 

}    


std::vector<std::string> Setup::split(std::string str, std::string delim) 
{
    
   std::vector<std::string> items;
   std::size_t dlm_idx;
   
   if(str.npos == (dlm_idx = str.find_first_of(delim))) {
     items.push_back(str.substr(0, dlm_idx));
   }
   
   while(str.npos != (dlm_idx = str.find_first_of(delim))) {
   
     if(str.npos == str.find_first_not_of(delim)) break;
     
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
  extraLineOptions   = setup_ExArgv;
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

  
  config["Parallel.Decomposition"] = setup_decomposition;
  config["Helios.Process_ID"]      =  number2string(process_id);
     

  
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
      
  
  if(flags & GKC_READ_STDIN)  {
    while(std::getline(file, line)) parseOption(line);
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
    
  parser_constants = eraseCharacter(get("Setup.Constants", ""), "{}");

}



int Setup::parseOption(std::string line, bool fromFile) 
{

   
  if((line[0] == '#') || (line[0] == '!')) return GKC_SUCCESS;
  if(line.empty() == true) return GKC_SUCCESS;
  
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
  
  return GKC_SUCCESS;
}


FunctionParser Setup::getFParser() 
{

   FunctionParser parser;

   parser.AddConstant("pi", M_PI);
   parser.AddConstant("Lx", Lx); parser.AddConstant("Nx", (double) Nx);
   parser.AddConstant("Ly", Ly); parser.AddConstant("Nky", (double) Nky);
   parser.AddConstant("Lz", Lz); parser.AddConstant("Nz", (double) Nz);
   parser.AddConstant("Lv", Lv); parser.AddConstant("Nv", (double) Nv);
   parser.AddConstant("Lm", Lm); parser.AddConstant("Nm", (double) Nm);
   parser.AddConstant("Ns", (double) Ns);
    
   //  BUG : Crashes if parser_constants is empty. Why ?
   
   if (!parser_constants.empty()) {

      std::vector<std::string> const_vec = split(parser_constants, ",");
      
      for(int s = 0; s < const_vec.size(); s++) { 

            std::vector<std::string> key_value = split(const_vec[s],"=");
            parser.AddConstant(trimLower(key_value[0], false), string_to_double(key_value[1]));

      };
   }
   
   return parser;

}


int Setup::getSecondsFromTimeString(std::string time_string) 
{
   
  std::vector<std::string> const_vec = split(time_string, ":");
  int seconds = 0;

  for(unsigned int s = 0; s < const_vec.size(); s++) { 
            
    std::string token = const_vec[s];

    
    int c2sec = 1;
    char id   = token[token.length()-1];
            
    
    switch(id) {
    
      case('s') : c2sec = 1       ; break;
      case('m') : c2sec = 60      ; break; 
      case('h') : c2sec = 3600    ; break; 
      case('d') : c2sec = 86400   ; break; 
      case('w') : c2sec = 604800  ; break; 
      case('l') : c2sec = 2592000 ; break;
      case('a') : c2sec = 31536000; break;
      
      default : check(-1, DMESG("No such token for time string"));
      
    };
    
    seconds += c2sec*atoi(token.substr(0, token.length()-1).c_str());
    
  };
    
        return seconds;
};
   
int Setup::check_config() 
{
   
   if(config_check.size() != 0) {
  
      vector<std::string>::const_iterator mode;
    
      for(mode =  config_check.begin(); mode != config_check.end(); mode++) { 
            std::cout << "Element not accessed : " << *mode << std::endl;
      }
      check(-1, DMESG("Parsing Error, Elements not accessed"));
    
   }
  
   return GKC_SUCCESS;

};
   
void Setup::printOn(ostream &output) const 
{
  output << "Parser     | Parsing Config  File" << std::endl;
  
  if(commandLineOptions != "")   output << "Parser   | with " << commandLineOptions << std::endl;
  if(extraLineOptions   != "")   output << "Extra    | with " << commandLineOptions << std::endl;
  if(parser_constants   != "")   output << "         | FParser Constants : " << parser_constants << std::endl; 
}
