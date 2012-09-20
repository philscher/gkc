/*
 * =====================================================================================
 *
 *       Filename:  ReaderXYV.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  07/19/2011 03:06:30 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef READER_XYV_H
#define READER_XYV_H

#include <iostream>
#include <fstream>

#include "Global.h"
#include "Reader/Reader.h"


#include "Special/Interpolate/Interpolate.h"


/**
 *     Reads in structure 
 *     X
 *     Y
 *     Data (Y=0)
 *     Data (Y=1)
 *     Data (Y=2)
 *     .....
 *
 *     Pretty picky, only one white space allowed between subsequent numbers !  
 *
 *
 *
 * */



class ReaderXYV : public Reader {

    Array1R X;
    Array1R Y;

    Interpolate *interpol;

    Array1R parseNumberLine(std::string line, const char delimiter=' ') {

       std::vector<std::string> token = Setup::split(Setup::trimLower(line), &delimiter);
           
       Array1R V(blitz::Range(0,token.size()-1));


       for(int n = 0; n < token.size();n++)  V(n) = Setup::string_to_double(token[n]);

      return V;
    };


public:
    ReaderXYV(std::string fileName, Setup *setup) : Reader(setup) {
        
       ifstream inputFile(fileName.c_str(), ifstream::in);
   
       if(!inputFile.is_open()) check(-1, DMESG("Coudl not open Reader file, Path correct ?"));

       std::string line; 
       // get first line count
       int lc=0;
       while( getline(inputFile, line) ) lc++;

       // rewind and analyis;
       inputFile.clear() ; 
       inputFile.seekg(0) ;

       getline(inputFile, line);
       X.reference(parseNumberLine(line));
       getline(inputFile, line);
       Y.reference(parseNumberLine(line));

       std::cout << " Counted data lines : " << lc -2 << " Y size : " << Y.numElements() << std::endl;
       // Simple check for consistency, Y len is the number of data rows
       if((lc-2) != Y.numElements()) check(-1, DMESG("Number of lines and size in Y does not match"));

       Array2R V(blitz::Range(0,X.numElements()-1), blitz::Range(0, Y.numElements()-1));

       // Read in main data file
       for(int ny = 0; ny < Y.numElements(); ny++) {

            getline(inputFile, line);
            V(blitz::Range::all(), ny) = parseNumberLine(line);

       }

      inputFile.close();
      interpol = new Interpolate(X,Y,V, setup, setup->get("Reader.Interpolate", "Linear"));

    };
   
    double getValue(const double x, const double y, const double z) {
        
        if(z != 0.) check(-1, DMESG("Reader allowes only 2-dimensional values"));

        return interpol->getValue(x,y,z);


    };

   



};




#endif // READERXYV_H

