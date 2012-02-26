/*
 * =====================================================================================
 *
 *       Filename:  ReaderXV.h
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


#include "Interpolate/Interpolate1D.h"


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



class ReaderXV : public Reader {

    Array1d X;

    Interpolate *interpol;


public:
    ReaderXV(std::string fileName, Setup *setup) : Reader(setup) {
        
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

       std::cout << " Counted data lines : " << lc -2 << " Y size : " << Y.numElements() << std::endl;
       // Simple check for consistency, Y len is the number of data rows
       if((lc-2) != Y.numElements()) check(-1, DMESG("Number of lines and size in Y does not match"));

       Array1d V(Range(0,X.numElements()-1));

       // Read in main data file
       for(int ny = 0; ny < Y.numElements(); ny++) {

            getline(inputFile, line);
            V(Range::all(), ny) = parseNumberLine(line);

       }

      inputFile.close();
      interpol = new Interpolate1D(X,V, setup, setup->get("Reader.Interpolate", "Linear"));

    };
   
    double getValue(const double x, const double y=0, const double z=0) {
        
        if((z != 0.) || (y != 0)) check(-1, DMESG("Reader allowes only 2-dimensional values"));

        return interpol->getValue(x);


    };

   



};




#endif // READERXV_H

