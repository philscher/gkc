/*
 * =====================================================================================
 *
 *       Filename:  FieldsDST.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  01/18/2012 05:57:35 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef __FIELDS_DST_H
#define __FIELDS_DST_H

#include "Global.h"

#include "Setup.h"
#include "Grid.h"
#include "Parallel.h"
#include "Geometry.h"

#include "FieldsFFT.h"

class FieldsDST : public FieldsFFT {

private:

public:

  //* Constructor
  FieldsDST(Setup *setup, Grid *grid, Parallel *parallel, FileIO * fileIO, Geometry<HELIOS_GEOMETRY> *geo, FFTSolver *fftsolver); 
  //* Destructor
  ~FieldsDST();
  
  Array3z virtual solvePoissonEquation(Array3z rho, Timing timing);
        
  virtual void printOn(ostream &output) const {
        Fields::printOn(output);
         output   << "Fields     |  DST "  << std::endl;
  }

};



#endif // __FIELDS_DST_H

