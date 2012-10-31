/*
 * =====================================================================================
 *
 *       Filename:  FileOutput.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/09/2011 10:57:37 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */


#ifndef __DATAOUTPUT_H
#define __DATAOUTPUT_H

#include <hdf5.h>
#include <hdf5_hl.h>

#include <string>

#include "FileAttr.h"
#include "Input.h"

typedef _Complex double cmplxd;

//typedef std::complex<double> cmplxd;

class DataOutput {
double OutputTimeField, OutputTimePower;
  // Data Output stuff
  hid_t fileOutput;
  FileAttr *FA_Psi_Time, *FA_Psi, *FA_Phi, *FA_PowerKy_Phi, *FA_PowerKy_Psi, *FA_PowerKy_Time;
  
  public:

   DataOutput(Input *input, int Nx, int Nky, double Lx, double Ly, double *X, double *Psi0, double Viscosity, double Resistivity);
   ~DataOutput();
   
   void output(int Step, double Time, double dt, int Nky, int Nx, cmplxd Phi[Nky][Nx], cmplxd Psi[Nky][Nx]);

};


#endif // __DATAOUTPUT_H
