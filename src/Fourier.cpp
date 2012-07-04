/*
 * =====================================================================================
 *
 *       Filename: Fourier3D.cpp
 *
 *    Description: Helper functions for Fourier space
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include "Fourier.h"

#include <cmath>
#include <cfloat>
#include<iostream>
#include<vector>
#include<string>
#include<random/uniform.h>



void parseSuppressMode(const std::string &value, std::vector<int> &suppressMode) {
    if(value == "") return;

    //int start_s = value.find_first_of("(");
    //int   end_s = value.find_first_of(")",start_s+1);
    //std::string sub_str = (value.substr(start_s+1,end_s-1));
    std::string sub_str = Setup::eraseCharacter(value, " ") ;
    std::vector<std::string> result = Setup::split(sub_str, ",");
                                              
    vector<std::string>::const_iterator it;
                                               for(it = result.begin(); it != result.end(); it++) {
                                                   if((*it).find("-") == string::npos) {
                                                   //if((*it).find("-") != -1) {
                                                    std::vector<std::string> res2 = Setup::split(*it, "-");
                                                    for(int i = atoi(res2[0].c_str()); i < atoi(res2[1].c_str()); i++) {
                                                        suppressMode.push_back(i);
                                                    }} else
                                                    suppressMode.push_back(std::atoi((*it).c_str()));
                                                }
              
              }


Fourier3D::Fourier3D(Setup *setup, Grid *grid, FFTSolver *_fftsolver)  : fft(_fftsolver)
{
       parseSuppressMode(setup->get("SuppressModeX", ""), suppressModeX);  
       parseSuppressMode(setup->get("SuppressModeY", ""), suppressModeY);  
       parseSuppressMode(setup->get("SuppressModeZ", ""), suppressModeZ);  

} 

Fourier3D::~Fourier3D() {}

// BUG : Input SuppressModeX = 0 does not work (Crashes)
int Fourier3D::suppressModes(Array4z A, const int field) 
{

  // suppress x,y, z - Modes
  vector<int>::const_iterator mode;

  for(mode = suppressModeX.begin(); (!suppressModeX.empty()) && (mode <= suppressModeX.end()) ; mode++) { 
	if((*mode < fft->K1xLlD) || (*mode > fft->K1xLuD)) continue;
    //std::cout << " Mode  X : " << (*mode) << std::endl << std::flush;
     A(*mode, RkyLD, RzLD, RFields) = 1.e-128;
  }
//  for(mode = suppressModeY.begin(); mode <= suppressModeY.end(); mode++) {
  for(mode = suppressModeY.begin(); (!suppressModeY.empty()) && (mode <= suppressModeY.end()) ; mode++) { 
	    if((*mode < NkyLlD) || (*mode > NkyLuD)) continue;
    //std::cout << " Mode  Y : " << (*mode) << std::endl << std::flush;
        A(fft->Rk1xL, *mode, RzLD, RFields) = 1.e-128;
  }


  return HELIOS_SUCCESS;
}


/*
// we should investigate the effect of the filter in more detail... !!
// should we have a user modificatable filter, e.g. using a mapping ?
double Fourier3D::filter(int x_k, int y_k) {

// Lets use Dorland Filtering, which means multiplications with exp(-k_p square)


// wrong filtering 
// BUG, use sinh(x/0.4)
   const double kx = fft->kx(x_k);
   const double ky = fft->ky(y_k);
   return      ((kx != 0.) ? pow3(tanh(0.4*kx)/kx) : 1.) * ((ky != 0.) ? pow3(tanh(0.4*ky)/ky) : 1.);

}
*/



/**

\delta n_{0e} = n_0 \exp_0 \left( - \frac{(x/L_x - \tfrac{1}{2})^2 + (y/L_y - \tfrac{1}{2})^2 + (z/L_z - \tfrac{1}{2} ) ^2} {2 \sigma^2} \right)

*/

