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

    int start_s = value.find_first_of("(");
    int   end_s = value.find_first_of(")",start_s+1);
    std::string sub_str = (value.substr(start_s+1,end_s-1));
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
       
       parseSuppressMode(setup->get("ConvolveModeX", ""), convolveModeX);  
       parseSuppressMode(setup->get("ConvolveModeY", ""), convolveModeY);  
       parseSuppressMode(setup->get("ConvolveModeZ", ""), convolveModeZ);  
       
       epsilon_0 = setup->get("Init.Epsilon0", 1.e-10); 
       sigma     = setup->get("Init.Sigma"   , 1.e-1); 

} 

Fourier3D::~Fourier3D() {}


int Fourier3D::suppressModes(Array4z phi_k, const int field) 
{

  // suppress x,y, z - Modes
  vector<int>::const_iterator mode;
  for(mode = suppressModeX.begin(); mode != suppressModeX.end(); mode++) { 
	if((*mode < fft->K2xLlD) || (*mode > fft->K2xLuD)) continue;
     phi_k(*mode, Range::all(), Range::all(), 1) = 1.e-30;
  }
  for(mode = suppressModeY.begin(); mode != suppressModeY.end(); mode++) {
	    if((*mode < fft->K2yLlD) || (*mode > fft->K2yLuD)) continue;
        phi_k(Range::all(), *mode, Range::all(), 1) = 1.e-30;
  }


  return HELIOS_SUCCESS;
}

/*
double sinConv(int x, int q, int Lp, int Up, int mode) {

          const int p = (x-q < Lp) ? Up - q + 1 : x - q; 
             
          return  sin(mode * 2. * M_PI * X(p)/Lx)/((double) Nx/2.);
}
 * \ */


int Fourier3D::suppress3DMode(Array4d A) {
    
      // this is dummy for flux average
   if((fft->getFlags() & FFT_XYZ) == false) return HELIOS_SUCCESS;

   //for(int nfield=1;nfield <= plasma->nfields; nfield++) { 
   for(int nfield=1;nfield <= plasma->nfields; nfield++) { 

    fft->r33In(RxLD, RyLD, RzLD) = fft->r3Out(RxLD, RyLD, RzLD, nfield);
    fft->solve(FFT_XYZ, FFT_FORWARD);
    fft->k3In = fft->k3Out/fft->Norm_XYZ;
  
    // For z we use convolution tp suppress unwanted modes
    vector<int>::const_iterator mode;
  
    for(mode = convolveModeX.begin(); mode != convolveModeX.end(); mode++) { 
	if((*mode < fft->K3xLlD) || (*mode > fft->K3xLuD) || (*mode > Nx/2)) continue;
        fft->k3In(*mode, fft->Rk3yL, fft->Rk3zL) = 0.e0;
        //fft->k3In(*mode, fft->Rk3yL, fft->Rk3zL) = DBL_MIN;
        // suppress conjugate mode
        int cnjg_mode = (Nx - *mode);
	if((cnjg_mode < fft->K3xLlD) || (cnjg_mode > fft->K3xLuD) || (cnjg_mode == 0)) continue;
        //fft->k3In(cnjg_mode, fft->Rk3yL, fft->Rk3zL) = DBL_MIN;
        fft->k3In(cnjg_mode, fft->Rk3yL, fft->Rk3zL) = 0.e0;

    }

    for(mode = convolveModeY.begin(); mode != convolveModeY.end(); mode++) { 
	
        if((*mode < fft->K3yLlD) || (*mode > fft->K3yLuD) || (*mode > Nky-1)) continue;
        fft->k3In(fft->Rk3xL, *mode, fft->Rk3zL) = 0.e0;
        //fft->k3In(fft->Rk3xL, *mode, fft->Rk3zL) = DBL_MIN;
    
        // suppress conjugate mode
        //int cnjg_mode = (Nky - *mode);
	//if((cnjg_mode < fft->K3yLlD) || (cnjg_mode > fft->K3yLuD) || (cnjg_mode == 0)) continue;
        //fft->k3In(fft->Rk3xL, cnjg_mode, fft->Rk3zL) = 0.e0;
    }
  
    for(mode = convolveModeZ.begin(); mode != convolveModeZ.end(); mode++) { 
	if((*mode < fft->K3zLlD) || (*mode > fft->K3zLuD) || (*mode > Nz/2)) continue;
        fft->k3In(fft->Rk3xL, fft->Rk3yL, *mode) = 0.e0;
        // suppress conjugate mode
        int cnjg_mode = (Nz - *mode);
	if((cnjg_mode < fft->K3zLlD) || (cnjg_mode > fft->K3zLuD) || (cnjg_mode == 0)) continue;
        fft->k3In(fft->Rk3xL, fft->Rk3yL, cnjg_mode) = 0.e0;
    }

    fft->solve(FFT_XYZ, FFT_BACKWARD);
    fft->r3Out(RxLD, RyLD, RzLD, nfield) = fft->r33Out(RxLD, RyLD, RzLD);
   
   } 
  /* 
  // X
  for(mode = convolveModeX.begin(); mode != convolveModeX.end(); mode++) {
          fft->r3In = fft->r3Out; fft->r3Out = 0.;
          
          for(int x=NxGlD; x <= NxGuD; x++) {
       
            
                for(int q=NxGlD-3; q <= NxGuD-3; q++) fft->r3Out(x, RyLD, RzLD) += fft->r3In(x, RyLD, RzLD) *  sinConv(x,q, NxGlD, NxGuD, *mode);
          //for(int q=NxGlD; q <= NxGuD; q++) fft->r3Out(x, RyLD, RzLD) += fft->r3In((x-q < NxGlD) ? NxGuD - q + 1 : x - q, RyLD, RzLD) *  sin(*mode * 2. * M_PI * X(q)/Lx) / ((double) Nx);

          }
        }


  // For z we use convolution tp convolve unwanted modes
  
  for(mode = convolveModeZ.begin(); mode != convolveModeZ.end(); mode++) {
          fft->r3In = fft->r3Out;
          
          for(int z=NzGlD; z <= NzGuD; z++) {
          fft->r3Out(RxLD, RyLD, z) = 0.;
          
          for(int q=NzGlD; q <= NzGuD; q++) fft->r3Out(RxLD, RyLD, z) += fft->r3In(RxLD, RyLD, (z-q < NzGlD) ? NzGuD - q + 1 : z- q) *  sin(*mode * 2. * M_PI * Z(q)/Lz) / ((double) Nz);


        }
  }
  * */ 
  
  
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

