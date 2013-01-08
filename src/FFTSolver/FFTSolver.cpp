/*
 * =====================================================================================
 *
 *       Filename: FFTSolver.cpp
 *
 *    Description: Interface for various (Fast) Fourier Solvers
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include "FFTSolver/FFTSolver.h"

int FFTSolver::X_NkxL;


void FFTSolver::setNormalizationConstants() {


  /////////////////   Find normalization constants for Y-transformation    //////////
  double   rY[NyLD ][NxLD]; 
  CComplex kY[NkyLD][NxLD]; 


  // Real -> Complex (Fourier space) transform
  rY[:][:] = 1.; 
  solve(FFT_Type::Y_NL, FFT_Sign::Forward, rY, kY);
  Norm_Y_Forward = creal(kY[0][0]);
  
  // Complex (Fourier sapce) -> Real transform
  kY[:][:] = 0. ; kY[0][0] = 1.;
  solve(FFT_Type::Y_NL, FFT_Sign::Backward, kY, rY);
  Norm_Y_Backward = rY[0][0];

  /////////////////   Find normalization constants for X-transformation    //////////
  [&](CComplex kXIn [1][NzLD][Nky][FFTSolver::X_NkxL],
      CComplex kXOut[1][NzLD][Nky][FFTSolver::X_NkxL])
  {
     CComplex rXIn  [1][NzLD][Nky][NxLD];

     // Real -> Complex (Fourier sapce) transform
     rXIn[:][:][:][:] = 0.;
     rXIn[:][:][0][:] = 1.;
     solve(FFT_Type::X_FIELDS, FFT_Sign::Forward, ((CComplex *) &rXIn[0][0][0][0]));

     Norm_X_Forward = (K1xLlD == 0) ? creal(kXOut[0][NzLlD][0][0]) : 0.;

     Norm_X_Backward = Norm_X / Norm_X_Forward;
    
  } ((A4zz) kXIn, (A4zz) kXOut); 
  
  // broadcast normalization to all nodes
  //parallel->send(Norm_Y_Forward, parallel->Coord[DIR_ALL] == 0); parallel->send(Norm_Y_Backward, parallel->Coord[DIR_ALL] == 0);
  //parallel->send(Norm_X_Forward, parallel->Coord[DIR_ALL] == 0); parallel->send(Norm_X_Backward, parallel->Coord[DIR_ALL] == 0);
  //parallel->barrier();
};



FFTSolver::FFTSolver(Setup *setup, Parallel *_parallel, Geometry *_geo, double _Norm_XYZ, double _Norm_XY, double _Norm_X, double _Norm_Y) :
      parallel(_parallel), Norm_XYZ(_Norm_XYZ), Norm_XY(_Norm_XY), Norm_X(_Norm_X), Norm_Y(_Norm_Y), geo(_geo)
   {
     // is this really necessary ?
     //flags = FFT_X | FFT_Y;
     //if(setup->get("Vlasov.useAA" , 0) == 1)   flags |= FFT_AA;
       
     parseSuppressMode(setup->get("SuppressModeX", ""), suppressModeX);  
     parseSuppressMode(setup->get("SuppressModeY", ""), suppressModeY);  
   };


FFTSolver::~FFTSolver() {};


void FFTSolver::parseSuppressMode(const std::string &value, std::vector<int> &suppressMode) {
    if(value == "") return;


    std::string sub_str = Setup::eraseCharacter(value, " ") ;
    std::vector<std::string> result = Setup::split(sub_str, ",");
                     
    // needs cleanup !
    std::vector<std::string>::const_iterator it;
                                               for(it = result.begin(); it != result.end(); it++) {
                                                   if((*it).find("-") == std::string::npos) {
                                                   //if((*it).find("-") != -1) {
                                                    std::vector<std::string> res2 = Setup::split(*it, "-");
                                                    for(int i = atoi(res2[0].c_str()); i < atoi(res2[1].c_str()); i++) {
                                                        suppressMode.push_back(i);
                                                    }} else
                                                    suppressMode.push_back(std::atoi((*it).c_str()));
                                                }
              
              }

/*
// BUG : Input SuppressModeX = 0 does not work (Crashes)
void FFTSolver::suppressModes(CComplex kXOut[Nq][NzLD][NkyLD][FFTSolver::X_NkxL]) 
{

  // suppress x,y, z - Modes
  //vector<int>::const_iterator mode;

  for(auto mode_x = suppressModeX.begin(); (!suppressModeX.empty()) && (mode_x <= suppressModeX.end()) ; mode_x++) { 
   if((*mode_x < K1xLlD) || (*mode_x > K1xLuD)) continue;
    int x_k = *mode_x;
    kXOut[:][:][:][x_k] = 0.
  }
//  for(mode = suppressModeY.begin(); mode <= suppressModeY.end(); mode++) {
  for(mode = suppressModeY.begin(); (!suppressModeY.empty()) && (mode <= suppressModeY.end()) ; mode++) { 
       if((*mode < NkyLlD) || (*mode > NkyLuD)) continue;
        kXOut[:][:][y_k][:] = 0.
  }


  return;
}
*/

