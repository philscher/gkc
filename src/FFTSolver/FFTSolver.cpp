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
  solve(FFT_Type::Y_NL, FFT_SIGN::Forward, rY, kY);
  Norm_Y_Forward = creal(kY[0][0]);
  
  // Complex (Fourier sapce) -> Real transform
  kY[:][:] = 0. ; kY[0][0] = 1.;
  solve(FFT_Type::Y_NL, FFT_SIGN::Backward, kY, rY);
  Norm_Y_Backward = rY[0][0];

  /////////////////   Find normalization constants for X-transformation    //////////

  /* 
  // Real -> Complex (Fourier sapce) transform
  rXIn(RxLD, NkyLlD, NzLlD, 1) = 1.;
  solve(FFT_X_FIELDS, FFT_SIGN::Forward, rXIn.data());
  Norm_X_Forward = (K1xLlD == 0) ? real(kXOut(0, NkyLlD, NzLlD, 1)) : 0.;      


  // Complex (Fourier sapce) -> Real transform
  kXIn(Rk1xL, RkyLD, RzLD, 1) = 0.; if(K1xLlD == 0) kXIn(0, NkyLlD, NzLlD, 1) = 1.;
  solve(FFT_X_FIELDS, FFT_SIGN::Backward, kXOut.data());
  Norm_X_Backward = real(rXOut(NxLlD, NkyLlD, NzLlD, 1));
      
  // broadcast normalization to all nodes
  parallel->send(Norm_Y_Forward, parallel->Coord[DIR_ALL] == 0); parallel->send(Norm_Y_Backward, parallel->Coord[DIR_ALL] == 0);
  parallel->send(Norm_X_Forward, parallel->Coord[DIR_ALL] == 0); parallel->send(Norm_X_Backward, parallel->Coord[DIR_ALL] == 0);
  parallel->barrier();
   * */

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


/* 
// BUG : Input SuppressModeX = 0 does not work (Crashes)
void FFTSolver::suppressModes(CComplex kXOut[Nq][NzLD][NkyLD][FFTSolver::X_NkxL], const int field) 
{

  // suppress x,y, z - Modes
  vector<int>::const_iterator mode;

  for(mode = suppressModeX.begin(); (!suppressModeX.empty()) && (mode <= suppressModeX.end()) ; mode++) { 
   if((*mode < K1xLlD) || (*mode > K1xLuD)) continue;
    A(*mode, RkyLD, RzLD, RFields) = 1.e-128;
  }
//  for(mode = suppressModeY.begin(); mode <= suppressModeY.end(); mode++) {
  for(mode = suppressModeY.begin(); (!suppressModeY.empty()) && (mode <= suppressModeY.end()) ; mode++) { 
       if((*mode < NkyLlD) || (*mode > NkyLuD)) continue;
        A(Rk1xL, *mode, RzLD, RFields) = 1.e-128;
  }


  return;
}
 * */


