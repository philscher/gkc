/*
 * =====================================================================================
 *
 *       Filename: Auxiliary.cpp
 *
 *    Description: Auxiliary diagnostic functions for gkc
 *
 *         Author: Paul P. Hilscher (2013), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef _GKC__AUXILIARY_H__
#define _GKC__AUXILIARY_H__

#include "Global.h"
#include "Setup.h"

#include "Grid.h"
#include "Vlasov/Vlasov.h"
#include "Fields/Fields.h"

/**
*    @brief Auxiliary functions
*
*    Currently includes calculation of zonal flow production rates
*
*
*
*
**/
class Auxiliary
{
  Vlasov *vlasov;
  Fields *fields;
  FFTSolver *fft;
  Parallel *parallel;

  FileAttr *FA_ZFProd  ,
           *FA_ZFProdTime;
   
  CComplex *ZF_Gyro_In, *ZF_Gyro_Out;
   
  nct::allocate ArrayZF,    // Zonal flow production array
                ArrayZFGy;
 
  CComplex *ZFProd;
  bool zonalFlow; ///< calculate Zonal flow production rate
  
  Timing dataOutputZF; ///< Timing to output whole phase distribution function
   
  hid_t auxGroup;


 public:

  Auxiliary(Setup *setup, FileIO *fileIO, Parallel *parallel, Fields *fields, Vlasov *vlasov, Grid *grid, FFTSolver *fft);
 ~Auxiliary();

  void writeData(const Timing &timing, const double dt); 
   
  void printOn(std::ostream &output) const;
  
 private:

  /**
  *    Calculates the zonal flow production rate through the 
  *    non-linear terms.
  *
  *
  **/
  void calculateZonalFlowProduction(const double dt);
};

#endif // _GKC__AUXILIARY_H__
