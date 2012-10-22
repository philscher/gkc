/*
 * =====================================================================================
 *
 *       Filename:  gkc-iCode.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  02/16/2012 10:23:13 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef GKC_ICODE
#define GKC_ICODE


#include <complex>


void setupMatrix(std::complex<double> w, double ky, double* X, double* K,  double Ls,  double Ln, int Nx, 
                std::complex<double>* A,  double h, double q, double mass, double T, double eta, double rho); 




#endif // GKC_ICODE
