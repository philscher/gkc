/*
 * =====================================================================================
 *
 *       Filename:  MHD.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/10/2011 06:24:21 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef __MHD_H
#define __MHD_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <complex>
#include <string>

#include<omp.h>

#include "DataOutput.h"
#include "Input.h"




//typedef std::complex<double> cmplxd;
typedef _Complex double cmplxd;
typedef cmplxd(*A1z)[];

class MHD {

  public:
    static const int Nx  = 2048;   ///< Define number of x points
    static const int Nky = 10  ;   ///< Define number of poloidal modes
     
  private:  
    double _kw_2_dx,               ///<  \f$ \frac{1}{2 dx} \f$
           _kw_dx2;                ///<  \f$ \frac{1}{dx^2} \f$

    bool doNonLinear;
    DataOutput *data;
    Input      *input;

    __declspec(align(64)) cmplxd  
              Vor[Nky][Nx], ///< Vorticity
              Cur[Nky][Nx], ///< Current
              Phi[Nky][Nx], ///< Kinetic
              Psi[Nky][Nx], ///< Magnetic
             dVor[Nky][Nx], ///< Derivative Magnetic
             dCur[Nky][Nx], ///< Derivative of Current
             dPhi[Nky][Nx], ///< Derivative of Flow
             dPsi[Nky][Nx], ///< Derivative of Flux
            ddPhi[Nky][Nx], ///< Derivative (2nd) of Flow
            ddPsi[Nky][Nx], ///< Derivative (2nd) of Flux
           VorOLD[Nky][Nx], ///<
           PsiOLD[Nky][Nx], ///< 
           PsiRHS[Nky][Nx], ///<
           VorRHS[Nky][Nx]; ///< 

    __declspec(align(64)) double 
                       X[Nx], ///< Domain in X
                    Psi0[Nx], ///< Equilibrium Flux
                   dPsi0[Nx], ///< Derivative of Equilibrium Flux
                    Cur0[Nx], ///< Current  (d2 psi_0 / dx^2)
                   dCur0[Nx]; ///< Derivative of Current 

    __declspec(align(64)) cmplxd 
                     ky[Nky]; ///< Poloidal wave number \f$ k_y = \frac{2\pi}{L_y}m \f$
   


    double Viscosity,         ///< \f$ \nu  \f$
           Resistivity;       ///< \f$ \eta \f$
  
  public:   

    MHD(std::string setup_filename, std::string setup_Xoptions);
    
    void startMainLoop();

    ~MHD();
  
  private:
  
	 void calculateRHSLinear(const int y_k);
	 void calculateRHSNonLinear();

    void initVariables(double Lx, double Ly);
    void calculateDerivatives(const int y_k);
    void calculateVor(const int y_k, const double dt) ;
    void swap_Old_New(const int y_k) ;
    void setBoundary(const int y_k);
   
    /**
	 *  using Tomas algorithm (source Wikipedia http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm)
    *  having fixed boundary conditions
    **/
    void solveTriDiagonalMatrix(cmplxd Sb[Nx], cmplxd D[Nx], cmplxd Sp[Nx], 
                                cmplxd  b[Nx], cmplxd X[Nky][Nx], const int ky);

};


#endif // __MHD_H

