/*
 * =====================================================================================
 *
 *       Filename:  Interpolate.cpp
 *
 *    Description:  i
 *
 *        Version:  1.0
 *        Created:  07/15/2011 03:15:12 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */


#ifndef INTERPOLATION1D_H
#define INTERPOLATION1D_H


#include "Global.h"

// 2D interpolation

class Interpolate1D {
    Array1R A;
    Array1R Xa;

    std::string type;
  public:
    Interpolate1D(Array1R _X, Array1R _A, Setup *setup, std::string interpol = "Linear") : type(interpol) {



//         A.resize(Range(0,X.size),  Range(0, Y.size), Range(0, ZSize));


         A.reference(_A);
         Xa.reference(_X);


    };
  private:
    /**
     *
     *
     *
     *  Get position of nearest value
     */
    int getPosX(const double xv) { 
      int nx = -1;
      
      for(int x = 0; x < Xa.numElements(); x++) if((xv >= Xa(x)) && (xv <= Xa(x+1))) nx = x;
      if(nx == -1) check(-1, DMESG("No such point"));
      return nx;
    };



  public:
    double getValue(const double x, const double y, const double z=0.) {


          const int nx = getPosX[x];

    
          // Bi-Linear Interpolation
          
          const double  xL = Xa(nx);
          const double  xU = Xa(nx+1);


          return A(nx) + (x - xL) * (A(nx+1)-A(nx))/(xU-xL);
    };

};



#endif // INTERPOALTION1D_H
