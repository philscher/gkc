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


#ifndef INTERPOLATION_H
#define INTERPOLATION_H


#include "Global.h"

// 2D interpolation

class Interpolate {
    Array2R A;
    Array1R Xa,Ya;

    std::string type;
  public:
    Interpolate(Array1R _X, Array1R _Y, Array2R _A, Setup *setup, std::string interpol = "Linear") : type(interpol) {



//         A.resize(Range(0,X.size),  Range(0, Y.size), Range(0, ZSize));


         A.reference(_A);
         Xa.reference(_X);
         Ya.reference(_Y);


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


    int getPosY(const double yv) {
      int ny = -1;
      for(int y = 0; y < Ya.numElements(); y++) if((yv >= Ya(y)) && (yv <= Ya(y+1))) ny = y;
      if(ny == -1) check(-1, DMESG("No such point"));
      return ny;
    };

  public:
    double getValue(const double x, const double y, const double z=0.) {


          const int nx = getPosX[x];
          const int ny = getPosY[y];

    
          // Bi-Linear Interpolation
          
          const double  xL = Xa(nx);
          const double  xU = Xa(nx+1);
          const double  yL = Ya(ny);
          const double  yU = Ya(ny+1);

          const double dyULdxUL = (xU - xL) * (yU - yL);

          return
             A(nx  , ny  )/dyULdxUL * (xU - x ) * (yU - y )
          +  A(nx+1, ny  )/dyULdxUL * (x  - xL) * (yU - y )
          +  A(nx  , ny+1)/dyULdxUL * (xU - x ) * (y  - yL)
          +  A(nx+1, ny+1)/dyULdxUL * (x  - xL) * (y  - yL);
          







    };

};



#endif // INTERPOALTION_H
