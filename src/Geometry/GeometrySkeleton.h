/*
 * =====================================================================================
 *
 *       Filename: GeometrySkeleton.h
 *
 *    Description: Skeleton for Geometry
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#ifndef GEOMETRY_H
#define GEOMETRY_H

#include<iostream>

#include "Global.h"
#include "Setup.h"
#include "FFTSolver.h"




/** 
 *   class Geometry
 *
 *   Base class for the Geometry abstraction. For efficiency, we don't
 *   use virtual inheritance but CRTP. The latter enables to resolve all
 *   inline functions directly and thus should avoid any overhead due to 
 *   Geometry abstraction, see
     
 *
 *
 *  This part mostly relyies on Gorles, PhD, Eq. (2.51)  and goemetric parts
 *
 **/
struct ShearB {
  ShearB(const int _ly, const int _uy, const double _y0=0.) : ly(_ly),uy(_uy),  y0(_y0) {};
    const int ly, uy;
    const double   y0;
};



template<typename T> class Geometry : public IfaceGKC
{
public:
  const double LoCB
// C has remainder modulo operator, we need real one
int realmod(double x,double y) { return static_cast<int>(fmod((fmod(x,y) + y), y)); };
  /** The ratio between the parallel and the perpendicular scales */
  void printOn(ostream& o) { static_cast<T*>(this)->printOn (o); } ;

  bool isSymmetric;
  double eps_hat;

  Geometry(Setup *_setup, bool _isSymmetric=true) : isSymmetric(_isSymmetric), eps_hat(1.) , LoCB(1.){};
  virtual ~Geometry() {};

  // Jacobian
  inline const double J(const int x, const int y, const int z) { return static_cast<T*>(this)->J(x,y,z); };
   
  // ** Nabla 2 operator  // /
  inline const double k2_p(const int x_k, const int y_k, const int z) { return static_cast<T*>(this)->k2_p(x_k,y_k,z); };
  // BUG missing LB
  inline const double Kx(const int x, const int y, const int z) { return - LoCB * (dB_dy + g2/g1 * dB_dz); };
  inline const double Ky(const int x, const int y, const int z) { return   LoCB * (dB_dx + g3/g1 * dB_dz); };
  

  // is metric always symmetric or not ???!!! otherwise modify K3 term g_xy = g_yx 
  inline const double g1(const int x, const int y, const int z) { return g_xx(x,y,z) * g_yy(x,y,z) - pow2(g_xy(x,y,z))          ; };
  inline const double g2(const int x, const int y, const int z) { return g_xx(x,y,z) * g_yz(x,y,z) - g_xy(x,y,z) * g_xz(x,y,z)  ; };
  inline const double g3(const int x, const int y, const int z) { return g_xy(x,y,z) * g_yz(x,y,z) - g_yy(x,y,z) * g_xz(x,y,z)  ; };

  // define metric elements
  inline const double g_xx(const int x, const int y, const int z) { return static_cast<T*>(this)->g_xx(x,y,z); };
  inline const double g_xy(const int x, const int y, const int z) { return static_cast<T*>(this)->g_xy(x,y,z); };
  inline const double g_xz(const int x, const int y, const int z) { return static_cast<T*>(this)->g_xz(x,y,z); };
  inline const double g_yy(const int x, const int y, const int z) { return static_cast<T*>(this)->g_yy(x,y,z); };
  inline const double g_yz(const int x, const int y, const int z) { return static_cast<T*>(this)->g_yz(x,y,z); };
  inline const double g_zz(const int x, const int y, const int z) { return static_cast<T*>(this)->g_zz(x,y,z); };
 
  // define magnetic field and its variations
  inline const double B   (const int x, const int y, const int z) { return 1.;};
  
  inline const double dB_dx  (const int x, const int y, const int z) { return static_cast<T*>(this)->g_zz(x,y,z); };
  inline const double dB_dy  (const int x, const int y, const int z) { return static_cast<T*>(this)->g_zz(x,y,z); };
  inline const double dB_dz  (const int x, const int y, const int z) { return static_cast<T*>(this)->g_zz(x,y,z); };


  // Virtual Friend Function Idiom.
//   friend ostream& operator<<(ostream& output, Geometry const& g) {
//       g.printOn(output);  
//       return output;  // for multiple << operators.
///  };


  // for sheared magnetic fields, we have special boundary conditions
  inline const ShearB getYPos(const int x, const int y) { return static_cast<T*>(this)->getYPos(x,y); };

};



#endif // GEOMETRY_H
