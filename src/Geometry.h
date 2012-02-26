/*
 * =====================================================================================
 *
 *       Filename:Geometry.h
 *
 *    Description: Geometry interface
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#ifndef GEOMETRY_H
#define GEOMETRY_H


#include "Global.h"
#include "Setup.h"

#include "FileIO.h"


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
 * BUG : if u compile and derive and doe not implement one of the function you
 *        will get a self-recursion. Any idea how to fix this  ? 
 *        See
 *        http://stackoverflow.com/questions/4721047/this-is-crtp-usage-for-static-polymorphism-but-without-implementation-of-a-derive
 *
 **/
struct ShearB {
  ShearB(const int _ly, const int _uy, const double _y0=0.) : ly(_ly),uy(_uy),  y0(_y0) {};
    const int ly, uy;
    const double   y0;
};



template<typename T> class Geometry : public IfaceHelios
{

public:
  // C has remainder modulo operator, we need real one
  int realmod(double x,double y) { return static_cast<int>(fmod((fmod(x,y) + y), y)); };
  /** The ratio between the parallel and the perpendicular scales */
//  void printOn(ostream& o) { static_cast<T*>(this)->printOn (o); } ;

  bool isSymmetric;
  double eps_hat;
  const double LoCB;

  Geometry(Setup *setup, FileIO *fileIO, bool _isSymmetric=true) : isSymmetric(_isSymmetric) , LoCB(1.){

        eps_hat = setup->get("Geometry.eps_hat", 1.);

  
  };

  virtual ~Geometry() {};

  // Jacobian
  inline  double J(const int x, const int y, const int z)  { return static_cast<T*>(this)->J(x,y,z); };
   
  // BUG missing LB // we should pre-caluclate this value // write defnition of these and reference
  // Definition
  inline  double Kx(const int x, const int y, const int z)  { return static_cast<T*>(this)->Kx(x,y,z); };
  inline  double Ky(const int x, const int y, const int z)  { return static_cast<T*>(this)->Ky(x,y,z); };
  
  
  /**
   *
   *   Get the value of shear at position x 
   *
   */
  inline  double getShear(const int x)  { return static_cast<T*>(this)->getShear(x); };
  

  // is metric always symmetric or not ???!!! otherwise modify K3 term g_xy = g_yx
  //Define metric helper gterms gamma
  
  /**  Helper function
   *
   *   \f[ \gamma_1 = g_{xx} g_{yy} - g_{xy} g_{yx} \f]
   */ 
  inline  double g1(const int x, const int y, const int z) { return g_xx(x,y,z) * g_yy(x,y,z) - pow2(g_xy(x,y,z))          ; };
  /**  Helper function
   *
   *   \f[ \gamma_2 = g_{xx} g_{yz} - g_{xy} g_{xz} \f]
   */ 
  inline  double g2(const int x, const int y, const int z) { return g_xx(x,y,z) * g_yz(x,y,z) - g_xy(x,y,z) * g_xz(x,y,z)  ; };
  /**  Helper function
   *
   *   \f[ \gamma_3 = g_{xy} g_{yz} - g_{yy} g_{xz} \f]
   */ 
  inline  double g_3(const int x, const int y, const int z) { return g_xy(x,y,z) * g_yz(x,y,z) - g_yy(x,y,z) * g_xz(x,y,z)  ; };

  // define metric elements
  inline  double g_xx(const int x, const int y, const int z) { return static_cast<T*>(this)->g_xx(x,y,z); };
  inline  double g_xy(const int x, const int y, const int z) { return static_cast<T*>(this)->g_xy(x,y,z); };
  inline  double g_xz(const int x, const int y, const int z) { return static_cast<T*>(this)->g_xz(x,y,z); };
  inline  double g_yy(const int x, const int y, const int z) { return static_cast<T*>(this)->g_yy(x,y,z); };
  inline  double g_yz(const int x, const int y, const int z) { return static_cast<T*>(this)->g_yz(x,y,z); };
  inline  double g_zz(const int x, const int y, const int z) { return static_cast<T*>(this)->g_zz(x,y,z); };
  
  // minor simplification of k2_p
  inline  double g_xx(const int z) { return static_cast<T*>(this)->g_xx(-1,-1,z); };
  inline  double g_xy(const int z) { return static_cast<T*>(this)->g_xy(-1,-1,z); };
  inline  double g_yy(const int z) { return static_cast<T*>(this)->g_yy(-1,-1,z); };
 
  // define magnetic field and its variations
  
  /**  Defined Magnetic field strength \f[ B_0(\vec{x}) \f] at position \f[ \vec{x} = (x,y,z) \f]
   *
   */ 
  inline  double B   (const int x, const int y, const int z)    { return static_cast<T*>(this)->B(x,y,z); };
  
  /**  Defined Magnetic field strength \f[ \frac{\partial B_0}{\partial x}(\vec{x}) \f] at position \f[ \vec{x} = (x,y,z) \f]
   *
   */ 
  inline  double dB_dx  (const int x, const int y, const int z) { return static_cast<T*>(this)->dB_dx(x,y,z); };

  /**  Defined Magnetic field strength \f[\frac{\partial B_0}{\partial y}(\vec{x}) \f] at position \f[ \vec{x} = (x,y,z) \f]
   *
   */ 
  inline  double dB_dy  (const int x, const int y, const int z) { return static_cast<T*>(this)->dB_dy(x,y,z); };
  /**  Defined Magnetic field strength \f[\frac{\partial B_0}{\partial z}(\vec{x}) \f] at position \f[ \vec{x} = (x,y,z) \f]
   *
   */ 
  inline  double dB_dz  (const int x, const int y, const int z) { return static_cast<T*>(this)->dB_dz(x,y,z); };


  // Virtual Friend Function Idiom.
//   friend ostream& operator<<(ostream& output, Geometry const& g) {
//       g.printOn(output);  
//       return output;  // for multiple << operators.
///  };


  // for sheared magnetic fields, we have special boundary conditions
  inline  ShearB getYPos(const int x, const int y) { return static_cast<T*>(this)->getYPos(x,y); };
  
 
  inline  double getY (const int x, const int y, const int z) { return static_cast<T*>(this)->getY(x,y,z); };

  // std::string getGeometryName() { return static_cast<T*>(this)->getGeometryName(); };

   // dataOutput Stuff
  // void printOn(ostream& o) { static_cast<T*>(this)->printOn (o); } ;
   
   void initDataOutput(FileIO *fileIO) {
        hid_t geometryGroup = check(H5Gcreate(fileIO->getFileID(), "/Geometry",H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), DMESG("Error creating group for Geometry : H5Gcreate"));
  //    check(H5LTset_attribute_string(geometryGroup, ".", "Type", static_cast<T*>(this)->getGeometryName()), DMESG("H5LTset_attribute"));
        static_cast<T*>(this)->initDataOutput(fileIO, geometryGroup); 
        H5Gclose(geometryGroup);
   } ;
        

//     };
     virtual void writeData(Timing *timing) {};
     virtual void closeData() {};


};


        //static_cast<T*>(this)->initDataOutput(fileIO, geometryGroup); 

#endif // GEOMETRY_H
