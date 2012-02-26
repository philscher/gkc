/*
 * =====================================================================================
 *
 *       Filename: GeometrySlab.h
 *
 *    Description: Definition of shearless slab geomtry
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#ifndef GEOMETRY_SLAB_H
#define GEOMETRY_SLAB_H

#include<iostream>

#include "Global.h"
#include "Setup.h"
#include "FFTSolver.h"
#include "Geometry.h"


#include "FileIO.h"


/** GeometrySlab 
 *
 *
 * Basic Geometry class for basic shearless slab geometry. In this case
 * the magnetic field is written as
 *
 *  \f[ \vec{B}_0 / B_0 = \vec{b} = \left(0,0,1 \right) \f]
 *  
 *  thus the spatial terms becomes ....
 *
 *   The perpendicular gradient squared is 
 *   \f[ k_\perp^2 = k_x^2 + k_y^2 \f]
 *
 * */
class GeometrySlab : public Geometry<GeometrySlab>
{
public:

   void printOn(ostream& o)  const {
         o  << "Geometry   |  Shearless Slab" << std::endl;
   };


  GeometrySlab(Setup *_setup, FileIO *fileIO) : Geometry<GeometrySlab>(_setup, fileIO, true) { };
  inline  double J(const int x, const int y, const int z) const { return 1.;};
   
  // define metric elements
  inline  double g_xx(const int x, const int y, const int z) { return 1.0; };
  inline  double g_xy(const int x, const int y, const int z) { return 0.0; }; 
  inline  double g_xz(const int x, const int y, const int z) { return 0.0; };
  inline  double g_yy(const int x, const int y, const int z) { return 1.0; };
  inline  double g_yz(const int x, const int y, const int z) { return 0.0; };
  inline  double g_zz(const int x, const int y, const int z) { return 1.0; };
  
  // define magnetic field and its variations
  inline  double B      (const int x, const int y, const int z) { return 1.;};
  
  inline  double dB_dx  (const int x, const int y, const int z) { return 0.;};
  inline  double dB_dy  (const int x, const int y, const int z) { return 0.;};
  inline  double dB_dz  (const int x, const int y, const int z) { return 0.;};
 
  inline  double Kx(const int x, const int y, const int z) const { return 0.; };
  inline  double Ky(const int x, const int y, const int z) const { return 0.; };

  // no shearing, so simply y
  inline ShearB getYPos(const int x, const int y) { return ShearB(y, y, 0.); };

 std::string getGeometryName() { return "Shearless Slab"; };

    void initDataOutput(FileIO *fileIO, hid_t geometryGroup) {
          check(H5LTset_attribute_string(geometryGroup, ".", "Type", "Shearless Slab"), DMESG("H5LTset_attribute"));
          check(H5LTset_attribute_double(geometryGroup, ".", "eps_hat"   ,  &eps_hat, 1), DMESG("H5LTset_attribute"));
    }


// transformation
  double getY(const int x, const int y, const int z) {
	return Y(y);
 }
//   void initDataOutput(FileIO *fileIO) {
 //  } ;
     virtual void initDataOutput(FileIO *fileIO) {};
     virtual void writeData(Timing *timing) {};
     virtual void closeData() {};
};


#endif // GEOMETRY_SLAB_H
