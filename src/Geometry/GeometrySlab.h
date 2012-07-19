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


/** 
*   @brief Shearless Slab geometry definitions
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
**/
class GeometrySlab : public Geometry<GeometrySlab>
{
public:

   void printOn(ostream& o)  const {
         o  << "Geometry   |  Shearless Slab" << std::endl;
   };


   GeometrySlab(Setup *_setup, FileIO *fileIO) : Geometry<GeometrySlab>(_setup, fileIO, true) 
   { 
  
  
   };


   /// \f$ J = 1 \f$
   inline  double J(const int x, const int y, const int z) const { return 1.;};
   
   /**  
   *    @name The metric coefficient the individual metric components
   **/
   ///@{
   ///  \f$ g_{xx} = 1 \f$
   inline  double g_xx(const int x, const int y, const int z) { return 1.0; };
   ///  \f$ g_{xy} = 0 \f$
   inline  double g_xy(const int x, const int y, const int z) { return 0.0; }; 
   ///  \f$ g_{xz} = 0 \f$
   inline  double g_xz(const int x, const int y, const int z) { return 0.0; };
   ///  \f$ g_{yy} = 1 \f$
   inline  double g_yy(const int x, const int y, const int z) { return 1.0; };
   ///  \f$ g_{yz} = 0 \f$
   inline  double g_yz(const int x, const int y, const int z) { return 0.0; };
   ///  \f$ g_{zz} = 1 \f$
   inline  double g_zz(const int x, const int y, const int z) { return 1.0; };
   ///@}
 

   /**  
   *    @name Defines the magnetic field and its variations
   **/
   ///@{
   /// \f$  B = 1 \f$
   inline  double B      (const int x, const int y, const int z) { return 1.;};
   /// \f$  \partial_x B = 0 \f$
   inline  double dB_dx  (const int x, const int y, const int z) { return 0.;};
   /// \f$  \partial_y B = 0 \f$
   inline  double dB_dy  (const int x, const int y, const int z) { return 0.;};
   /// \f$  \partial_z B = 0 \f$
   inline  double dB_dz  (const int x, const int y, const int z) { return 0.;};
   ///@}
 
   /**  
   *    @name The metric coefficient the individual metric components
   **/
   ///@{
   /// \f$ K_x = 0 \f$
   inline  double Kx(const int x, const int y, const int z) const { return 0.; };
   /// \f$ K_y = 0 \f$
   inline  double Ky(const int x, const int y, const int z) const { return 0.; };
   ///@}

   // no shearing, so simply y
   inline ShearB getYPos(const int x, const int y) { return ShearB(y, y, 0.); };

   std::string getGeometryName() { return "Shearless Slab"; };

   void initDataOutput(FileIO *fileIO, hid_t geometryGroup) 
   {
          check(H5LTset_attribute_string(geometryGroup, ".", "Type", "Shearless Slab"), DMESG("H5LTset_attribute"));
          check(H5LTset_attribute_double(geometryGroup, ".", "eps_hat"   ,  &eps_hat, 1), DMESG("H5LTset_attribute"));
   }


   // transformation
   double getY(const int x, const int y, const int z) {
	   return Y(y);
   }

   //   void initDataOutput(FileIO *fileIO) {
   virtual void initDataOutput(FileIO *fileIO) {};
   virtual void writeData(Timing *timing) {};
   virtual void closeData() {};

};


#endif // GEOMETRY_SLAB_H
