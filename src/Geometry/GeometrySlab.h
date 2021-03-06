/*
 * =====================================================================================
 *
 *       Filename: GeometrySlab.h
 *
 *    Description: Definition of three dimensional shearless slab geomtry
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#ifndef GEOMETRY_SLAB_H
#define GEOMETRY_SLAB_H

#include "Global.h"
#include "Geometry.h"


/** 
*   @brief Three dimensional shearless slab geometry
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
*   \f$ k_{1\parallel} \f$ cannot be fixed but is defined through the
*   parallel length and number of discretization points, such that
*   \f$ k_{1\parallel} = \frac{2\pi}{L_z} n \f$.
*
**/
class GeometrySlab : public Geometry
{
public:

  GeometrySlab(Setup *setup, Grid *grid, FileIO *fileIO) : Geometry(setup, grid, fileIO) 
  { 
  
  
  };


  /// \f$ J = 1 \f$
  double get_J(const int x, const int z)  { return 1.;};
   
  /**  
  *    @name The metric coefficient the individual metric components
  **/
  ///@{
  ///  \f$ g_{xx} = 1 \f$
  inline  double g_xx(const int x, const int z) { return 1.0; };
  ///  \f$ g_{xy} = 0 \f$
  inline  double g_xy(const int x, const int z) { return 0.0; }; 
  ///  \f$ g_{xz} = 0 \f$
  inline  double g_xz(const int x, const int z) { return 0.0; };
  ///  \f$ g_{yy} = 1 \f$
  inline  double g_yy(const int x, const int z) { return 1.0; };
  ///  \f$ g_{yz} = 0 \f$
  inline  double g_yz(const int x, const int z) { return 0.0; };
  ///  \f$ g_{zz} = 1 \f$
  inline  double g_zz(const int x, const int z) { return 1.0; };
  ///@}


  /**  
  *    @name Defines the magnetic field and its variations
  *
  *    Constant magnetic field through the whole domain,
  *    thus gradients are zero.
  *
  **/
  ///@{
  /// \f$  B = 1 \f$
  inline  double get_B      (const int x, const int z) { return 1.;};
  /// \f$  \partial_x B = 0 \f$
  inline  double get_dB_dx  (const int x, const int z) { return 0.;};
  /// \f$  \partial_y B = 0 \f$
  inline  double get_dB_dy  (const int x, const int z) { return 0.;};
  /// \f$  \partial_z B = 0 \f$
  inline  double get_dB_dz  (const int x, const int z) { return 0.;};
  ///@}

  // how to connect the field lines ?
  double nu (const int x) { return 0.; };


  void printOn(std::ostream& output)  const {
         output  << "Geometry   |  Shearless Slab" << std::endl;
  };


  void initData(hid_t geometryGroup) 
  {
         check(H5LTset_attribute_string(geometryGroup, ".", "Type", "Shearless Slab"), DMESG("H5LTset_attribute"));
  }



  //virtual void writeData(Timing *timing) {};
  //virtual void closeData() {};

};


#endif // GEOMETRY_SLAB_H
