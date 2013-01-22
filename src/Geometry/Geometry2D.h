/*
 * =====================================================================================
 *
 *       Filename: Geometry2D.h
 *
 *    Description: Two dimensional geometry (const theta, shear)
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#ifndef __GKC_GEOMETRY_2D_H__
#define __GKC_GEOMETRY_2D_H__

#include "Global.h"
#include "Geometry.h"

/**
*  @brief 2-d sheared slab geometry definition
*
*  ToDO : Parallelize the parallel boundary condition.
*
*
*  In the two-dimensional geometry, we set the \f$ k_\parallel \f$
*  is a function of, namely we distinguish between
*
*
*  \f$ k_\parallel = \theta k_y    \f$ - constant theta , local simulations with constant ange
*  \f$ k_\parallel = k_z           \f$ - sheareless slab, local simulations with constant ange
*  \f$ k_\parallel = \hat{s} x k_y \f$ - sheared slab   , non-local simulations with constant ange
*  
*  and other combinations (including extension to three dimensions).
*
*
*
*
**/
class Geometry2D : public Geometry
{

  double *By;
  std::string By_str;

  nct::allocate ArrayBy;

 public:

  double theta, ///< Angle to the magnetic field
         shear, ///< Magnetic field shear
         kz;    ///< Constant wavenumber factor

  Geometry2D(Setup *setup, Grid *grid, FileIO *fileIO) : Geometry(setup, grid, fileIO)  {
 
    // Parse magnetic field 
    theta               = setup->get("Geometry.Theta"   , 0.);
    shear               = setup->get("Geometry.Shear"   , 0.);
    By_str              = setup->get("Geometry.By"   , "shear*x+theta");
    kz                  = setup->get("Geometry.kz"   ,  0.);
  
    check(Lz == 1. ? 1 : -1, DMESG("In two dimensional simulations Lz should be 1"));

    FunctionParser By_parser = setup->getFParser();
    By_parser.AddConstant("shear", shear);
    By_parser.AddConstant("theta", theta);
    check(((By_parser.Parse(By_str, "x") == -1) ? 1 : -1), DMESG("Parsing error of shear condition"));

    ArrayBy = nct::allocate(nct::Range(NxGlB, NxGB))(&By);

    for(int x=NxGlB; x <= NxGuB; x++) By[x] = By_parser.Eval(&X[x]) ;

    Geometry::initData(fileIO);
    
    setupArrays();

  }

 
  /// \f$ J(\vec{x}) = 0 \f$
  inline double get_J(const int x, const int z) { return 1.; };

  /**  
  *    @name Defines the magnetic field and its variations
  *
  **/
  ///@{
  ///  \f$ g_{xx} = 1 \f$
  inline double g_xx(const int x, const int z) { return 1.0; };
  ///  \f$ g_{xy} = 0 \f$
  inline double g_xy(const int x, const int z) { return 0.0; }; 
  ///  \f$ g_{xz} = 0 \f$
  inline double g_xz(const int x, const int z) { return 0.0; };
  ///  \f$ g_{yy} = 1 \f$
  inline double g_yy(const int x, const int z) { return 1.0; };
  ///  \f$ g_{yz} = 0 \f$
  inline double g_yz(const int x, const int z) { return 0.0; };
  ///  \f$ g_{zz} = 0 \f$
  inline double g_zz(const int x, const int z) { return 0.0; };
  ///@}
  
  /**  
  *    @name Defines the magnetic field and its variations
  *
  **/
  ///@{
  /// \f$  B = 1 \f$
  inline double get_B      (const int x, const int z) { return 1.; };
  /// \f$  \partial_x B = 0 \f$ 
  inline double get_dB_dx  (const int x, const int z) { return 0.; };
  /// \f$  \partial_y B = 0 \f$
  inline double get_dB_dy  (const int x, const int z) { return 0.; };
  /// \f$  \partial_z B = 0 \f$
  inline double get_dB_dz  (const int x, const int z) { return 0.; };
  ///@}
 
  /**
  *
  *   Get the value of shear at position x, which is constant for sheared slab geometry 
  *
  **/
  inline CComplex get_kp(const int x, const CComplex ky, const int z) const  { 
    return ky * By[x] + ((CComplex (0. + 1.j)) * kz); 
  };
  
  double nu (const int x) { return 0.; };

  void printOn(std::ostream& output) const {
    output   << "Geometry   |  Sheared Slab   By : " << By_str << " Shear : " << shear << 
                                        " Theta : "  <<  theta << " kz : "    << kz    << std::endl;
  };

  void initData(hid_t geometryGroup) {
     
    check(H5LTset_attribute_string(geometryGroup, ".", "Type" , "2D"), DMESG("H5LTset_attribute"));

    check(H5LTset_attribute_double(geometryGroup, ".", "By"   ,  &By[NxGlD], Nx), DMESG("Setting attribute"));
    check(H5LTset_attribute_double(geometryGroup, ".", "Shear",  &shear    , 1 ), DMESG("Setting attribute"));
    check(H5LTset_attribute_double(geometryGroup, ".", "Theta",  &theta    , 1 ), DMESG("Setting attribute"));
    check(H5LTset_attribute_double(geometryGroup, ".", "kz"   ,  &kz       , 1 ), DMESG("Setting attribute"));
  };

};

#endif // __GKC_GEOMETRY_2D_H__

