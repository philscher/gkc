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


#ifndef GEOMETRY_2D_H
#define GEOMETRY_2D_H

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
  Array1d By;
  std::string shear_str;
  std::string By_str;


 public:
  double theta, shear, kz;




  Geometry2D(Setup *setup, FileIO *fileIO) : Geometry(setup, fileIO)  {
 
    // Parse magnetic field 
    theta               = setup->get("Geometry.Theta"   , 0.);
    shear               = setup->get("Geometry.Shear"   , 0.);
    By_str              = setup->get("Geometry.By"   , "shear*x+theta");
    kz                  = setup->get("Geometry.kz"   ,  0.);
   

   //std::cout << "SubString : " << setup->get("Geometry.By" , "0.").substr(0,4) << std::endl;
  // if (setup->get("Geometry.By" , "0.").substr(0,4) == "File") setFieldFromDataFile(setup, fields->Field0, Field::Ap, setup->get("Init.FixedAp", "0."));
  // else if (plasma->nfields >= 2) setFieldFromFunction(setup, fields->Field0, Field::Ap, setup->get("Init.FixedAp", "0."));
    
    FunctionParser By_parser = setup->getFParser();
    By_parser.AddConstant("shear", shear);
    By_parser.AddConstant("theta", theta);
    check(((By_parser.Parse(By_str, "x") == -1) ? 1 : -1), DMESG("Parsing error of shear condition"));

    By.resize(Range(NxGlD, NxGuD));
    for(int x=NxGlD; x <= NxGuD; x++) By(x) = By_parser.Eval(&X(x)) ;


        // BUG
    Geometry::initDataOutput(fileIO);
    
    setupArrays();

  };

  
  inline double get_J(const int x, const int z) { return 1.; };



   /**  
   *    @name Defines the magnetic field and its variations
   *
   **/
   ///@{
   ///  \f$ g_{xx} = 1 \f$
   inline double g_xx(const int x, const int z) { return 1.0;                 };
   ///  \f$ g_{xy} = 0 \f$
   inline double g_xy(const int x, const int z) { return 0.0;             }; 
   ///  \f$ g_{xz} = 0 \f$
   inline double g_xz(const int x, const int z) { return 0.0;                 };
   ///  \f$ g_{yy} = 1 \f$
   inline double g_yy(const int x, const int z) { return 1.0; };
   ///  \f$ g_{yz} = 0 \f$
   inline double g_yz(const int x, const int z) { return 0.0;                 };
   ///  \f$ g_{zz} = 0 \f$
   inline double g_zz(const int x, const int z) { return 0.0;                 };
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
   inline cmplxd get_kp(const int x, const cmplxd ky, const int z) const  { return ky * By(x) + cmplxd(0., kz); };
  
  
   double nu (const int x) { return 0.; };


   void printOn(ostream& output) const {
         output   << "Geometry   |  Sheared Slab   By : " << By_str  << " Shear : " << shear_str << " Theta : "  <<  theta << std::endl;
   };

   void initDataOutput(hid_t geometryGroup) {
          check(H5LTset_attribute_string(geometryGroup, ".", "Type", "2D"), DMESG("H5LTset_attribute"));
        //check(H5LTset_attribute_double(geometryGroup, ".", "By"   ,  By.data(), Nx), DMESG("H5LTset_attribute"));
          check(H5LTset_attribute_double(geometryGroup, ".", "Shear"   ,  &shear, 1), DMESG("H5LTset_attribute"));
          check(H5LTset_attribute_double(geometryGroup, ".", "Theta"   ,  &theta, 1), DMESG("H5LTset_attribute"));
          check(H5LTset_attribute_double(geometryGroup, ".", "kz"   ,  &kz, 1), DMESG("H5LTset_attribute"));
    }


//     virtual void writeData(Timing *timing) {};
//     virtual void closeData() {};

private:

  
};


#endif // GEOMETRY_2D_H
