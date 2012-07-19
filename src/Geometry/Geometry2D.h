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

#include "Geometry.h"
#include "Global.h"
#include "Setup.h"

#include "FileIO.h"
/**
*  @brief 2-d sheared slab geometry definition
*
*  ToDO : Parallelize the parallel boundary condition.
*
*
**/
class Geometry2D : public Geometry<Geometry2D>
{
  Array1d By;
  std::string shear_str;
  std::string By_str;
  double theta, shear, kz;


 public:

   void printOn(ostream& output) const {
         output   << "Geometry  |  Sheared Slab   By : " << By_str  << " Shear : " << shear_str << " Theta : "  <<  theta << std::endl;
   };



  Geometry2D(Setup *setup, FileIO *fileIO) : Geometry<Geometry2D>(setup, fileIO, true)  {
 
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
  };

  
  inline double J(const int x, const int y, const int z) { return 1.; };



   /**  
   *    @name Defines the magnetic field and its variations
   *
   **/
   ///@{
   ///  \f$ g_{xx} = 1 \f$
   inline double g_xx(const int x, const int y, const int z) { return 1.0;                 };
   ///  \f$ g_{xy} = 0 \f$
   inline double g_xy(const int x, const int y, const int z) { return 0.0;             }; 
   ///  \f$ g_{xz} = 0 \f$
   inline double g_xz(const int x, const int y, const int z) { return 0.0;                 };
   ///  \f$ g_{yy} = 1 \f$
   inline double g_yy(const int x, const int y, const int z) { return 1.0; };
   ///  \f$ g_{yz} = 0 \f$
   inline double g_yz(const int x, const int y, const int z) { return 0.0;                 };
   ///  \f$ g_{zz} = 0 \f$
   inline double g_zz(const int x, const int y, const int z) { return 0.0;                 };
   ///@}
  
   /**  
   *    @name Defines the magnetic field and its variations
   *
   **/
   ///@{
   /// \f$  B = 1 \f$
   inline double B      (const int x, const int y, const int z) { return 1.; };
   /// \f$  \partial_x B = 0 \f$
   inline double dB_dx  (const int x, const int y, const int z) { return 0.; };
   /// \f$  \partial_y B = 0 \f$
   inline double dB_dy  (const int x, const int y, const int z) { return 0.; };
   /// \f$  \partial_z B = 0 \f$
   inline double dB_dz  (const int x, const int y, const int z) { return 0.; };
   ///@}
 
   ///@{
   /// \f$  K_x = 0 \f$
   inline double Kx(const int x, const int y, const int z)  { return 0.; };
   /// \f$  K_y = 0 \f$
   inline double Ky(const int x, const int y, const int z)  { return 0.; };
   ///@}
  
   /**
   *
   *   Get the value of shear at position x, which is constant for sheared slab geometry 
   *
   **/
   inline cmplxd get_kp(const int x, const cmplxd ky, const int z)  { return ky * By(x) + cmplxd(0., kz); };
  


   // for sheared magnetic fields, we have special boundary conditions, because
   // we need to take care of possible shear
   inline ShearB getYPos(const int x,const int y) {

            return ShearB(y, y, 0.0  );//(yv - Y(yi-1))/dy);
              
   };


   std::string getGeometryName() { return "Sheared Slab 2D"; };

   void initDataOutput(FileIO *fileIO, hid_t geometryGroup) {
          check(H5LTset_attribute_string(geometryGroup, ".", "Type", "2D"), DMESG("H5LTset_attribute"));
        //check(H5LTset_attribute_double(geometryGroup, ".", "By"   ,  By.data(), Nx), DMESG("H5LTset_attribute"));
          check(H5LTset_attribute_double(geometryGroup, ".", "Shear"   ,  &shear, 1), DMESG("H5LTset_attribute"));
          check(H5LTset_attribute_double(geometryGroup, ".", "Theta"   ,  &theta, 1), DMESG("H5LTset_attribute"));
          check(H5LTset_attribute_double(geometryGroup, ".", "kz"   ,  &kz, 1), DMESG("H5LTset_attribute"));
    }


// transformation
  double getY(const int x, const int y, const int z) {
	return Y(y);
 }
   
     virtual void writeData(Timing *timing) {};
     virtual void closeData() {};

private:

  
};


#endif // GEOMETRY_2D_H
