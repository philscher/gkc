/*
 * =====================================================================================
 *
 *       Filename: GeometryCHEASE.h
 *
 *    Description: Geometry using Chease numerical MHD equilibrium
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#ifndef GEOMETRY_CHEASE_H
#define GEOMETRY_CHEASE_H

#include "Geometry.h"
#include "Global.h"

/**
*  @brief Geometry defintions using numerical MHD equilibrium
*         from CHEASE code 
*
*  @warning not working
*  @todo  add link to Chease
*
**/
class GeometryCHEASE : public Geometry
{


  std::string chease_file;
 public:



  GeometryCHEASE(Setup *setup, FileIO *fileIO) : Geometry(setup, fileIO) 
  {

    chease_file = setup->get("Geometry.MHDFile", "");

  }

   /// \f$ J = 1 \f$
  inline double get_J(const int x, const int z) { return 1.; };

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
   **/
   ///@{
   /// \f$  B = 1 \f$
   inline  double B      (const int x, const int z) { return 1.;};
   /// \f$  \partial_x B = 0 \f$
   inline  double dB_dx  (const int x, const int z) { return 0.;};
   /// \f$  \partial_y B = 0 \f$
   inline  double dB_dy  (const int x, const int z) { return 0.;};
   /// \f$  \partial_z B = 0 \f$
   inline  double dB_dz  (const int x, const int z) { return 0.;};
   ///@}

   // how to connect the field lines ?
   double nu (const int x) { return 0.; };


private:

  // check out GS2 -> license GPLv2+ ?
  // check out GKW -> license GPLv2+ ?
   void loadCheaseFile() {};

protected:

   void printOn(std::ostream& output) const {
         output   << "Geometry  |  CHEASE (Numerical Equilibrium) Data Input : " << "mysterious file" << std::endl;
   };


   void initDataOutput(hid_t geometryGroup) 
   {
          check(H5LTset_attribute_string(geometryGroup, ".", "Type", "CHEASE"), DMESG("H5LTset_attribute"));
   }


   //virtual void writeData(Timing *timing) {};
   //virtual void closeData() {};

};


#endif // GEOMETRY_CHEASE_H


