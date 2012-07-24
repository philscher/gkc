/*
 * =====================================================================================
 *
 *       Filename: GeometrySA.h
 *
 *    Description: s-a Geomtry definitions
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#ifndef GEOMETRY_FIELDALIGNED_H
#define GEOMETRY_FIELDALIGNED_H

#include "Global.h"
#include "Setup.h"
#include "Geometry.h"
#include "FFTSolver.h"


/**
*   @brief s-ɑ geometry definitions
*
*    Gives the Geometric coefficients for a s-a geometry geometry with
*    a magnetic field accroding to.
*    
*    \f[ \vec{B}_0 / B_0 = \vec{b} = ???? \f].
*
*    Reference : Dannert, PhD , Chapter 2.2 Geometrie and
*                Goerler, PhD, Chaper 2.2 Geometry
*
*
*    \f$q_0, R_0, \hat{s} \f$
*
*    \f[
*          g^{ij} = 
*              \left( \begin{array}{lll}   
*     1                         &  \frac{\hat{s}} {q_0 R_0} z                                  & 0 \\
*     \frac{hat{s}}{q_0}{R_0} z &  1 + \frac{z^2}{L_s^2} + \left( \frac{r_0}{q_0 R_0} \right)^2 & 0 \\
*                       0       &  \frac{q_0 R_0}{r_0}                                          & \left( \frac{q_0 R_0}{r_0} \right)^2 
*                       \end{array} \right)
*    \f]
*
*    This results in the spatial operators of the form
*
*    \f[ \vec{b} \cdot \nabla = ? \f]
*    \f[ \nabla_\perp^2 = ? \f]
*    \f[ \vec{b} \times \nabla A \cdot \nabla = \f]
*
*  also, with
*
*
*  ToDO : Parallelize the parallel boundary condition.
*
*
**/
class GeometrySA : public Geometry
{
  public:

   double shear, ///< Magnetic shear
          eps,   ///< Aspect ratio 
          q0,    ///< Safety factor
          R0,    ///< Major radius
          alpha; ///< Shafranov Shift

    double Bhat2;

   GeometrySA(Setup *setup, FileIO *fileIO) : Geometry(setup, fileIO) {
    
   
      Bhat2 = 1.;
      // set default values to Cyclon Base Case (CBC)
      shear = setup->get("Geometry.SA.Shear"         , 0.16);
      q0    = setup->get("Geometry.SA.SafetyFactor"  , 1.5 );
      R0    = setup->get("Geometry.SA.MajorRadius"   , 1.5 );
      eps   = setup->get("Geometry.SA.AspectRatio"   , 0.19);
      alpha = setup->get("Geometry.SA.ShafranovShift", 0.0 );
  

      setupArrays();
    };

  
   /// \f$ J = \left( g_{zz} \left[ g_{xx} g_{yy} - g_{xy}^2 \right] \right)^{-\frac{1}{2}} \f$
   double get_J(const int x, const int z) {

      // where does this comes from
      return 1./sqrt( g_zz(x,z) * ( g_xx(x,z) * g_yy(x,z) - pow2(g_xy(x,z))));

    };

/**
    /// Nabla 2 operator /
    const double k2_p(const int x_k, const int y_k, const int z) {
      
//      const double kx = k(Nx,Lx,x_k);
//      const double ky = k(Ny,Ly,y_k);

//      return pow2(kx2) + g_xx(x,z) * pow2(ky) + 2.* g_xy(x,z) * kx * ky;

    };
*/

  private:
  
   /**  
   *    @name The metric coefficient the individual metric components
   *
   **/
   ///@{
   ///  \f$ g_{xx} = 1 \f$
   inline double g_xx(const int x, const int z) { return 1.0; };
   ///  \f$ g_{xy} = \frac{\hat{s}}{q_0 R_0} z \f$
   inline double g_xy(const int x, const int z) { return shear/(q0*R0)*Z(z); }; 
   ///  \f$ g_{xx} = 0 \f$
   inline double g_xz(const int x, const int z) { return 0.0; };
   ///  \f$ g_{xx} = 1 + \left(\frac{\hat{s}}{q_0 R_0} \right)^2 z + \left( \epsilon q_0 \right)^2 \f$
   inline double g_yy(const int x, const int z) { return 1.  + pow2(shear/(q0*R0) * Z(z)) + pow2(eps*q0); };
   ///  \f$ g_{xx} = q_0 \epsilon \f$
   inline double g_yz(const int x, const int z) { return q0 * eps; };
   ///  \f$ g_{xx} = \left( q_0 \epsilon \right) \f$
   inline double g_zz(const int x, const int z) { return pow2(q0*eps); };
   ///@}
  
   /**  
   *    @name Defines the magnetic field and its variations
   *
   **/
   ///@{
   /// \f$  B = 1 \f$
   double get_B      (const int x, const int z) { return 1.;};
   /// \f$  \partial_x B = \frac{-\hat{B}}{R_0} \cos{(z)}  \f$
   double get_dB_dx (const int x , const int z) { return - Bhat2 / R0 * cos(Z(z));            };
   /// \f$  \partial_y B = 0 \f$
   double get_dB_dy (const int x, const int z) { return 0.;                                  };
   /// \f$  \partial_x B = -\hat{B} \frac{\epsilon}{q_0 R_0} \cos{(z)}  \f$
   double get_dB_dz (const int x, const int z) { return - Bhat2 * eps/(q0 * R0) * cos(Z(z)); };
   ///@}
   
   inline double nu (const int x) { return 0.; };

  protected:

   virtual void printOn(ostream& output) const {
         output   << "Geometry  |  s-alpha Geometry Shear : " << shear << " Safety Factor   : " << q0 << std::endl;
         output   << "          |  Aspect Ratio           : " << eps   << " Shfaranov-Shift : " << alpha << std::endl;
   };
    
   void initDataOutput(hid_t GeometryGroup) {

          check(H5LTset_attribute_string(GeometryGroup, ".", "Geometry", "s-alpha"), DMESG("H5LTset_attribute"));
          check(H5LTset_attribute_double(GeometryGroup, ".", "Shear"         ,  &shear, 1), DMESG("H5LTset_attribute"));
          check(H5LTset_attribute_double(GeometryGroup, ".", "SafetyFactor"  ,  &q0   , 1), DMESG("H5LTset_attribute"));
          check(H5LTset_attribute_double(GeometryGroup, ".", "AspectRatio"   ,  &eps  , 1), DMESG("H5LTset_attribute"));
          check(H5LTset_attribute_double(GeometryGroup, ".", "ShafranovShift",  &alpha, 1), DMESG("H5LTset_attribute"));
    }


};


#endif // GEOMETRY_H
