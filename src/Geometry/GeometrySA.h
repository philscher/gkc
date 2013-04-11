/*
 * =====================================================================================
 *
 *       Filename: GeometrySA.h
 *
 *    Description: Implementation of the s-a Geometry
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#ifndef GEOMETRY_SA_H
#define GEOMETRY_SA_H

#include "Global.h"
#include "Geometry.h"


/**
*
*  @brief s-a geometry definitions
*
*  Gives the Geometric coefficients for a s-a geometry geometry with
*  a magnetic field according to.
*    
*  \f[ \vec{B}_0 / B_0 = \vec{b} = ???? \f].
*
*  Reference : Dannert, PhD, Chapter 2.2 Geometrie
*              Görler , PhD, Chapter 2.2 Geometry
*
*
*  \f$q_0, R_0, \hat{s} \f$
* 
*  The corresponding metric is 
*  
*  \f[
*
*    g^{ij} = 
*       \left( \begin{array}{lll}   
*  1                          &  \frac{\hat{s}} {q_0 R_0} z                                   & 0 \\
*  \frac{\hat{s}}{q_0}{R_0} z &  1 + \frac{z^2}{L_s^2} + \left( \frac{r_0}{q_0 R_0} \right)^2 & 0 \\
*                     0       &  \frac{q_0 R_0}{r_0}                                          & \left( \frac{q_0 R_0}{r_0} \right)^2 
*       
*   \end{array} \right)
*  \f]
*
*  r_0 is the minor radius, R_0 the major radius 
*
*  For the parallel boundary condition, we have to connect the poloidal modes
*  with the radial modes, thus we have to fulfill following quantization condition 
*    
*  \f[
*    \mathcal{N}{} = frac{2\pi\hat{s} L_x}{L_y}
*  \f]
*
**/
class GeometrySA : public Geometry
{
 
 public:

  double shear, ///< Magnetic shear
         eps,   ///< Aspect ratio    \f$ \epsilon = \frac{r_0}{R_0} f\$
         q0,    ///< Safety factor   
         R0,    ///< Major radius
         alpha; ///< Shafranov Shift

  double Bhat2;

  GeometrySA(Setup *setup, Grid *grid, FileIO *fileIO) : Geometry(setup, grid, fileIO) {
    
    Bhat2 = 1.;
    // set default values to Cyclon Base Case (CBC)
    shear = setup->get("Geometry.SA.Shear"         , 0.78 );
    q0    = setup->get("Geometry.SA.SafetyFactor"  , 1.40 );
    R0    = setup->get("Geometry.SA.MajorRadius"   , 1.0  );
    eps   = setup->get("Geometry.SA.AspectRatio"   , 0.19 );
    alpha = setup->get("Geometry.SA.ShafranovShift", 0.0  );
  
    setupArrays();
      
    Geometry::initData(fileIO);
   
  };

  /// \f$ J = \left( g_{zz} \left[ g_{xx} g_{yy} - g_{xy}^2 \right] \right)^{-\frac{1}{2}} \f$
  double get_J(const int x, const int z) {

    // where does this comes from
    return 1./sqrt( g_zz(x,z) * ( g_xx(x,z) * g_yy(x,z) - pow2(g_xy(x,z))));

  };

 private:
  
  /**  
  *  @name The metric coefficients
  *
  **/
  ///@{
  ///  \f$ g_{xx} = 1 \f$
  inline double g_xx(const int x, const int z) { return 1.0; };
  ///  \f$ g_{xy} = \frac{\hat{s}}{q_0 R_0} z \f$
  inline double g_xy(const int x, const int z) { return shear/(q0*R0)*Z[z]; }; 
  ///  \f$ g_{xx} = 0 \f$
  inline double g_xz(const int x, const int z) { return 0.0; };
  ///  \f$ g_{xx} = 1 + \left(\frac{\hat{s}}{q_0 R_0} \right)^2 z + \left( \epsilon q_0 \right)^2 \f$
  inline double g_yy(const int x, const int z) { return 1.  + pow2(shear/(q0*R0) * Z[z]) + pow2(eps/q0); };
  ///  \f$ g_{xx} = q_0 \epsilon \f$
  inline double g_yz(const int x, const int z) { return q0 / eps; };
  ///  \f$ g_{xx} = \left( q_0 \epsilon \right) \f$
  inline double g_zz(const int x, const int z) { return pow2(q0/eps); };
  ///@}
  
  /**  
  *    @name Defines the magnetic field and its variations
  *
  **/
  ///@{
  /// \f$  B = 1 \f$
  double get_B    (const int x, const int z) { return 1.;};
  /// \f$  \partial_x B = \frac{-\hat{B}}{R_0} \cos{(z)}  \f$
  double get_dB_dx(const int x , const int z) { return - Bhat2 / R0 * cos(Z[z]);            };
  /// \f$  \partial_y B = 0 \f$
  double get_dB_dy(const int x, const int z) { return 0.;                                  };
  /// \f$  \partial_x B = -\hat{B} \frac{\epsilon}{q_0 R_0} \cos{(z)}  \f$
  double get_dB_dz(const int x, const int z) { return - Bhat2 * eps/(q0 * R0) * cos(Z[z]); };
  ///@}
  
   
  /**
  *  \f[
  *    K_x = - 2 \frac{L_\perp}{R_0} \sin{(z)}           
  *  \f]
  *
  *  see Dannert PHD
  */
  virtual double get_Kx(const int x, const int z)  { 

    return - 2. * 1./R0 * sin(Z[z]);
  };
   
  /**
  *  \f[
  *    K_y = - 2 \frac{L_\perp}{R_0} \left( \cos{(z)} + \hat{s} z \sin{(z)} \right)           
  *  \f]
  *
  *  see Dannert PHD
  */
  virtual double get_Ky(const int x, const int z)  { 
    
    return -2. * 1./R0 * ( cos(Z[z]) + ( shear * Z[z]  - alpha * sin(Z[z])) * sin(Z[z]) );

  };


  /**
  * 
  *  \f[ 
  *    \nu = q_0 \left( 1 + \hat{s} \frac{x}{R_0} \right) 
  *  \f]
  *  
  *  as by first order Taylor approximation
  *  \f[
  *    q(x) \approx q_0 + x \partial_x q = q_0 \left( 1 + x/R_0 \hat{s} \right)
  *  \f]
  *
  */  
  inline double nu (const int x) { return  q0 * (1. + shear * (X[x]/R0)); };

 protected:

  virtual void printOn(std::ostream& output) const {

    output   << "Geometry   |  s-alpha Geometry Shear : " << shear << " Safety Factor   : " << q0 << std::endl;
    output   << "           |  Aspect Ratio           : " << eps   << " Shfaranov-Shift : " << alpha << std::endl;

  };
    
  void initData(hid_t GeometryGroup) {

    check(H5LTset_attribute_string(GeometryGroup, ".", "Geometry"      , "s-alpha"), DMESG("H5LTset_attribute"));
    check(H5LTset_attribute_double(GeometryGroup, ".", "Shear"         , &shear, 1), DMESG("H5LTset_attribute"));
    check(H5LTset_attribute_double(GeometryGroup, ".", "SafetyFactor"  , &q0   , 1), DMESG("H5LTset_attribute"));
    check(H5LTset_attribute_double(GeometryGroup, ".", "MajorRadius"   , &R0   , 1), DMESG("H5LTset_attribute"));
    check(H5LTset_attribute_double(GeometryGroup, ".", "AspectRatio"   , &eps  , 1), DMESG("H5LTset_attribute"));
    check(H5LTset_attribute_double(GeometryGroup, ".", "ShafranovShift", &alpha, 1), DMESG("H5LTset_attribute"));

  }

};

#endif // GEOMETRY_SA_H
