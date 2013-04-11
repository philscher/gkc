/*
 * =====================================================================================
 *
 *       Filename: GeometryShear.h
 *
 *    Description: Definition of three dimensional sheared slab geometry
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#ifndef GEOMETRY_SHEAR_H
#define GEOMETRY_SHEAR_H

#include "Global.h"
#include "Geometry.h"


/**
*
* @brief Sheared Slab geometry definitions
*
* Gives the Geometric coefficients for a sheared slab geometry with
* a magnetic field according to.
*
* \f[ 
*   \vec{B}_0 / B_0 = \vec{b} = \left(0,-x/L_s,1 \right) 
* \f].
*
*  Reference : @cite Jenko_2001:PhDThesis , Section 3.5 Flussschlauchgeometry
*
*  where $L_s$ is defined as the connection length. The metric takes
*  the simple form 
*
*  \f[
*     
*     g^{ij} = 
*     
*       \left( \begin{array}{lll}
*                1     &  z/L_s                 & 0 \\
*                z\L_s &  1 + \frac{z^2}{L_s^2} & 0 \\
*                0     &  0                     & 1
*              \end{array} \right)
*  \f]
*
*  This results in the spatial operators of the form
*
*  \f[ 
*           \vec{b} \cdot \nabla = \partial_z 
*  \f]
*  and
*  \f[ 
*           \nabla_\perp^2 = \frac{\partial^2}{\partial x^2} + \left(1 + \frac{z^2}{L_s^2} \right)
*              \frac{\partial^2}{\partial y^2} + 2 \frac{z}{L_s} \frac{\partial^2}{\partial x \partial y} 
*  \f]
*
*  \f[
*          \vec{b} \times \nabla A \cdot \nabla = \frac{\partial A}{\partial x}{\partial}{\partial y} -
*                  \frac{\partial A}{\partial y}{\partial }{\partial x}  
*  \f]
*
*  also, with
*
*  \f[
*    L_s = 
*    L_c = 2 \pi q R
*    L_s =  \frac{q R}{\hat{s}}
*  \f]
*  with \f$ \hat{s} \f$ the shear. For the parallel length L_z, we need to choose the connection length 
*  L_c. So once L_s and L_z is given we calculate the shear to
*  \f[ \hat{s} = \frac{L_z }{2 \pi L_s} \f]
*  
*  for consistency with the periodic parallel boundary conditions we need to fulfill the 
*  relation
*  \f[ 
*    2 \pi \hat{s} L_x = n_s L_y  or \frac{L_z}{L_s} = n * L_y 
*  \f]
*  where n_s is an integer number, or an interpolation method needs to be used. 
*
*  ToDO : Parallelize the parallel boundary condition.
*
**/
class GeometryShear : public Geometry
{
  double Ls,    ///< Shearing length  
         shear; ///< Shearing

  bool connectFieldLines,   ///< True if field lines are connected at z-boundary
       roundShearToConnect;

 public:

  GeometryShear(Setup *setup, Grid *grid, FileIO *fileIO) : Geometry(setup, grid, fileIO)  {
     

    shear               = setup->get("Geometry.Shear"   , 0.);
    connectFieldLines   = setup->get("Geometry.ConnectFieldLines", 1);
    roundShearToConnect = setup->get("Geometry.RoundShearToConnect", 1);

    //check connection length, otherwise parallel boundary will fail
    // modf :  extract signed integral and fractional values from floating-point number
    if((abs(fmod(2. * M_PI * shear * Lx/Ly, 1.)) > 1.e-5) && roundShearToConnect) {

      std::cout << "Warning : Rounding shear to ensure magnetic field line connection" << std::endl;
      const double di = round(2. * M_PI * shear * Lx/Ly);
      shear = di / (2. * M_PI) * Ly/Lx;
    }

    // if shear is zero, set connection length to high number (what about FPE)
    Ls = (shear != 0.) ?   Lz / (2. * M_PI * shear) : 1.e99;
   
  };

 
  /// \f$ J = 1 \f$
  inline double get_J(const int x, const int z) { return 1.; };

  /**  
  *    @name The metric coefficient the individual metric components
  **/
  ///@{
  ///  \f$ g_{xx} = 1 \f$
  inline double g_xx(const int x, const int z) { return 1.0;                 };
  ///  \f$ g_{xy} = z / L_s \f$
  inline double g_xy(const int x, const int z) { return Z[z]/Ls;             }; 
  ///  \f$ g_{xz} = 0 \f$
  inline double g_xz(const int x, const int z) { return 0.0;                 };
  ///  \f$ g_{yy} = 1 + \left( \frac{z}{L_s} \right)^2 \f$
  inline double g_yy(const int x, const int z) { return 1.0 + pow2(Z[z]/Ls); };
  ///  \f$ g_{yz} = 0 \f$
  inline double g_yz(const int x, const int z) { return 0.0;                 };
  ///  \f$ g_{zz} = 1 \f$
  inline double g_zz(const int x, const int z) { return 1.0;                 };
  ///@}
  
  /**  
  *    @name Defines the magnetic field and its variations
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
 
 
  // how to connect the field lines ?
  double nu (const int x) { return -shear * X[x]; };
   


  void printOn(std::ostream& output) const {
         
     output   << "Geometry  |  Sheared Slab   shear : " << shear << " (Ls : " << Ls << " ) "  << " eps_hat : " << eps_hat << std::endl;
     output   << "          |  Connect Field Lines (" << (connectFieldLines ? "true" : "false") << 
                                   ")   RoundShear (" << (roundShearToConnect ? "true" : "false") << ")" << std::endl;
  };

  void initData(hid_t geometryGroup) {

    check(H5LTset_attribute_string(geometryGroup, ".", "Type", "Sheared Slab"), DMESG("H5LTset_attribute"));
    check(H5LTset_attribute_string(geometryGroup, ".", "ConnectFieldLines", connectFieldLines ? "true" : "false"), DMESG("H5LTset_attribute"));
    check(H5LTset_attribute_string(geometryGroup, ".", "RoundShearToConnect", roundShearToConnect ? "true" : "false"), DMESG("H5LTset_attribute"));
      
    check(H5LTset_attribute_double(geometryGroup, ".", "Shear"   ,  &shear, 1), DMESG("H5LTset_attribute"));
    check(H5LTset_attribute_double(geometryGroup, ".", "eps_hat"   ,  &eps_hat, 1), DMESG("H5LTset_attribute"));

  };

};

#endif // GEOMETRY_SHEAR_H
