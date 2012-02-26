/*
 * =====================================================================================
 *
 *       Filename: GeometryCircular.h
 *
 *    Description: Circular geometry definition
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
 *
 *    Gives the Geometric coefficients for a s-a geometry geometry with
 *    a magnetic field accroding to.
 *    \f[ \vec{B}_0 / B_0 = \vec{b} = ???? \f].
 *
 *    Reference : Dannert, PhD , Chapter 2.2 Geometrie \and
 *                Goerler, PhD, Chaper 2.2 Geometry
 *
 *
 *    \f $q_0, R_0, \hat{s}$ \f
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
 * */
class GeometrySA : public Geometry<GeometrySA>
{
 public:

  /** Magnetic shear */
  double shear, eps, q;
  
  /** safety factor */
  /** Curvature term 1 */
  
  
  inline const double Kx(const int x, const int y, const int z) const { return - 2. * 2. * sin(z); };
  inline const double Ky(const int x, const int y, const int z) const { return - (cos(Z(z)) + shear * Z(z) * sin(Z(z))); };

  GeometrySA(Setup *_setup) : Geometry<GeometrySA>(setup, true, true) {
    
    
    // set default values to Cyclon Base Case (CBC)
    shear = setup->get("GeometrySA.Shear"         , 0.16);
    q0    = setup->get("GeometrySA.SafetyFactor"  , 1.5 );
    eps   = setup->get("GeometrySA.AspectRatio"   , 0.19);
    alpha = setup->get("GeometrySA.ShafranovShift", 0.0 );
   
  };

 
  inline const double B(const int x, const iny y, const int z) {
    return 1./(1 + eps * cos(Z(z))); 
  };

  inline const double dB_dz(const int x, const iny y, const int z) { return - Bhat2 * eps/(q0 * R0) * cos(Z(z)); };
  inline const double dB_dx(const int x, const iny y, const int z) { return - Bhat2 / R0 * cos(Z(z)); };

  inline const double J(const int x, const int y, const int z) {

//    const double g_33 = 1.;

    return 1.;
    // where does this comes from
    return 1./sqrt( g_33(z) * ( g_11() * g_22(z) - pow2(g_12(z))));

  };


   /** Nabla 2 operator  */
  const double k2_p(const int x_k, const int y_k, const int z) {
      
      const double kx2 = pow2(k(Nx,Lx,x_k));
      const double ky2 = pow2(k(Ny,Ly,y_k));

      return g_11() * kx2 + g_22(z) * ky2 + 2.* g_12(z) * sqrt(kx2 + ky2);
  };

private:
  
  /** The metric coefficient the individual metric components
   *
   *    The componets of the s-a metric is following
   *
   *    NOTE : It's a symmetric metric !
   *
   * */
  inline const double g_11()      { return 1.;};
  inline const double g_12(int z) { return shear/(q0*R0)*Z(z);};
  inline const double g_22(int z) { return 1.  + pow2(shear/(q0*R0) * Z(z)) + pow2(eps*q0); };
  inline const double g_23(int z) { return q0*eps; };
  inline const double g_33(int z) { return pow2(q0*eps); };


protected:

   virtual void printOn(ostream& output) const {
         output   << "Geometry  |  s-alpha Geometry shear : "  << shear  << " q0 : " << q0 << " R0 : " << R0 << std::endl;
         output   << "          |  Shfaranov-Shift : " << alpha << "Epsilon_trapped : " << epsilon  << std::endl;
   };


};


#endif // GEOMETRY_H
