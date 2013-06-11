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
*  where \f$L_s\f$ is defined as the connection length. The metric takes
*  the simple form 
*
*  \f[
*     
*     g^{ij} = 
*     
*       \left( \begin{array}{lll}
*                1     &  z / L_s                 & 0 \\
*                z / L_s &  1 + \frac{z^2}{L_s^2} & 0 \\
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
*  with \f$ \hat{s} \f$ the shear. For the parallel length \f$ L_z \f$, we need to choose the connection length 
*  L_c. So once L_s and \f$ L_z \f$ is given we calculate the shear to
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
*  @section Matrix Derivation of the matrix coefficients
*
*  In the usual case the safety profile is not constant but changes along the radial direction. This leads
*  to a shearing of the magnetic field lines as shown in Fig.(\ref{Fig:Geometry_FluxSurfaces}).
*  The magnetic field lines are not straight anymore, but tilted as given by
*  \f[
*    \frac{\vec{B}_0}{B_0} = \vec{b} = \left( \begin{array}{l} 0 \\ -\frac{x}{L_s} \\ 1 \end{array} \right) \quad, 
*  \f]
*  where \f$ L_s \gg x \f$ is the corresponding shearing length and \f$ L_z \gg L_x,L_y \f$. 
*  Here straight lines are now only find at \f$ x_0 \f$, where along radial \f$ x \f$ the tilting increases. 
*  Although we could, simply set angle of the magnetic field one - 
*  \cite{Antonsen_1980:KELowFreq} discussed that the turbulence structure develops 
*  along the magnetic field line.
*  Thus as discussed by \cite{Jenko_1998:PhDThesis} it is recommended to change our coordinate system
*  in such a way, that 
*  to align our coordinated system to the magnetic field so that it appears straight in the new coordinates,
*  so that the sheared magnetic field in Eq.(\ref{Eq:Geometry_ShearedSlab}) can be written as 
*  \f[
*    \frac{\vec{B'}_0}{B'_0} = \vec{b}' = \left( \begin{array}{l} 0 \\ 0 \\ 1 \end{array} \right) \quad,
*  \f]
*  in the new coordinate system. This is favorable, as the turbulence is know to extend along the magnetic
*  field lines. With a field-aligned coordinate system, we can greatly reduced then number of required grid points
*  needed in numerical simulations in order to resolve  turbulence. The desired coordinate transformation 
*  from slab geometry ( Cartesian) to sheared slab geometry is
*  \f{align}{
*    x' &= x \\
*    y' &= y + \frac{xz}{L_s} \\
*    z' &= z \quad,
*  \f}
*  where \f$ x,y,z \f$ denotes the old coordinate system and \f$ x',y',z' \f$ denotes the new coordinate system. 
*
*  @subsection Metric The metric
*
*  The shape of the Laplacian in arbitrary local coordinates can be derived from the metric tensor (roughly speaking the 
*  metric tensor tells us how distances are measured in a coordinate system according to 
*  \f$ ds^2 = g^{\mu\nu} dx_\mu dx_\nu \f$), which is calculated according to the following rule
*  \f[
*    g_{\text{new}}^{ab} = \frac{\partial u^a}{\partial v^k} g_{cart}^{kl} \frac{\partial u^b}{\partial v^l} \quad,
*  \f]
*  where \f$ g_{\text{cart}} \f$, the metric tensor in Cartesian coordinates, is defined as 
*  \f[
*    g_{cart} =  \left( \begin{array}{lll} 1 & 0 & 0 \\ 0 & 1 & 0  \\ 0 & 0 & 1  \end{array} \right)  = \delta^{kl} \quad.
*  \f]
*  Thus we need to calculate
*  \f[
*    g_{\text{new}}^{ab} = \frac{\partial u^a}{\partial v^k} \frac{\partial u^b}{\partial v^k} \quad,
*  \f]
*  where \f$ (u^x,u^y,u^z)=(x',y',z') \f$ and \f$ (v^2,v^2,v^3)=(x,y,z) \f$.
*  It can bee seen that the metric is symmetric, thus we need to calculate only the following terms
*  \f{align}{
*  g^{xx} &= \frac{\partial x'}{\partial x}\frac{\partial x'}{\partial x} +  \frac{\partial x'}{\partial y}\frac{\partial x'}{\partial y} +  \frac{\partial x'}{\partial z}\frac{\partial x'}{\partial z}  = 1 \\
*  g^{xy} &= \frac{\partial x'}{\partial x}\frac{\partial y'}{\partial x} +  \frac{\partial x'}{\partial y}\frac{\partial y'}{\partial y} +  \frac{\partial x'}{\partial z}\frac{\partial y'}{\partial z}  = \frac{z}{L_s} \\
*  g^{xz} &= \frac{\partial x'}{\partial x}\frac{\partial z'}{\partial x} +  \frac{\partial x'}{\partial y}\frac{\partial z'}{\partial y} +  \frac{\partial x'}{\partial z}\frac{\partial z'}{\partial z}  = 0 \\
*  g^{yy} &= \frac{\partial y'}{\partial x}\frac{\partial y'}{\partial x} +  \frac{\partial y'}{\partial y}\frac{\partial y'}{\partial y} +  \frac{\partial y'}{\partial z}\frac{\partial y'}{\partial z}  = \left( \frac{z}{L_s} \right)^2 + 1 + \left( \frac{x}{L_s} \right)^2 \\
*  g^{yz} &= \frac{\partial y'}{\partial x}\frac{\partial z'}{\partial x} +  \frac{\partial y'}{\partial y}\frac{\partial z'}{\partial y} +  \frac{\partial y'}{\partial z}\frac{\partial z'}{\partial z}  = 0 \\
*  g^{zz} &= \frac{\partial z'}{\partial x}\frac{\partial z'}{\partial x} +  \frac{\partial z'}{\partial y}\frac{\partial z'}{\partial y} +  \frac{\partial z'}{\partial z}\frac{\partial z'}{\partial z}  = 1 \quad,
*  \f}
*  neglecting higher order terms results in our metric tensor in sheared slab geometry.
*  \f[ 
*    g_{shear} =  \left( \begin{array}{lll} 1 & \frac{z}{L_s} & 0 \\ \frac{z}{L_s} & 1 + z^2  & 0  \\ 0 & 0 & 1  \end{array} \right) \quad.
*  \f]
*  The metric tensor above is non-orthogonal, as can be seen by the non-diagonal elements.
*
*  @subsection Laplacian
*
*  In slab geometry the Laplacian takes the well known form of \f$ \nabla_\perp^2 = \partial_{xx}^2 + \partial_{yy} \f$. The form of
*  the Laplacian is not coordinate invariant, and needs to be calculated.
*  For a general  curvilinear coordinate system, the Laplacian operator has the form 
*  \f[
*    \Delta f = \frac{1}{\sqrt{|g|}} \partial_i \left( \sqrt{|g|} g^{ij} \partial_j f \right) \quad.
*  \f]
*  For the above case Einstein's summation rule apply (sum over equal indices) and \f$ |g| \f$ is the determinant of the metric. For polar coordinates
*  the metric is defined as 
*  \f[
*    g_{pol}^{ab} = \left( \begin{array}{ll} 1 & 0 \\ 0 & r^{-2} \end{array} \right) \quad.
*  \f]
*  Calculating the determinant \f$ |g_{pol}| = r^{-2} \f$ and using above rules, we get for the Laplacian in polar coordinates as
*  \f[
*    \nabla_{pol} f = \partial_i \left(g^{ij} \partial_i \right) = g^{xx} \partial_x^2 + g^{yy} \partial_y^2 + 2 g^{xy} partial_x \partial_y \quad.
*  \f]
*
*  For our sheared coordinate system, using the sheared slab metric, with the determinant  \f$ |g_{\text{shear}}|=1 \f$ and above rule
*  we get
*  \f{align}{
*    \nabla_{sh} f &= \partial_i \left(g^{ij} \partial_i \right) = g^{xx} \partial_x^2 + g^{yy} \partial_y^2 + 2 g^{xy} partial_x \partial_y  \\
*                  &= \partial_{x}^2 + \left( 1+ (\frac{z}{L_s})^2 \right) \partial_{yy}^2 + 2 \frac{z}{L_s} \partial_x \partial_y
*  \f}
*  or in Fourier space as
*  \f[
*    k_\perp^2 = k_x^2 + \left( 1+ (\frac{z}{L_s})^2 \right) k_y^2 + 2 \frac{z}{L_s} k_x k_y
*  \f]
*  Thus for our model
*  \f{align}{
*    \vec{b} \cdot \nabla &= \frac{\partial}{\partial z} \\
*    \nabla_\perp^2       &=  \frac{\partial^2}{\partial^2 x} + \left(1 + \frac{z^2}{L_s^2} \right) \frac{\partial^2}{\partial^2 y}  + 2 \frac{z}{L_s} \frac{\partial^2}{\partial x}{\partial y}\\
*    \vec{b} \times \nabla A \cdot \nabla &= \frac{\partial A}{\partial x} \frac{\partial}{\partial y} - \frac{\partial A}{\partial y}{\partial}{\partial x}
*  \f}
*
*  @subsection boundary Boundary condition
*  For the boundary condition we need to connect the field lines. So that the field lines at \f$ z=L_z \f$ connects to the same field lines
*  at \f$ z=0 \f$. This is not trivial as can be seen in Fig. (\ref{Fig:Geometry_3DShearedSlab}), instead a shift is necessary in \f$ y \f$ direction.
*  Thus boundary condition becomes,
*  \f{align}{
*    A(x,y,z) &= A(x+L_x, y, z)  \quad,\\ 
*    A(x,y,z) &= A(x,y + L_y, z) \quad,\\
*    A(x,y,z) &= A(x,y-2\pi\hat{s}x, z+L_z) \quad,
*  \f}
*  where we used the definition of the magnetic shear defined as \f$ \hat{s} = \frac{1}{2\pi L_s} \f$. Note that as 
*  as \cite{Jenko_1998:PhDThesis} and \cite{Dorland_1988:PhDThesis} noted that this requires a
*  quantization equation \f$ \hat{s} = \frac{n_s}{2\pi} \frac{L_y}{L_x} \f$ of the box sizes so that the magnetic field
*  lines connect to themselves.
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
    if((std::abs(fmod(2. * M_PI * shear * Lx/Ly, 1.)) > 1.e-5) && roundShearToConnect) {

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
