/*
 * =====================================================================================
 *
 *       Filename: Geometry.h
 *
 *    Description: Geometry interface
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#ifndef GEOMETRY_H
#define GEOMETRY_H


#include "Global.h"

#include "Setup.h"
#include "Grid.h"
#include "FileIO.h"

/** 
*   @brief Base class for the Geometry abstraction
*
*
**/
class Geometry : public IfaceGKC
{

 public:

  double *Kx,     ///< Curvature term
         *Ky,     ///< Curvature term 
         *B,      ///< Equilibrium magnetic field strength
         *dB_dx,  ///< $\frac{\partial B}{\partial x}$
         *dB_dy,  ///< $\frac{\partial B}{\partial y}$ 
         *dB_dz,  ///< $\frac{\partial B}{\partial z}$
         *J    ;      ///< $J$
//         *C    ; ///< C


  nct::allocate ArrayGeom;

  double eps_hat, C;

  Geometry(Setup *setup, Grid *grid, FileIO *fileIO) {

    ArrayGeom = nct::allocate(grid->RzLD, grid->RxLD);
    //ArrayGeom(&Kx, &Ky, &B, &dB_dx, &dB_dy, &dB_dz, &J, &C);
    ArrayGeom(&Kx, &Ky, &B, &dB_dx, &dB_dy, &dB_dz, &J);
    
    eps_hat = setup->get("Geometry.eps", 1.);

    C = 1./(dx * dy * dz);
  };

  virtual ~Geometry() {};

 protected:

  /**
  *    @brief setup the geometry arrays, called from derived class
  *
  *    @note This function is necessary, as calling pure virtual functions
  *          from the base constructor is not allowed. 
  *          See http://stackoverflow.com/questions/99552/
  *
  *    Set virtual, as for efficiency reason can be set directly when
  *    numerical MHD equilibrium is used
  *
  **/ 
  virtual void setupArrays() {

    // Setup Arrays
    [=](double Kx   [NzLD][NxLD], double Ky   [NzLD][NxLD], double B[NzLD][NxLD], double dB_dx[NzLD][NxLD], 
        double dB_dy[NzLD][NxLD], double dB_dz[NzLD][NxLD], double J[NzLD][NxLD])
    {
        // set metric elements
      for(int x=NxLlD; x <= NxLuD; x++) {  for(int z=NzLlD; z <= NzLuD; z++) { 

//      C    [z][x] = get_C    (x,z);
        Kx   [z][x] = get_Kx   (x,z);
        Ky   [z][x] = get_Ky   (x,z);
        B    [z][x] = get_B    (x,z);
        dB_dx[z][x] = get_dB_dx(x,z);
        dB_dy[z][x] = get_dB_dy(x,z);
        dB_dz[z][x] = get_dB_dz(x,z);
        J    [z][x] = get_J    (x,z);

      } }
    } ((A2rr) Kx, (A2rr) Ky, (A2rr) B, (A2rr) dB_dx, (A2rr) dB_dy, (A2rr) dB_dz, (A2rr) J);

  }

 public:

  // Jacobian
  virtual double get_J(const int x, const int z)  = 0;
   
  /**
  *   @name Normalized Curvature terms
  *   Get the value of shear at position x
  *   what is the definition of it ?
  **/
  ///@{
  
  /**
  *   Helper functions defined as
  *   \f[
  *     K_x = -\frac{1}{C} \frac{L_{ref}}{B_{ref}} 
  *             \left( 
  *             \frac{\partial B_0}{\partial_y} 
  *             - \frac{\gamma_2}{\gamma_1}  \frac{\partial B_0}{\partial z} 
  *             \right)
  *   \f]
  * 
  *   see Goerler, PhD
  *   see Lapillone, PhD Thesis, Eq.(2.135)
  *
  **/ 
  virtual  double get_Kx(const int x, const int z)  
  { 
    return - 1./C * ( g_2(x,z)/g_1(x,z) * get_dB_dz(x,z));
  };

  /**
  *   Helper functions defined as
  *
  *   \f[
  *     K_y = \frac{1}{C} \frac{L_{ref}}{B_{ref}} 
  *             \left( 
  *             \frac{\partial B_0}{\partial_x} 
  *             - \frac{\gamma_3}{\gamma_1}  \frac{\partial B_0}{\partial z} 
  *             \right)
  *   \f]
  * 
  *   see Goerler, PhD
  *   see Lapillone, PhD Thesis, Eq.(2.135)
  *
  **/ 
  virtual double get_Ky(const int x, const int z)  
  {
    return 1./C * ( get_dB_dx(x,z) - g_3(x,z)/g_1(x,z) * get_dB_dz(x,z));
  };
  ///@} 
  
  /**
  *   @name Metric Gamma helper functions
  *
  **/
  ///@{
  /**  Helper function
   *
   *   \f[ \gamma_1 = g_{xx} g_{yy} - g_{xy} g_{yx} \f]
   */ 
  inline  double g_1(const int x, const int z) { return g_xx(x,z) * g_yy(x,z) - g_xy(x,z) * g_xy(x,z); };

  /**  Helper function
   *
   *   \f[ \gamma_2 = g_{xx} g_{yz} - g_{xy} g_{xz} \f]
   */ 
  inline  double g_2(const int x, const int z) { return g_xx(x,z) * g_yz(x,z) - g_xy(x,z) * g_xz(x,z); };

  /**  Helper function
   *
   *   \f[ \gamma_3 = g_{xy} g_{yz} - g_{yy} g_{xz} \f]
   */ 
  inline  double g_3(const int x, const int z) { return g_xy(x,z) * g_yz(x,z) - g_yy(x,z) * g_xz(x,z); };
  ///@}
  
  /**
  *   @name Metric elements
  *
  *   Defines the metric elements as
  *   
  *   \f[
  *      g = \left( \begin{array}{ccc}
  *           
  *           g_{xx} & g_{xy} & g_{xz} \\
  *                  & g_{yy} & g_{yz} \\
  *                  &        & g_{zz} \\
  *          
  *          \end{array} \right)
  *   \f]
  * 
  *   Note, that the metric is symmetric
  **/
  /// @{
  virtual  double g_xx(const int x, const int z) = 0;
  virtual  double g_xy(const int x, const int z) = 0;
  virtual  double g_xz(const int x, const int z) = 0;
  virtual  double g_yy(const int x, const int z) = 0;
  virtual  double g_yz(const int x, const int z) = 0;
  virtual  double g_zz(const int x, const int z) = 0;
  /// @}
  
  
  /**
  *   @name Magnetic field variations
  *
  **/ 
  ///@{
  /// Defined Magnetic field strength \f[ B_0(\vec{x}) \f] at position \f[ \vec{x} = (x,y,z) \f]
  virtual  double get_B    (const int x, const int z) = 0;  
  
  ///  Defined Magnetic field strength \f[ \frac{\partial B_0}{\partial x}(\vec{x}) \f] at position \f[ \vec{x} = (x,y,z) \f]
  virtual  double get_dB_dx(const int x, const int z) = 0;

  ///  Defined Magnetic field strength \f[\frac{\partial B_0}{\partial y}(\vec{x}) \f] at position \f[ \vec{x} = (x,y,z) \f]
  virtual  double get_dB_dy(const int x, const int z) = 0;
  
  ///  Defined Magnetic field strength \f[\frac{\partial B_0}{\partial z}(\vec{x}) \f] at position \f[ \vec{x} = (x,y,z) \f]
  virtual  double get_dB_dz(const int x, const int z) = 0;
  ///@}
  
  
  /**
  *   @brief factor to connect parallel field line at border
  *
  *   \f[
  *        A(x,k_y,z) = f(x,k_y,z+2\pi) e^{i 2\pi \nu(x)} 
  *   \f]
  *
  **/
  virtual  double nu (const int x) = 0;

  
  /**
  *    @brief creates Geometry group and delegates to child
  *
  *
  **/ 
  void initData(FileIO *fileIO) {

    hid_t geometryGroup = fileIO->newGroup("Geometry");

    initData(geometryGroup); 
    H5Gclose(geometryGroup);

  };
        
   
  virtual void initData(hid_t geometryGroup) = 0;

  // is it necessary for later time steps to 
  virtual void writeData(Timing *timing) {};
  virtual void closeData() {};

};

#endif // GEOMETRY_H
