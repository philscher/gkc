/*
 * =====================================================================================
 *
 *       Filename:Geometry.h
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

#include "FileIO.h"

/** 
*   @brief Base class for the Geometry abstraction using CRTP
*
*   For efficiency, we don't
*   use virtual inheritance but CRTP. The latter enables to resolve all
*   inline functions directly and thus should avoid any overhead due to 
*   Geometry abstraction, see
*
*
*  This part mostly relyies on Gorles, PhD, Eq. (2.51)  and goemetric parts
*
* BUG : if u compile and derive and doe not implement one of the function you
*        will get a self-recursion. Any idea how to fix this  ? 
*        See
*        http://stackoverflow.com/questions/4721047/this-is-crtp-usage-for-static-polymorphism-but-without-implementation-of-a-derive
*  @note Do we really see the overhead from virtual functions ? Test implementation
*
*
**/
class Geometry : public IfaceGKC
{

protected:


public:
  Array2d Kx, Ky, B, dB_dx, dB_dy, dB_dz, J;
double eps_hat, C;

  Geometry(Setup *setup, FileIO *fileIO) {

        allocate(RxLD, RzLD, Kx, Ky, B, dB_dx, dB_dy, dB_dz, J);
    
        eps_hat = 1.;
        C       = 1.;

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
  **/ 
  void setupArrays() {
        
        // set metric elements
        for(int x=NxLlD; x <= NxLuD; x++) {  for(int z=NzLlD; z <= NzLuD; z++) { 

             Kx(x,z)    = get_Kx   (x,z);
             Ky(x,z)    = get_Ky   (x,z);
             B(x,z)     = get_B    (x,z);
             dB_dx(x,z) = get_dB_dx(x,z);
             dB_dy(x,z) = get_dB_dy(x,z);
             dB_dz(x,z) = get_dB_dz(x,z);
             J(x,z)     = get_J    (x,z);

       } }


  }

public:

  // Jacobian
  virtual double get_J(const int x, const int z)  = 0;
   
  /**
  *   @name K values  
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
  *
  **/ 
  inline  double get_Kx(const int x, const int z)  { 
    
        return - 1./C * ( dB_dy(x,z) + g_1(x,z)/g_2(x,z) * dB_dz(x,z));
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
  *
  **/ 
  inline  double get_Ky(const int x, const int z)  
  {
    
      return 1./C * ( dB_dx(x,z) - g_3(x,z)/g_1(x,z) * dB_dz(x,z));
  
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
  inline  double g_1(const int x,  const int z) { return g_xx(x,z) * g_yy(x,z) - pow2(g_xy(x,z))          ; };
  /**  Helper function
   *
   *   \f[ \gamma_2 = g_{xx} g_{yz} - g_{xy} g_{xz} \f]
   */ 
  inline  double g_2(const int x, const int z) { return g_xx(x,z) * g_yz(x,z) - g_xy(x,z) * g_xz(x,z)  ; };
  /**  Helper function
   *
   *   \f[ \gamma_3 = g_{xy} g_{yz} - g_{yy} g_{xz} \f]
   */ 
  inline  double g_3(const int x, const int z) { return g_xy(x,z) * g_yz(x,z) - g_yy(x,z) * g_xz(x,z)  ; };
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
  virtual  double get_B   (const int x, const int z   ) = 0;  
  
  ///  Defined Magnetic field strength \f[ \frac{\partial B_0}{\partial x}(\vec{x}) \f] at position \f[ \vec{x} = (x,y,z) \f]
  virtual  double get_dB_dx  (const int x,  const int z) = 0;

  ///  Defined Magnetic field strength \f[\frac{\partial B_0}{\partial y}(\vec{x}) \f] at position \f[ \vec{x} = (x,y,z) \f]
  virtual  double get_dB_dy  (const int x, const int z) = 0;
  
  ///  Defined Magnetic field strength \f[\frac{\partial B_0}{\partial z}(\vec{x}) \f] at position \f[ \vec{x} = (x,y,z) \f]
  virtual  double get_dB_dz  (const int x, const int z) = 0;
  ///@}
  
  
  virtual  double nu (const int x) = 0;

  
   /**
   *    @brief creates Geometry groupd and delegates to child
   *
   *
   **/ 
   void initDataOutput(FileIO *fileIO) {
        hid_t geometryGroup = check(H5Gcreate(fileIO->getFileID(), "/Geometry",H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), DMESG("Error creating group for Geometry : H5Gcreate"));
        initDataOutput(geometryGroup); 
        H5Gclose(geometryGroup);
   } ;
        
   
     virtual void initDataOutput(hid_t geometryGroup) = 0;

     // is it necessary for later timesteps to 
     virtual void writeData(Timing *timing) {};
     virtual void closeData() {};


};



#endif // GEOMETRY_H
