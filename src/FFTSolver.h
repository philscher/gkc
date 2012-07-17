/*
 * =====================================================================================
 *
 *       Filename: FFTSolver.h
 *
 *    Description: Interface for various (Fast) Fourier Solvers
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __FFTSolver_H
#define __FFTSolver_H

#include "Global.h"
#include "Parallel.h"
#include "Setup.h"
#include "config.h"
#include "Geometry.h"

#include "GeometrySlab.h"
#include "GeometryShear.h"
#include "Geometry2D.h"



/**
*  @enum sign of transform (forward or backward)
*  @todo avoid global scope 
**/
enum FFT_DIR    {FFT_NO=0, FFT_FORWARD=1, FFT_BACKWARD=2 };

/**
*  @enum type of transform
*  @todo avoid global scope 
**/
enum FFT_FLAGS  {FFT_DUMMY=0, FFT_XYZ=1, FFT_X=2, FFT_XY=4, FFT_Y=16, FFT_AA=32, FFT_FIELDS=64};


/** 
*
*  @brief Interface to various FFT transform libraries
*
*  @note  only real-to-complex transforms are supported, as well
*         as in-place transformation. Take care, that in this
*         case e.g. rXIn = kXOut and kXIn = rXOut will point
*         to the same memory region.
*
**/
class FFTSolver : public IfaceHelios {

  std::vector<int> suppressModeX, suppressModeY;
    
   /**  We can suppress various modes, this is set in the setup of the fields.
   *   and sets the Fourier mode to zero. Move to FFT solver.
   */
   int suppressModes(Array4z k2Out, const int field=1);
   void parseSuppressMode(const std::string &value, std::vector<int> &suppressMode);

  protected:

   Geometry<HELIOS_GEOMETRY> *geo;
   Parallel *parallel;

   int flags;

   /**
   *  @brief check normalization for forward and backward transform
   *    
   *  Get Normalization of FFT Solver. Most of the solvers give only normalization
   *  factors for forward AND backward transformation. However, we require to know
   *  these terms seperately in case of multiplication of two terms in real space.
   *  (e.g. due to non-linear term).
   *  
   *  This is simply checked by performing a forward transformation, and
   *  checking the values, and vice vera.
   *
   *  @todo describe how it does it
   *
   **/
   virtual void setNormalizationConstants();

  public:
   /**    Get perpendicular gradient in Fourier space. To to non-rectangular 
   *     coordiantes (shear) we have to include the non-diagonal metric component
   *     g12, g21 too
   *
   *     \f[ k_\perp^2 = k_x^2 + k_y^2 \f]  
   **/
   inline double k2_p(const int x_k, const int y_k, const int z) 
   {
      
      const double kx_ = kx(x_k);
      const double ky_ = ky(y_k);

      return geo->g_xx(z) * pow2(kx_) + geo->g_yy(z) * pow2(ky_) + 2. * geo->g_xy(z) * kx_ * ky_;
   };

   /**
   *  Gives \f[ k_x = \frac{2\pi}{L_x} * x_k \f] Fourier wavenumber.
   *  We use fttw ordering ([k_0, k_1, ..., k_Nyq, -k_{Ny-1}, -k_{Ny-2}, ...., -k_{-1} ]
   **/
   inline double  kx(const int x_k) { return 2.*M_PI/Lx * ((x_k <= Nx/2) ? x_k : x_k - Nx); }

   /**
   *  Gives \f[ k_y = \frac{2\pi}{L_y} * y_k \f] Fourier wavenumber
   *  We use fttw ordering ([k_0, k_1, ..., k_Nyq, -k_{Ny-1}, -k_{Ny-2}, ...., -k_{-1} ]
   *  Half modes in k_y due to hermitian symmetry. We use c2r transform thus only positive
   *  values are reqtuied
   */
   inline double  ky(const int y_k) { return 2.*M_PI/Ly * y_k; };


   /** Normalization factors in N3(Nx,Ny,Nz), N2(Nx,Ny) and Norm_X(Nx)  */
   double Norm_XYZ, Norm_XY, Norm_X, Norm_Y;
 
   /** Normalization factors for single-X Forward/Backward transformation
    *  
    *  We need to know this, in order to properly rescale a transformation to
    *  real-space, with a multiplication followed by backtransformation.
    *  (multiplty)
    *  
    * */
   double Norm_X_Forward, Norm_X_Backward;
   double Norm_Y_Forward, Norm_Y_Backward;

   int getFlags() const { return flags; };
  
   // @{
   /**
   *  @brief Arrays for FFT in Y-direcction
   *
   *  Transforms the field equations \f$ A(x,k_y, z, n) rightarrow A(x, y, z,n) \f$,
   *  which is needed to calculate the non-linearity in real space.
   *
   *  @note that this is uses a real-to-complex and complex-to-real transform
   **/
   Array4d rYIn, rYOut;
   Array4z kYOut, kYIn;
   int Y_kyLlD, Y_kyLuD , Y_NyLlD , Y_NyLuD;
   Range Y_RkyL , Y_RyLD;
   // @}
   
   // @{
   /**
   *  @brief Arrays for FFT in X-direcction
   *
   *  Transforms the field equations \f$ A(x,k_y, z, n) rightarrow A(k_x, k_y, z,n) \f$,
   *  which is needs on e.g. the Poisson's equation.
   **/
   Array4z rXIn, rXOut, kXOut, kXIn;

   int K1xLlD, K1xLuD;
   int K1yLlD, K1yLuD;
   Range Rk1xL, Rk1yL;
   // @}

   /**
   *   @brief the contstructor
   *
   **/
   FFTSolver(Setup *setup, Parallel *_parallel, Geometry<HELIOS_GEOMETRY> *_geo, double _Norm_XYZ, double _Norm_XY, double _Norm_X, double _Norm_Y); 

   virtual ~FFTSolver() ;

   /**
   *
   *
   *   
   *
   *   @param type       FFT_FLAG enum, e.g. FFT_X, FFT_Y
   *   @param direction  flag either FFT_FORWARD or FFT_BACKWARD
   *   @param nstacked   depreceated
   *
   **/
   virtual int solve(const int type, const int direction, const int nstacked=1) = 0;

   /**
   *   @brief  returns the FFT library name in use
   *
   *   @return the name with version of the FFT library used
   **/
   virtual string getLibraryName()  = 0;

   /**
   *  @multiplies 2-arrays by transforming to real space
   *  A Input Array 1
   *  B
   *  out : R Results
   *  
   *  The final result is properly rescaled.
   *
   *  @param Array A(x,y_k,z)
   *  @param Array B(x_y_k,z)
   *
   *  @returns C(x,y_k,z)
   *
   **/
   virtual Array3z multiply(Array3z &A, Array3z &B, Array3z  &R) = 0;
  
};


#endif // __FFTSolver_H
