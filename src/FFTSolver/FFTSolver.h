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
#include "Parallel/Parallel.h"
#include "Setup.h"
#include "config.h"
#include "Geometry/Geometry.h"



/**
*  @brief sign of transform (forward or backward)
*  @todo  avoid global scope 
**/
enum class FFT_Sign : int { No=0, Forward=1, Backward=2 };

/**
*  @brief type of transform
*  @todo  avoid global scope 
**/
enum class FFT_Type : int {DUMMY=0, XYZ=1, X=2, XY=4, Y=16, AA=32, FIELDS=64, X_FIELDS=128, Y_FIELDS=256, Y_PSF=512, Y_NL=1024};


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
class FFTSolver : public IfaceGKC {

  std::vector<int> suppressModeX, suppressModeY;
    
   /**  We can suppress various modes, this is set in the setup of the fields.
   *   and sets the Fourier mode to zero. Move to FFT solver.
   */
//   void suppressModes(CComplex kXOut[Nq][NzLD][NkyLD][FFTSolver::X_NkxL]); 
   void parseSuppressMode(const std::string &value, std::vector<int> &suppressMode);

  protected:

   Geometry *geo;
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

      return geo->g_xx(0, z) * pow2(kx_) + geo->g_yy(0, z) * pow2(ky_) + 2. * geo->g_xy(0, z) * kx_ * ky_;
   };

   /**
   *  Gives \f[ k_x = \frac{2\pi}{L_x} * x_k \f] Fourier wavenumber.
   *  We use fttw ordering ([k_0, k_1, ..., k_Nyq, -k_{Ny-1}, -k_{Ny-2}, ...., -k_{-1} ]
   **/
   static inline double  kx(const int x_k) { return 2.*M_PI/Lx * ((x_k <= Nx/2) ? x_k : x_k - Nx); }

   /**
   *  Gives \f[ k_y = \frac{2\pi}{L_y} * y_k \f] Fourier wavenumber
   *  We use fttw ordering ([k_0, k_1, ..., k_Nyq, -k_{Ny-1}, -k_{Ny-2}, ...., -k_{-1} ]
   *  Half modes in k_y due to hermitian symmetry. We use c2r transform thus only positive
   *  values are reqtuied
   */
   static inline double  ky(const int y_k) { return 2.*M_PI/Ly * y_k; };


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
   *  @brief Arrays for FFT in X-direcction
   *
   *  Transforms the field equations \f$ A(x,k_y, z, n) rightarrow A(k_x, k_y, z,n) \f$,
   *  which is needs on e.g. the Poisson's equation.
   * 
   *  Output of the FFT Solver. Note that the FFT solver is responsible for
   *  correct (re-) ordering of the output data, as well as input data,
   *  e.g. transpositions.
   *
   *  Real space data (in x) is directly passed as the coresponding Fields0, Q, Qm.
   *  
   *
   **/
   CComplex *kXIn, *kXOut;

   static int X_NkxL;
   int K1xLlD, K1xLuD;
   int K1yLlD, K1yLuD;
   // @}

   /**
   *   @brief the contstructor
   *
   **/
   FFTSolver(Setup *setup, Parallel *_parallel, Geometry *_geo, double _Norm_XYZ, double _Norm_XY, double _Norm_X, double _Norm_Y); 

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
   virtual void solve(const FFT_Type type, const FFT_Sign direction, void *in=nullptr, void *out=nullptr) = 0;

   /**
   *   @brief  returns the FFT library name in use
   *
   *   @return the name with version of the FFT library used
   **/
   virtual std::string getLibraryName()  = 0;

   /**
   *  @brief multiplies 2-arrays by transforming to real space
   *  A Input Array 1
   *  B
   *  out : R Results
   *  
   *  The final result is properly rescaled.
   *
   *  @param A  Array A(x,y_k,z)
   *  @param B Array B(x_y_k,z)
   *  @param R Result C(x,y_k,z) = A B
   *
   **/
   virtual void multiply(const CComplex A[NkyLlD][NxLD], const CComplex B[NkyLlD][NxLD],
                               CComplex R[NkyLlD][NxLD]) = 0;
  
};


#endif // __FFTSolver_H
