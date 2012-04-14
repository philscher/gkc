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


enum FFT_DIR    {FFT_NO=0, FFT_FORWARD=1, FFT_BACKWARD=2 };
enum FFT_FLAGS  {FFT_DUMMY=0, FFT_XYZ=1, FFT_X=2, FFT_XY=4, FFT_Y=16, FFT_AA=32, FFT_FIELDS=64};


// note we only support r2c, you can add simply suport for c2c as parallelization efficiency
// is reduced due to unneccesary data, we don't recommend.



/** Wrapper to various FFT transform libraries
 * 
 *
 *  NOTE : Rea2COmplex are only supported (due to efficenciy)
 *         in-place transformation are allowed, in this case
 *         set r3In = k2Out and k2In = r3Out to point at the
 *         same memory region.
 *
 *
 */


class FFTSolver : public IfaceHelios {

protected:
  Geometry<HELIOS_GEOMETRY> *geo;
  Parallel *parallel;

  int flags;

  /**
   *
   *    Get Normalization of FFT Solver. Most of the solvers give only normalization
   *    factors for forward AND backward transformation. However, we require to know
   *    these terms seperately in case of multiplication of two terms in real space.
   *    (e.g. due to non-linear term).
   *
   *
   *
   *
   *
   * */
  void checkNormalization();

public:

  /**    Get perpendicular gradient in Fourier space. To to non-rectangular 
   *     coordiantes (shear) we have to include the non-diagonal metric component
   *     g12, g21 too
   *
   *     \f[ k_\perp^2 = k_x^2 + k_y^2 \f]  
  */
  inline double k2_p(const int x_k, const int y_k, const int z) 
  {
      
      const double kx_ = kx(x_k);
      const double ky_ = ky(y_k);

      return geo->g_xx(z) * pow2(kx_) + geo->g_yy(z) * pow2(ky_) + 2. * geo->g_xy(z) * kx_ * ky_;
  };

  /**
   *  Gives \f[ k_x = \frac{2\pi}{L_x} * x_k \f] Fourier wavenumber.
   *  We use fttw ordering ([k_0, k_1, ..., k_Nyq, -k_{Ny-1}, -k_{Ny-2}, ...., -k_{-1} ]
   */
  inline double  kx(const int x_k) { return 2.*M_PI/Lx * ((x_k <= Nx/2) ? x_k : x_k - Nx); }

  	//else          return ( (x_k <= Nx  ) ? ((double) x_k)*2.e0*M_PI/Lx : ((double) (2*Nx - x_k))*2.e0*M_PI/Lx);
  
  /**
   *  Gives \f[ k_y = \frac{2\pi}{L_y} * y_k \f] Fourier wavenumber
   *  We use fttw ordering ([k_0, k_1, ..., k_Nyq, -k_{Ny-1}, -k_{Ny-2}, ...., -k_{-1} ]
   *  Half modes in k_y due to hermitian symmetry. We use c2r transform thus only positive
   *  values are reqtuied
   */
 inline double  ky(const int y_k) { return 2.*M_PI/Ly * y_k; };


  /**  Input Array for 3 dimensional data  */
  Array4z r3In, r3Out;
  Array3z r33In, r33Out;

// 3-Dimensional FFT
  Range Rk3xL, Rk3yL, Rk3zL;
  int k3NxL, k3NyL, k3NzL, k3NxG, k3NyG, k3NzG;
  int K3xLuD, K3yLuD, K3zLuD, K3xLlD, K3yLlD, K3zLlD;
  Array3z k3In, k3Out;

  /** 2 dimensional FFT needed e.g. for poisson equation \f[ \phi(x,y)  -> \phi(x_k, y_k) \f]
   * Note : Need also to support stacked transformations
   *        Arrays can share a common points as only one active FFT is allowed
   */ 
  Range Rk2xL, Rk2yL;
  int K2xLuD, K2yLuD, K2xLlD, K2yLlD;
  Array4z rXYIn, rXYOut, kXYIn, kXYOut;

  /** Normalization factors in N3(Nx,Ny,Nz), N2(Nx,Ny) and Norm_X(Nx)  */
  double Norm_XYZ, Norm_XY, Norm_X, Norm_Y;
 
  double Norm_X_Forward, Norm_X_Backward;
  double Norm_Y_Forward, Norm_Y_Backward;

  int getFlags() const { return flags; };
  
 /*  FFT in Y - direction */ 
 Array4d rYIn, rYOut;
 Array4z kYOut, kYIn;
 int Y_kyLlD, Y_kyLuD , Y_NyLlD , Y_NyLuD;
 Range Y_RkyL , Y_RyLD;




  Array4z rXIn, rXOut, kXOut, kXIn;
  Array4z r2In, r2Out, k2In, k2Out;
  int K1xLlD, K1xLuD;
  int K1yLlD, K1yLuD;
  Range Rk1xL, Rk1yL;


   FFTSolver(Setup *setup, Parallel *_parallel, Geometry<HELIOS_GEOMETRY> *_geo, double _Norm_XYZ, double _Norm_XY, double _Norm_X, double _Norm_Y) :
      parallel(_parallel), Norm_XYZ(_Norm_XYZ), Norm_XY(_Norm_XY), Norm_X(_Norm_X), Norm_Y(_Norm_Y), geo(_geo)
  {
     flags = FFT_X | FFT_Y;
     if(setup->get("FFTSolver.XYZ", 0 ) == 1)              flags |= FFT_XYZ; 
     if(setup->get("Vlasov.useAA", 0) == 1)              flags |= FFT_AA;
  };

  virtual ~FFTSolver() {};

  // this libray needs to support 1D FFT and XYZ-FFT
  virtual int solve(const int type, const int direction, const int nstacked=1) = 0;

  virtual string getLibraryName() { return "No"; }; // = 0;
 // virtual int    getDecomposition() { return DECOMP_NO; };

    /**
     *
     *  A Input Array 1
     *  B
     *  out : R Results
     *
     * */
 virtual Array3z multiply(Array3z &A, Array3z &B, Array3z  &R) = 0;

protected :
     virtual void initDataOutput(FileIO *fileIO) {};
     virtual void writeData(Timing *timing) {};
     virtual void closeData() {};
   


 void printOn(ostream &output) const {
//            output << "Libraries |  FFT : " << getLibraryName() << "      ";


    }
};


#endif // __FFTSolver_H
