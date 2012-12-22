/*
 * =====================================================================================
 *
 *       Filename: Grid.h
 *
 *    Description: Grid definitions incuding boundary
 *
 *         Author: Paul P. Hilscher (2009-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __GRID_H_
#define __GRID_H_

#include "Global.h"

#include "Setup.h"
#include "Parallel/Parallel.h"

#include "FileIO.h"

/** 
*    @brief Size of the computational domain for \f$x,k_y,z, v_\parallel, \mu, \sigma \f$
*    @image html Grid_BoundaryDomain.png
*
*    We treat X-direction in real space, while the poloidal direction is expanded into
*    Fourier harmonics, like i
*
*    --- Poloidal direction ---
*
*    The poloidal direction (y) is expanded into its Fourier harmonics
*
*    \f[
*       F(y_k) = sum_{n=0}^{N_y-1} f(y) exp(- 2\pi j  \frac{k}{N_y}{y_k} \quad.
*    \f]
*
*    (see http://en.wikipedia.org/wiki/Discrete_Fourier_transform). 
*    
*    espcially, we note that the DC comonent F(y_k=0)  = int_ f(y) dy.
*    The Nyquiest frequency (Nyquist : y_k=Nky/2) is removed as it has
*    no phase information. Also as the field quantities are real valued
*    only positive modes (y_k >= 0) are evolved, as the negative modes
*    are simply the complex conjugates and can easily be calculated.
*
*
*
**/
class Grid : public IfaceGKC {

 public:

  /// @name Number of ghost cells
  ///@{
  int NxGC, ///< Number of ghost cells in x-directions
      NyGC, ///< Number of ghost cells in y-direction
      NzGC, ///< Number of ghost cells in z-direction
      NvGC; ///< Number of ghost cells in v-direction
  ///@}

  double *dm; ///< Weights for \f$ \mu(m) \f$
  
  /// @name Ranges for computational domain
  ///@{
  nct::Range RzLD, RxLD, RvLD, RmLD, RsLD, RkyLD; 
  nct::Range RxLB, RzLB, RvLB, RmLB, RsLB;
   
  nct::Range RxGB, RzGB, RvGB, RmGB, RsGB;
  nct::Range RxLB4;

  ///@}
   
  /// @name Ranges for computational domain
  ///@{
  nct::Range RxGD, RkyGD, RzGD, RvGD, RmGD, RsGD; 
  ///@}

  /** \f$ textrm{d}x textrm{d}y textrm{d}z \f$
  *   \f$ textrm{d}x textrm{d}y textrm{d}z textrm{d}v_\parallel  \f$
  *
  *    @note will depend on geometry
  **/
  double dXYZ, dXYZV;

  /// Local Grid size (for Domain & Boundary)
  int NxGD, NyGD, NkyGD, NzGD, NvGD, NmGD, NsGD;

  /// Constructor
  Grid(Setup *setup, Parallel *parallel, FileIO *fileIO);

  /// Destructor
 ~Grid();

 protected:

  /// Class information
  virtual void printOn(std::ostream &output) const;

 public: 
   
  /**
  *    @brief saves grid information to output file
  *
  *
  *    @param fileIO class
  **/
  void initData(FileIO *fileIO);
   
  /// Class information
  virtual void writeData(const Timing &timing, const double dt) {};
   
  /// Class information
  virtual void closeData() {};

  ///  returns total global domain size (@todo accept DIR)
  int getGlobalSize() const { return Nx * Nky * Nz * Nv * Nm * Ns; };
   
  ///  returns total local domain size (@todo accept DIR)
  int getLocalSize() const { return NxLD * NkyLD * NzLD * NvLD * NmLD * NsLD;};

};

#endif // __GRID_H_

