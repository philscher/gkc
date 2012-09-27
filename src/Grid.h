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
   
   //Array1R dm; ///< Weights for \f$ \mu(m) \f$
   double *dm;
  
   /// @name Ranges for computational domain
   ///@{
   blitz::Range RxGB, RyGB, RzGB, RvGB, RmGB, RsGB; ///< The boundary domain
   ///@}
   
   /// @name Ranges for computational domain
   ///@{
   blitz::Range RxGD, RyGD, RzGD, RvGD, RmGD, RsGD; 
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
   void initDataOutput(FileIO *fileIO);
   
   /// Class information
   virtual void writeData(Timing *timing) {};
   
   /// Class information
   virtual void closeData() {};

   ///  returns total global domain size (@todo accept DIR)
   int getGlobalSize() const { return Nx * Nky * Nz * Nv * Nm * Ns; };
   
   ///  returns total local domain size (@todo accept DIR)
   int getLocalSize() const { return NxLD * NkyLD * NzLD * NvLD * NmLD * NsLD;};
};



#endif // __GRID_H_

