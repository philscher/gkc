/*
 * =====================================================================================
 *
 *       Filename: Visualization.h
 *
 *    Description: Writes only slides of potnetial to disk.
 *
 *         Author: Paul P. Hilscher (2011-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

// ToDO : Use Fiedls instead of phi , Ap, Bp


#ifndef VISUALIZATION_H_
#define VISUALIZATION_H_

#include "Global.h"

#include "Setup.h"
#include "Vlasov/Vlasov.h"
#include "Fields.h"
#include "Timing.h"
#include "Grid.h"


/**
*   
*   @brief Visualization of Variables
*
*   Interface for visualization, only a slice of some variables
*   can be visualized - or stored for post-processing.
*
**/
class Visualization {
 
  protected:

   Timing dataOutputVisual;

   Vlasov *vlasov;
   Fields *fields;
   Parallel *parallel;
  public:

   Visualization(Vlasov *_vlasov, Fields *_fields, Grid *_grid, Setup *setup, Parallel *_parallel) 
                : vlasov(_vlasov), fields(_fields), parallel(_parallel) {

    dataOutputVisual      = Timing( setup->get("DataOutput.Visualization.Step", -1),
                                    setup->get("DataOutput.Visualization.Time", -1.));
    
   }   
    
   virtual ~Visualization() {
    
   };

  /**
  *   @brief visualize the data
  *
  *   The arrays can all be accessed from vlasov and fields
  **/
  virtual void writeData(Timing timing, double dt, bool force=false)  = 0;

};


#endif // VISUALIZATION__H_
