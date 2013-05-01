/*
 * =====================================================================================
 *
 *       Filename: Event.h
 *
 *    Description: Initiates events started by triggers
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __EVENT_H
#define __EVENT_H

#include "Global.h"

#include "Setup.h"
#include "Parallel/Parallel.h"
#include "Grid.h"
#include "Geometry/Geometry.h"

#include "FileIO.h"

#include "Timing.h"
#include "Plasma.h"

#include "Vlasov/Vlasov.h"
#include "Vlasov/Vlasov_Cilk.h"
#include "Fields/Fields.h"


/**
*  
*  @brief support of timed events
*
*  checkEvent() is called every time step. In case of an event, e.g.
*  suppress zonal-flow from T > 100 , this should be written here.
*
*  @todo provide example
*
**/
class Event : public IfaceGKC {

  bool useEvent; //< Enabled in case events are enabled
public:
  Event(Setup *setup, Grid *grid, Parallel *parallel, FileIO *fileIO, Geometry *geo);


  void checkEvent(Timing timing, Vlasov *vlasov, Fields *fields);
        
  virtual void printOn(std::ostream &output) const {
         output   << "Event     |  " << ((useEvent ? "On" : "Off")) << std::endl;
  }
};


#endif // __EVENT_H
