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
#include "Parallel.h"
#include "Grid.h"
#include "Geometry.h"
#include "GeometrySlab.h"
#include "GeometryShear.h"
#include "Geometry2D.h"

#include "FileIO.h"

#include "Timing.h"
#include "Plasma.h"

#include "Vlasov.h"
#include "Vlasov/Vlasov_Cilk.h"
#include "Fields.h"


class Event : public IfaceHelios {
  bool useEvent;
public:
  Event(Setup *setup, Grid *grid, Parallel *parallel, FileIO *fileIO, Geometry<HELIOS_GEOMETRY> *geo);


  void checkEvent(Timing timing, Vlasov *vlasov, Fields *fields);
        
  virtual void printOn(ostream &output) const {
         output   << "Event     |  " << ((useEvent ? "On" : "Off")) << std::endl;
  }
};


#endif // __EVENT_H
