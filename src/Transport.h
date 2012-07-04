/*
 * =====================================================================================
 *
 *       Filename: Transport.h
 *
 *    Description: Transport model to update Maxwellian
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef TRANSPORT_H
#define TRANSPORT_H

#include "Global.h"

#include "FileIO.h"

#include "Timing.h"
#include "Plasma.h"

#include "Analysis.h"

#include "Geometry.h"
#include "GeometrySlab.h"
#include "GeometryShear.h"
#include "Geometry2D.h"

#include "Vlasov.h"

class Transport {

   Transport(Setup *setup, Grid *grid, Parallel *parallel, FileIO *fileIO, Analysis *analysis);

  
  void initDataOutput(Setup *setup, FileIO *fileIO) {};
  virtual void writeData(Timing timing, double dt) {};
  void closeData() {};

  void update(Vlasov *vlasov);
protected:

        virtual void printOn(ostream &output) const {

            output << "Transport " << std::endl;

        }



};



#endif
