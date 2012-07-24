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

#include "Vlasov.h"

/**
 *
 *   @brief Calculated global transport rates
 *
 *   @note Not working
 */
class Transport {

   Transport(Setup *setup, Grid *grid, Parallel *parallel, FileIO *fileIO, Analysis *analysis);

  
  
  void update(Vlasov *vlasov);
protected:

        virtual void printOn(ostream &output) const {

            output << "Transport " << std::endl;

        }


  void initDataOutput(Setup *setup, FileIO *fileIO) {};
  virtual void writeData(Timing timing, double dt) {};
  void closeData() {};

};



#endif
