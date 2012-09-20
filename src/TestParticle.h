/*
 * =====================================================================================
 *
 *       Filename: TestParticle.h
 *
 *    Description: Traces Test Particles
 *
 *         Author: Paul P. Hilscher (2011-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#ifndef __TESTPARTICLE_H
#define __TESTPARTICLE_H

#include "Global.h"
#include "Setup.h"
#include "FileIO.h"

#include "Parallel.h"
#include "Vlasov.h"
#include "Fields.h"

#include "Special/Vector3D.h"

/**
* @brief passive tracer particle
*
*
**/
class TestParticles : public IfaceGKC 
{
  Parallel *parallel;
 
  int Total_Tracer;

  struct Particle {
/*
 * Particle() {
            position = velocity =0.;
            mass = 1.;
            charge = 1.;
        }
        */
      int number;

      double mass;
      double charge;
      
      double p[3];
      double v[3];
   };
  
  Particle *particles;
  FileAttr *FA_X, *FA_V, *FA_Time;
  public:


TestParticles(FileIO *fileIO, Setup *setup, Parallel *parallel) ;
~TestParticles() ;
// Simplest model, only follow phi potential
void integrate(Vlasov *vlasov, Fields *fields, int step);


        virtual void printOn(std::ostream &output) const;
      
        virtual void initDataOutput(Setup *setup, FileIO *fileIO);
        virtual void writeData(Timing timing, double dt);
        virtual void closeData();


};




#endif //TESTPARTICLE
