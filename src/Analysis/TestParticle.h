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


#ifndef _GKC_TESTPARTICLE_H__
#define _GKC_TESTPARTICLE_H__

#include "Global.h"
#include "Setup.h"
#include "FileIO.h"

#include "Parallel/Parallel.h"
#include "Vlasov/Vlasov.h"
#include "Fields/Fields.h"

#include "Special/Vector3D.h"

/**
* @brief Passive trace 
*
* @todo this class has no function yet
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
      int id;

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
  
  virtual void initData(Setup *setup, FileIO *fileIO);
  virtual void writeData(Timing timing, double dt);
  virtual void closeData();
}

#endif // _GKC_TESTPARTICLE_H__
