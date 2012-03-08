/*
 * =====================================================================================
 *
 *       Filename: FieldsPETSc.h
 *
 *    Description: Interface to PETSc Matrix Solver routine
 *
 *         Author: Paul P. Hilscher (2012), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#include "config.h"


#ifndef FIELDS_PETSc_H
#define FIELDS_PETSc_H

#include "FieldsFFT.h"

#include "Global.h"
#include "Fields.h"


#include <petscksp.h>


static char help[] = "Help for PETSc Interface not available, please look up gkc & PETSc manual.";



class FieldsPETSc : public FieldsFFT {
  
  
  
  Vec            Vec_x,Vec_b;    /* approx solution, RHS, exact solution */
  Mat            A, F;        /* linear system matrix */
  KSP            ksp;      /* linear solver context */
  PC       pc;

  double tolerance;
  int order;

//    SuperMatrix    A, L, U;
//    SuperMatrix    vec_v, vec_x;

  void setupMatrix();
  Array3z rho_star;
protected:   
  Array3z virtual solvePoissonEquation(Array3z rho, Timing timing);
  Array4z virtual solveFieldEquations(Array4z Q, Timing timing) ;
    
private: 
  void setupMatrixPoissonEquation(int n_points=2);
    /**
     *
     *  calculates \f[ \rho_\star = \left( 1 +  \rho_{t\sigma} k_\perp^2 right) \rho \f]
     *
     * */
  void calcRhoStar(Array4z Q);

public:
 FieldsPETSc(Setup *setup, Grid *grid, Parallel *parallel, FileIO *fileIO, Geometry<HELIOS_GEOMETRY> *geo, FFTSolver *fftsolver);

 //* Destructor
 ~FieldsPETSc();
 
protected:
   virtual void printOn(ostream &output) const {
         output   << "Poisson   |  PETSc       Order " << order << " Tolerance : " << tolerance << std::endl;
   }

};


#endif // __FIELDS_PETSc_H



