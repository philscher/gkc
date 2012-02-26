/*
 * =====================================================================================
 *
 *       Filename: FieldsHyper.h
 *
 *    Description: Interface to Hypre Matrix Solver routine
 *
 *         Author: Paul P. Hilscher (2011), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#include "config.h"


#ifndef FIELDS_HYPRE_H
#define FIELDS_HYPRE_H

#include "FieldsFFT.h"

#include "Global.h"
#include "Fields.h"
#include "HYPRE_struct_ls.h"
#include "HYPRE_sstruct_ls.h"

/**
 *
 *
 * For the Poisson's equation, we need to solve the quasi-neuitrality
 * condition, which included $\f[ \Gamma_0 \f]$. We use  Pade` 
 * approximations, and split the calculation into 2 parts 
 * \f[ \Gamma_0 = \frac{b}{1+b} \f]
 *
 * See Gysela Proceesings 
 *
 * \f[ lambda_D^2 b + \sum_n \frac{\gamma_n}{b}{1+\gamma_n b} + \alpha \phi = \rho \f]
 *
 *  
 *
 *
 * */
class FieldsHypre : public FieldsFFT {

  double tolerance;

  HYPRE_StructGrid     poisson_grid;
  HYPRE_StructMatrix   A;
  HYPRE_StructVector   vec_b;
  HYPRE_StructVector   vec_x;
  HYPRE_StructSolver   solver;

  int NXky_LlD[2], NXky_LuD[2];
  int order;

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
 FieldsHypre(Setup *setup, Grid *grid, Parallel *parallel, FileIO *fileIO, Geometry<HELIOS_GEOMETRY> *geo, FFTSolver *fftsolver);

 //* Destructor
 ~FieldsHypre();
 
protected:
   virtual void printOn(ostream &output) const {
         output   << "Poisson   |  Hypre       Order " << order << " Tolerance : " << tolerance << std::endl;
   }

};


#endif // __FIELDS_HYPRE_H



