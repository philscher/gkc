/*
 * =====================================================================================
 *
 *       Filename: TimeIntegration.h
 *
 *    Description:Time Integration Interface for gkc. 
 *
 *         Author: Paul P. Hilscher (2011), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */



#include "config.h"
#include "Global.h"

#ifndef SCANPOLOIDALEIGEN_H__
#define SCANPOLOIDALEIGEN_H__

#include "Setup.h"
#include "Grid.h"
#include "Parallel/Parallel.h"
#include "TimeIntegration.h"
#include "Eigenvalue/Eigenvalue_SLEPc.h"
#include "Fields.h"
#include "Init.h"

#include "Matrix/MatrixSolver.h"
#include "Matrix/MatrixPETSc.h"

#include "external/Array.h"
#include "Visualization/Visualization.h"
#include "TimeIntegration/TimeIntegration.h"
#include "Control.h"
#include "Analysis/Diagnostics.h"

#pragma warning(disable : 1478)
//#include "elemental.hpp"

//typedef elem::Complex<double> C;


//using namespace elem;


/**
*  @brief Class which handles explicit time integration
*
*  This classed is used to investigate the damping through
*  linear modes. First the eigenvalue/eigenvector pairs
*  are extracted using SLEPc. Then, we evolve the gk-system
*  as an IVP, where every couple of time steps we decompose
*  the solution into their linear eigenmodes, to see the
*  perturbation of the sub-dominant modes - this is specially
*  interesting in the non-linear regime.
*
*  References :
*
*    Hatch et al.
*
*
**/

/* 
class ScanPoloidalEigen  : public IfaceGKC {
 protected:
   
   Parallel *parallel;
   Grid     *grid;
   FileIO   *fileIO;
   Fields   *fields;
   Vlasov   *vlasov;
   Control *control;
   Visualization *visual;
   Diagnostics *diagnostics;
   TimeIntegration *timeIntegration;
  
   Array1<MatrixSolver*> EigVecDecompSolver;


   Vec Vec_EigVec, 
       Vec_Pertb;

   CComplex *Eigvec_F1, 
            *EigSystem_F1;

   nct::allocate ArrayEigvec, ArrayEigSystem;


   int y_k_; ///< Poloidal mode to investigate

  hid_t     scanGroupID;  ///< HDF-5 reference to Eigenvalue group
  
  FileAttr *FA_eigvec, *FA_eigval,
    *FA_eigerr,
    *FA_pertb,
    *FA_EigTime;

  bool elemental_init; ///< Set if Elemental i initialiazed

  
  typedef elem::DistMatrix<C> Matrix;

  std::vector<Matrix> Matrix_ky; 

  Matrix Mat_Eigvec_ky,
         Eigvec_x, 
         Eigvec_y;

  void solveIVP(Eigenvalue *eigenvalue, int y_k_);
  void setupEigenDecompSolver(int y_k_);
    
  int nStepDecomp; /// Decompose every nStepDecomp

 public:
   //   Constructor
   ScanPoloidalEigen(Setup *setup, Grid *_grid, Parallel *_parallel, Vlasov *_vlasov, Fields *_fields, FileIO *fileIO,
    Control *_control, Visualization *_visual, Diagnostics *_diagnostics, TimeIntegration *_timeIntegration);
  ~ScanPoloidalEigen();
  
   void solveEigenSystem(Vlasov *vlasov, Fields *fields, Eigenvalue_SLEPc *eigen, int y_k_) ;
   void solve(Vlasov *vlasov, Fields *fields, TimeIntegration *timeIntegration, Eigenvalue *eigenvalue, Init *init, Visualization *visual);
   void getPerturbations(int y_k, Timing timing);
   
   virtual void printOn(std::ostream &output) const {} ;
  
   void initData(Setup *setup, FileIO *fileIO);
};

 * */
#endif // SCANPOLOIDALEIGEN_H__
