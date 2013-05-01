/*
 * =====================================================================================
 *
 *       Filename: ScanPoloidalEigen.h
 *
 *    Description:Time Integration Interface for gkc. 
 *
 *         Author: Paul P. Hilscher (2011), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include "ScanPoloidalEigen.h"
#include "Special/PETScMatrixVector.h"
#include "Tools/Tools.h"

/**
*   @brief Initialize Elemental
*
*   We need to initialize elemental as Class includes variables with standard constructor
*   which implicitly calls elem::DefaultGrid (which has to be initialized first)
*
*
**/ 
/* 
bool init_Elemental(Setup *setup) { 
  
  elem::Initialize(setup->argc, setup->argv);
  return true; 
}

int numEigv; /// total numbers of eigenvalues for each poloidal modes


ScanPoloidalEigen::ScanPoloidalEigen(Setup *setup, Grid *_grid, Parallel *_parallel, Vlasov *_vlasov, Fields *_fields, FileIO *_fileIO,
    Control *_control, Visualization *_visual, Diagnostics *_diagnostics, TimeIntegration *_timeIntegration)
     :
  elemental_init(init_Elemental(setup)), 

 parallel(_parallel), grid(_grid), vlasov(_vlasov), fields(_fields), fileIO(_fileIO),
 control(_control), visual(_visual), diagnostics(_diagnostics), timeIntegration(_timeIntegration)
{

  numEigv = Nx * Nz * Nv * Nm * Ns;

  ArrayEigSystem = nct::allocate(nct::Range(0, numEigv), grid->RxLD, grid->RvLD)(&EigSystem_F1);
  ArrayEigvec    = nct::allocate(grid->RxLD, grid->RvLD)(&Eigvec_F1);
  
  nStepDecomp = setup->get("EigenSystem.nStepDecomp", 200);
  y_k_        = setup->get("EigenSystem.PoloidalMode", 3 );
  

  initData(setup, fileIO); 

  ///////////////////////  Setup elemental package ////////////////////////////////
  
  const elem::Grid g(parallel->Comm[DIR_X]); // DIR_XV ?

  // Setup matrix and vectors
  Mat_Eigvec_ky = elem::DistMatrix<C>(numEigv, numEigv, g);
  Eigvec_x      = elem::DistMatrix<C>(numEigv, 1      , g);
  Eigvec_y      = elem::DistMatrix<C>(numEigv, 1      , g);
 

  parallel->print("Using y_k_ : " + Setup::num2str(y_k_) + "\n");
}

ScanPoloidalEigen::~ScanPoloidalEigen()
{
  delete FA_eigvec;
  delete FA_eigval;
  delete FA_eigerr;

  delete FA_pertb;
  delete FA_EigTime;
    
  H5Gclose(scanGroupID);
  
  elem::Finalize();
}



void ScanPoloidalEigen::setupEigenDecompSolver(int y_k_)
{

  parallel->print("Filling matrix");
 
  // fill matrix
  [&](CComplex EigSystem_F1[numEigv][NxLD][NvLD]) {

    for(int nev = 0; nev < numEigv; nev++) {
     
    for(int x = NxLlD, n = 0; x <= NxLuD; x++) { for(int v = NvLlD; v <= NvLuD; v++, n++) {

      Mat_Eigvec_ky.SetLocal(n, nev, C(EigSystem_F1[nev][x][v]));
    } } }
  
  } ((A3zz) EigSystem_F1);
 
  elem::Inverse(Mat_Eigvec_ky);
  
  parallel->print("Done solving LU decomposition");

}

void ScanPoloidalEigen::getPerturbations(int y_k_, Timing timing)
{
  
  ///////////  Set current value of f1 ////////////////////
  [=](CComplex f[NsLD][NmLD][NzLB][Nky][NxLB][NvLB]) {

    for(int x = NxLlD, n = 0; x <= NxLuD; x++) { for(int v = NvLlD; v <= NvLuD; v++, n++) {

      Eigvec_y.SetLocal(n, 0, C(f[NsLlD][NmLlD][NzLlD][y_k_][x][v]));
    } }

  } ((A6zz) vlasov->f);

  ///////////// Solve using inverse matrix x = A^-1 y /////////////////////
  elem::Gemv(elem::NORMAL, C(1.), Mat_Eigvec_ky, Eigvec_y, C(0.), Eigvec_x);
  
  CComplex Pertb[numEigv]; Pertb[:] = 0.;
 
  /////  Get perturbations of eigenvectors /////////////////
  int colShift  = Eigvec_x.ColShift();
  int colStride = Eigvec_x.ColStride();
  int colLength = Eigvec_x.LocalHeight();
  
  for( int n = 0; n < colLength; ++n ) {

    int idx = colShift + n * colStride;
    // WARNING ! breaks strict aliasing rules ?! (should use temporary ... ?)
    Pertb[idx] = *((CComplex *)  &Eigvec_x.GetLocal(n, 0));
  } 

  parallel->reduce(Pertb, Op::sum, DIR_ALL, numEigv);

  FA_pertb->write(Pertb);
  FA_EigTime->write(&timing);
  
  return;
}           
         

void ScanPoloidalEigen::solve(Vlasov *vlasov, Fields *fields, TimeIntegration *timeIntegration, Eigenvalue *eigenvalue, Init *init, Visualization *visual)
{
 // int ky[] = { 2,3,4 };
  //Matrix Mat_Eigvec_ky = elem::DistMatrix<C>(numEigv, numEigv, g);
   
  // Solve complete eigensystem
 // for(int y_k_ = 1; y_k_ < Nky; y_k_++) {
//    if(y_k_ != 1) continue;
  
    parallel->print("Solving eigensystem");
    solveEigenSystem(vlasov, fields, (Eigenvalue_SLEPc *) eigenvalue, y_k_);
    
    parallel->print("LU decomposing eigensystem");
    setupEigenDecompSolver(y_k_);
    
    parallel->print("Solving IVP");
    solveIVP(eigenvalue, y_k_);
  
//  }
}

void ScanPoloidalEigen::solveEigenSystem(Vlasov *vlasov, Fields *fields, Eigenvalue_SLEPc *eigen, int y_k_) 
{
  // PETSC_DECIDE for ncv does not always converged from decomposed problem
  EPSSetDimensions(eigen->EigvSolver, numEigv, numEigv, PETSC_DECIDE);
  EPSSetTolerances(eigen->EigvSolver, 1.e-12, 100000); 
  //EPSSetBalance(eigen->EigvSolver, EPS_BALANCE_ONESIDE, 4, 1.e-9);

  //////// Solve //////// (simpler is to pass directly vector<ky_modes>
  PETScMatrixVector pMV(vlasov, fields, false, y_k_);

  ////////////////////////////////// solve eigensystem
  EPSSolve(eigen->EigvSolver);
  int nconv = 0; EPSGetConverged(eigen->EigvSolver, &nconv);

  check(nconv == numEigv ? 1 : -1, DMESG("Eigensystem did not converged"));

  std::cout << "Converged : " << nconv << std::endl;
  
  //////////////////////////    Read Out Results /////////////////////
  Vec Vec_F1, Vec_F1_dummy;   
  MatGetVecs(eigen->A_F1, PETSC_NULL, &Vec_F1);

  for(int nev = 0; nev < numEigv; nev++) {
    
    Complex eigv, eigv_dummy;

   EPSGetEigenpair(eigen->EigvSolver, nev, (PetscScalar *) &eigv, (PetscScalar *) &eigv_dummy, Vec_F1, Vec_F1_dummy);
    double eigerr;
    EPSComputeRelativeError(eigen->EigvSolver, nev, &eigerr);
    
    CComplex  *x_F1; 
    VecGetArray(Vec_F1, (PetscScalar **) &x_F1);
    
    // Copy eigenvector to eigensystem matrix
    [=](CComplex Eigvec_F1[NxLD][NvLD], CComplex EigSystem_F1[numEigv][NxLD][NvLD]) {

      for(int x = NxLlD, n = 0; x <= NxLuD; x++) { for(int v = NvLlD; v <= NvLuD; v++, n++) { 

        Eigvec_F1        [x][v] = x_F1[n];
      
      } }
        
      // Eigenvectors can differ by complex factor a between each other (e.g. decomposition), 
      // thus we need to normalize 
      const double max_val = parallel->reduce(__sec_reduce_max(cabs(Eigvec_F1[NxLlD:NxLD][NvLlD:NvLD])), Op::max); 
      const double phase   = carg(parallel->reduce(__sec_reduce_add(Eigvec_F1[NxLlD:NxLD][NvLlD:NvLD]), Op::sum));
      
      // normalize to max(|phi|) = 1 and set total phase to 0
      Eigvec_F1        [NxLlD:NxLD][NvLlD:NvLD] *= cexp(- _imag * phase) / max_val ;
      
      EigSystem_F1[nev][NxLlD:NxLD][NvLlD:NvLD] = Eigvec_F1[NxLlD:NxLD][NvLlD:NvLD];

    } ((A2zz) Eigvec_F1, (A3zz) EigSystem_F1);
   
    // save eigenvector / eigenvalue / eigenvalue (error)
    FA_eigvec->write(ArrayEigvec.data(Eigvec_F1));
    FA_eigval->write(&eigv);
    FA_eigerr->write(&eigerr);

    VecRestoreArray(Vec_F1, (PetscScalar **) &x_F1);
  }
        
}

void ScanPoloidalEigen::initData(Setup *setup, FileIO *fileIO) 
{
  scanGroupID = fileIO->newGroup("EigenSystem");
  
  // Set 
 check(H5LTset_attribute_int(scanGroupID, ".", "nStepDecomp" , &nStepDecomp, 1), DMESG("Attribute"));
 check(H5LTset_attribute_int(scanGroupID, ".", "PoloidalMode", &y_k_       , 1), DMESG("Attribute"));
  
  // Phase space dimensions
  {
    hsize_t dim[]       = { Nx     , Nv     ,             1 };
    hsize_t maxdim[]    = { Nx     , Nv     , H5S_UNLIMITED };
    hsize_t chunkBdim[] = { NxLD   , NvLD   ,             1 };
    hsize_t chunkdim[]  = { NxLD   , NvLD   ,             1 };
    hsize_t offset[]    = { NxLlB-1, NvLlB-1,             0 };
//    hsize_t moffset[]   = { NxLlB-1, NvLlB-1,             0 };
    //hsize_t moffset[]   = {      2 , 2      ,             0 };
    hsize_t moffset[]   = {      0 , 0      ,             0 };
    
    
 //     hsize_t dim[]       = { Nky, Nx     , Nv     ,             1 };
//    hsize_t maxdim[]    = { Nky, Nx     , Nv     , H5S_UNLIMITED };
//    hsize_t chunkBdim[] = { Nky, NxLB   , NvLB   ,             1 };
//    hsize_t chunkdim[]  = { Nky, NxLD   , NvLD   ,             1 };
//    hsize_t offset[]    = { Nky, NxLlB-1, NvLlB-1,             0 };
//    hsize_t moffset[]   = { Nky,       0,       0,             0 };
    FA_eigvec  = new FileAttr("eigvec", scanGroupID, fileIO->file, 3, dim, maxdim, chunkdim, moffset,  chunkBdim, offset, true, fileIO->complex_tid);
  }
  
  bool writeEigval = parallel->myRank == 0; // only master writes

  // Eigenvalues
  {
    hsize_t dim[]    = {  1 };
    hsize_t maxdim[] = { H5S_UNLIMITED };
    hsize_t off0[]   = { 0 };
 

    FA_eigval  = new FileAttr("eigval", scanGroupID, fileIO->file, 1, dim, maxdim, dim, off0, dim, off0, writeEigval, fileIO->complex_tid);
    FA_eigerr  = new FileAttr("eigerr", scanGroupID, fileIO->file, 1, dim, maxdim, dim, off0, dim, off0, writeEigval);
    
  }

  {
    hsize_t dim[]    = { numEigv,             1 };
    hsize_t maxdim[] = { numEigv, H5S_UNLIMITED };
    hsize_t off0[]   = { 0      ,             0 };
  
    FA_pertb   = new FileAttr("Pertb", scanGroupID, fileIO->file, 2, dim, maxdim, dim, off0, dim, off0, writeEigval, fileIO->complex_tid);
    FA_EigTime = fileIO->newTiming(scanGroupID);

  }
  
}

void ScanPoloidalEigen::solveIVP(Eigenvalue *eigenvalue, int y_k_)
{
  
  ////////////////// Start IVP Loop /////////////////////

  Timing timing(0,0.) ;
      
  timeIntegration->setMaxLinearTimeStep(eigenvalue, vlasov, fields);

  bool isOK = true; 
         
  do {
   
    /////////////////// Standard IVP Loop ////////////////////
        
    
    // integrate for one time-step, give current dt as output
    const double dt = timeIntegration->solveTimeStep(vlasov, fields, NULL, timing);     

    vlasov->writeData(timing, dt);
    fields->writeData(timing, dt);
    visual->writeData(timing, dt);

    diagnostics->writeData(timing, dt);
    
    // flush in regular intervals in order to minimize HDF-5 file 
    // corruption in case of an abnormal program termination
    fileIO->flush(timing, dt);  
    
    isOK = control->checkOK(timing, timeIntegration->maxTiming);
    
    if(timing.step % nStepDecomp == 0) getPerturbations(y_k_, timing);
    
  } while(isOK);
}

 * */
