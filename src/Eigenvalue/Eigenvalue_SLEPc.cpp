/*
 * =====================================================================================
 *
 *       Filename: Eigenvalue_SLEPc.cpp
 *
 *    Description: Interface for SLEPc to calculate Eigenvalues
 *
 *         Author: Paul P. Hilscher (2011), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include "TimeIntegration.h"
#include "Eigenvalue_SLEPc.h"
#include "Grid.h"
#include "Timing.h"
#include "Tools/System.h"

#include "Special/PETScMatrixVector.h"

static char help[] = "SLEPc help and how to get rid of it ..... ?";


// We set our own PETSc error handler in order to have more control over
// tracebacks.
PetscErrorCode petsc_error_handler(MPI_Comm comm, int line, const char *fun, const char* file,
           const char *dir, PetscErrorCode n, PetscErrorType p, const char *mess, void *ctx)
{

  std::cout << "PETSc error detected " << std::endl 
            << "line     : "           << line << std::endl
            << "function : "           << std::string(fun)  << std::endl
            << "file     : "           << std::string(file) << std::endl
            << "Message  : "           << std::string(mess) << std::endl;
                
  System::printStackTrace();

  // PETSc errors are fatal
  abort();
  exit(0);
}

Eigenvalue_SLEPc::Eigenvalue_SLEPc(FileIO *fileIO, Setup *setup, Grid *grid, Parallel *_parallel) : Eigenvalue(fileIO, setup,grid,  _parallel)
{
  
  SlepcInitialize(&setup->argc, &setup->argv, (char *) 0,  help);
     
  // Don't let PETSc catch signals, we handle signals ourself in Control.cpp
  PetscPopSignalHandler();
  PetscPushErrorHandler(petsc_error_handler, NULL);

  // create Matrix Operations
  int subDiv = 0;

  // Note : We do not include Nyquist frequency ( and optionally Zonal Flow (ZF) per switch)
  int local_size  = NxLD * (includeZF ? Nky-1 : Nky-2) * NzLD * NvLD * NmLD * NsLD;
  int global_size = Nx   * (includeZF ? Nky-1 : Nky-2) * Nz   * Nv   * Nm   * Ns  ;

  MatCreateShell(parallel->Comm[DIR_ALL], local_size, local_size, global_size, global_size, &subDiv, &A_F1);
  MatSetFromOptions(A_F1);

  MatShellSetOperation(A_F1, MATOP_MULT, (void(*)()) PETScMatrixVector::MatrixVectorProduct);
      
  // Initialize eigenvalue solver
  EPSCreate(parallel->Comm[DIR_ALL], &EigvSolver);
  EPSSetProblemType(EigvSolver, EPS_NHEP);
  EPSSetOperators(EigvSolver, A_F1, PETSC_NULL);
    
  // Krylov-Schur most promising one, Power does not converge
  EPSSetType(EigvSolver, EPSKRYLOVSCHUR);
  
  // use e.g. -x "-st_shift 0.2" to set more properties
  EPSSetFromOptions(EigvSolver);

  initData(setup, fileIO);
}

Complex Eigenvalue_SLEPc::getMaxAbsEigenvalue(Vlasov *vlasov, Fields *fields) 
{
    
  EPSSetWhichEigenpairs(EigvSolver, EPS_LARGEST_MAGNITUDE);

  parallel->print("Using SLEPc to find maximum abolute eigenvalue ...");
  // Solve Eigenvalues
  PETScMatrixVector pMV(vlasov, fields);

  // For an estimate of time-step accuracy does not have to be large
  EPSSetTolerances(EigvSolver, 1.e-5, 10000); 
  //control->signalForceExit(true);
  EPSSolve(EigvSolver);
  //control->signalForceExit(false);

  // Get Eigenvalues
  int nconv = 0;
  EPSGetConverged(EigvSolver, &nconv);

  Vec    Vec_F1;   
   
  CComplex eigv, eigv_dummy;
  if(nconv > 0) MatGetVecs(A_F1, PETSC_NULL, &Vec_F1);
  EPSGetEigenpair(EigvSolver, 0, (PetscScalar *) &eigv, (PetscScalar *) &eigv_dummy, Vec_F1, Vec_F1);

  VecDestroy(&Vec_F1);

  std::stringstream msg; msg << "Maximum Absolute Eigenvalue : " << eigv << std::endl;
  parallel->print(msg.str());

  return eigv;
}

void Eigenvalue_SLEPc::solve(Vlasov *vlasov, Fields *fields, Visualization *visual, Control *control) 
{

  //EPSSetWhichEigenpairs(EigvSolver, EPS_LARGEST_REAL);
  //EPSSetWhichEigenpairs(EigvSolver, EPS_ALL);
    
  int n_eigv = Nx * (includeZF ? Nky-1 : Nky-2) * Nz * Nv * Nm * Ns;
  EPSSetDimensions(EigvSolver, n_eigv, PETSC_DECIDE, PETSC_DECIDE);

  // No need to set initial solution vector (found convergence does not improve) 
     
  //////// Solve ////////
  PETScMatrixVector pMV(vlasov, fields, includeZF);

  control->signalForceExit(true);
  EPSSolve(EigvSolver);
  control->signalForceExit(false);

  //EPSSetBalance(EigvSolver, EPS_BALANCE_ONESIDE, 16, 1.e-13);

  int nconv = 0;
  EPSGetConverged(EigvSolver, &nconv);

  EigenValue eigvTable;

  Vec Vec_F1, Vec_F1_dummy;   

  std::cout << "Number of Converged Eigenvalues : " << nconv << std::endl;
       
  CComplex  *x_F1; 
  if(nconv > 0) MatGetVecs(A_F1, PETSC_NULL, &Vec_F1);


  //////////////////////////    Read Out Results /////////////////////
  for(int nev = 0; nev < nconv; nev++) {
          
    Complex eigv, eigv_dummy;
    EPSGetEigenpair(EigvSolver, nev, (PetscScalar *) &eigv, (PetscScalar *) &eigv_dummy, Vec_F1, Vec_F1_dummy);

    // Get numerical error of eigensolution
    eigvTable.EigenValue     = eigv;
    EPSComputeRelativeError(EigvSolver, nev, &eigvTable.AbsoluteError);
            
    EVTable->append(&eigvTable);
          
    // Skip eigenvalue results if growth rates are smaller than minimum value
    if(real(eigv) > growth_min) { // only write eigenvalues with larger growth rates otherwise ignore
    
      // Get EigenVector (Phase Space function) and calculate corresponding potentials
      std::cout << "Saving eigenvector : " << eigv << std::endl;
         
      VecGetArray(Vec_F1, (PetscScalar **) &x_F1);
   
      
      // Copy back PETSc solution vector to array
      [=](CComplex fs[NsLD][NmLD][NzLB][Nky][NxLB][NvLB]) {

      
        // copy whole phase space function (important due to boundary conditions)
        int n = 0;
        for(int s = NsLlD; s <= NsLuD; s++) { for(int m = NmLlD; m <= NmLuD ; m++  ) { 
        for(int z = NzLlD; z <= NzLuD; z++) {
        for(int y_k = (includeZF ? 0 : 1); y_k <= NkyLuD-1; y_k++) { 
        for(int x = NxLlD; x <= NxLuD; x++) { for(int v = NvLlD; v <= NvLuD; v++) { 

          fs[s][m][z][y_k][x][v] = x_F1[n++];
        
        }}} }}}
      
      }((A6zz) vlasov->fs);
           
      fields->solve(vlasov->f0,  vlasov->fs); 
   
      VecRestoreArray    (Vec_F1, (PetscScalar **) &x_F1);
        
      // Write put Fields Eigenvector, set Imaginary part as Time 
      visual->writeData(Timing(nev, real(eigv)), 0., true);
    }
  }
}

Eigenvalue_SLEPc::~Eigenvalue_SLEPc() 
{
  MatDestroy(&A_F1);
  EPSDestroy(&EigvSolver);

  // check error code of SLEPc finalize
  SlepcFinalize();  
  H5Gclose(eigvGroupID);
  delete EVTable;
}

void Eigenvalue_SLEPc::printOn(std::ostream &output) const 
{
  int global_size = Nx * (includeZF ? Nky-1 : Nky-2) * Nz * Nv * Nm * Ns;
  output << "Eigenvalue |  using SLEPc interface " << std::endl;
  output << "           | Approx ~ " << global_size << " Iterations " << std::endl;
}

void Eigenvalue_SLEPc::initData(Setup *setup, FileIO *fileIO) 
{
  eigvGroupID = fileIO->newGroup("Eigenvalue");

  /////////// Table to store eigenvalues (@todo move to Eigenvalues.cpp) ///
  EigenValue EigVal_table;
         
#pragma warning (disable : 1875) // ignore warnings about non-POD types
  size_t EigVal_offsets[] = { HOFFSET(EigenValue, EigenValue), HOFFSET(EigenValue, AbsoluteError) };
#pragma warning (enable  : 1875)
  size_t EigVal_sizes  [] = { sizeof(EigVal_table.EigenValue) ,  sizeof(EigVal_table.AbsoluteError) };
  hid_t  EigVal_types  [] = { fileIO->complex_tid , H5T_NATIVE_DOUBLE };
  const char * EigVal_names  [] = {"Eigenvalue", "Absolute Error"};

  EVTable = new TableAttr(eigvGroupID, "EigenValues", 2, EigVal_names, EigVal_offsets, EigVal_types, EigVal_sizes, &EigVal_table);
}

