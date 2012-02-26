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

#include "Special/PETScMatrixVector.h"

static char help[] = "SLEPc help and how to get rid of it ..... ?";


Eigenvalue_SLEPc::Eigenvalue_SLEPc(Setup *setup, Grid *grid, Parallel *_parallel) : Eigenvalue(setup,grid,  _parallel)
{
    
      SlepcInitialize(&setup->argc, &setup->argv, (char *) 0,  help);
      
      // create Matrix Operations
      int subDiv=0;
      MatCreateShell(parallel->Comm[DIR_ALL], grid->getLocalSize(), grid->getLocalSize(), grid->getGlobalSize(), grid->getGlobalSize(), &subDiv, &A_F1);
      MatSetFromOptions(A_F1);

      MatShellSetOperation(A_F1, MATOP_MULT, (void(*)()) PETScMatrixVector::MatrixVectorProduct);
      
      // Initialize eigenvalue solver
      EPSCreate(parallel->Comm[DIR_ALL], &EigvSolver);
      EPSSetProblemType(EigvSolver, EPS_NHEP);
      EPSSetOperators(EigvSolver, A_F1, PETSC_NULL);
      
      // use e.g. -x "-st_shift 0.2" to set more properties
      EPSSetFromOptions(EigvSolver);

}



cmplxd Eigenvalue_SLEPc::getMaxAbsEigenvalue(Vlasov *vlasov, Fields *fields) 
{
    
    EPSSetWhichEigenpairs(EigvSolver, EPS_LARGEST_MAGNITUDE);

    // Krylov-Schur most promising one, Power does not converge
    EPSSetType(EigvSolver, EPSKRYLOVSCHUR);

    // Solve Eigenvalues
    PETScMatrixVector pMV(vlasov, fields);
    EPSSolve(EigvSolver);


    // Get Eigenvalues
    int nconv = 0;
    EPSGetConverged(EigvSolver, &nconv);

    Vec    Vec_F1;   
    cmplxd *x_Vec_F1 = PETScMatrixVector::getCreateVector(grid, Vec_F1);
   
    cmplxd eigv, eigv_dummy;
    if(nconv > 0) MatGetVecs(A_F1, PETSC_NULL, &Vec_F1);
    EPSGetEigenpair(EigvSolver, 0, &eigv, &eigv_dummy, Vec_F1, Vec_F1);

    VecDestroy(&Vec_F1);

    std::stringstream msg; msg << " Eigenvalue : " << eigv << std::endl;
    parallel->print(msg);

    return eigv;

}


void Eigenvalue_SLEPc::solve(Vlasov *vlasov, Fields *fields, Visualization *visual, Control *control) 
{


    EPSSetWhichEigenpairs(EigvSolver, EPS_LARGEST_REAL);

    // init intial solution vector 
    Vec Vec_init;
    cmplxd *init_x = PETScMatrixVector::getCreateVector(grid, Vec_init);
    
   for(int x = NxLlD, n = 0; x <= NxLuD; x++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { for(int z = NzLlD; z <= NzLuD; z++) {
   for(int v = NvLlD       ; v <= NvLuD; v++) { for(int m   = NmLlD ; m   <= NmLuD ; m++  ) { for(int s = NsLlD; s <= NsLuD; s++) {

                init_x[n++] = 1.e-5 * vlasov->f0(x,y_k,z,v,m,s);

     }}} }}}

     VecRestoreArray    (Vec_init, &init_x);
     EPSSetInitialSpace(EigvSolver, 0, &Vec_init);
    //////// Solve ////////
    EPSSolve(EigvSolver);
     VecDestroy(&Vec_init);
    int nconv = 0;
    EPSGetConverged(EigvSolver, &nconv);

    EigenValue eigvTable[nconv];

    // get Eigensolver type
    const EPSType type;
    EPSGetType(EigvSolver, &type);

    EPSPrintSolution(EigvSolver, PETSC_NULL);

    Vec    Vec_F1, Vec_F1_dummy;   

    std::cout << "Number of Converged Eigenvalues : " << nconv << std::endl;
	    
    cmplxd  *x_F1; 
    if(nconv > 0) MatGetVecs(A_F1, PETSC_NULL, &Vec_F1);


    //////////////////////////    Read Out Results /////////////////////
    for(int m = 0; m < nconv; m++) {
          
            cmplxd eigv, eigv_dummy;
            EPSGetEigenpair(EigvSolver, m, &eigv, &eigv_dummy, Vec_F1, Vec_F1_dummy);
            std::cout << "Eigenvalue : " << eigv << " " << std::endl;

            double error = 0.;

            eigvTable[m].EigenValue     = eigv;
            EPSComputeRelativeError(EigvSolver,m, &error);
            eigvTable[m].AbsoluteError  = error;
            
            
            // Get EigenVector (Phase Space function) and calculate corresponding potentials
  	    VecGetArray(Vec_F1, &x_F1);

            // copy whole phase space function (waste but starting point) (important due to bounday conditions
           // we can built wrapper around this and directly pass it
   for(int x = NxLlD, n = 0; x <= NxLuD; x++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { for(int z = NzLlD; z <= NzLuD; z++) {
   for(int v = NvLlD       ; v <= NvLuD; v++) { for(int m   = NmLlD ; m   <= NmLuD ; m++  ) { for(int s = NsLlD; s <= NsLuD; s++) {

                vlasov->fs(x,y_k,z,v,m,s) = x_F1[n++];

           }}} }}}
           
           fields->solve(vlasov->f0,  vlasov->fs); 
   
           VecRestoreArray    (Vec_F1, &x_F1);
        
           visual->writeData(Timing(m,0.), 0., true);
  //  	   EigVal_TableAttr->append(&eigvTable[m], 1); 
//    void *ctx; MatShellGetContext(A, &ctx);
//    int N =  *(int *) ctx;
    }


};



  
  
  
Eigenvalue_SLEPc::~Eigenvalue_SLEPc() {
 	
      MatDestroy(&A_F1);
      EPSDestroy(&EigvSolver);

      // check error code of slepc finalize
      SlepcFinalize();  
    delete EigVal_TableAttr;
};

void Eigenvalue_SLEPc::print2On(ostream &output) {
           output << "Eigenvalue |  using SLEPc interface " << std::endl;
 }


 /////////////////////////////////// Data I/O Stuff ////////////////




void Eigenvalue_SLEPc::initDataOutput(Setup *setup, FileIO *fileIO) 
{
//	auto groupID = fileIO->newGroup(bla, fileIO);
	hid_t eigvGroupID = fileIO->newGroup(fileIO->getFileID(), "Eigenvalue");

    // ********************* setup Table for EigenValues *****************
    EigenValue EigVal_table;
         
    size_t EigVal_offsets[] = { HOFFSET( EigenValue, EigenValue ), HOFFSET( EigenValue, AbsoluteError) };
    size_t EigVal_sizes  [] = { sizeof(EigVal_table.EigenValue) ,  sizeof(EigVal_table.AbsoluteError) };
    hid_t  EigVal_types  [] = { fileIO->complex_tid , H5T_NATIVE_DOUBLE };
    const char * EigVal_names  [] = {"Eigenvalue", "Absolute Error"};

    EigVal_TableAttr = new TableAttr(eigvGroupID, "EigenValues", 2, EigVal_names, EigVal_offsets, EigVal_types, EigVal_sizes, &EigVal_table);


    H5Gclose(eigvGroupID);



}
