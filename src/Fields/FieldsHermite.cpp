/*
 * =====================================================================================
 *
 *       Filename: FieldsHermite.cpp
 *
 *    Description: Solver Field equations using Hermite interpolation (x, ky) to
 *                 calculate gyro-averaging.
 *
 *         Author: Paul P. Hilscher (2012), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#include "Fields/FieldsHermite.h"

#include "Special/HermiteInterpolation.h"
#include "Special/Integrate/Integrate.h"
#include "Matrix/SpecialMatrix.h"

#include "Tools.h"

#include "System.h"

#include "petscmat.h" 

static char FieldsHermite_help[] = "Help for PETSc Interface not available, please look up gkc & PETSc manual.";



FieldsHermite::FieldsHermite(Setup *setup, Grid *grid, Parallel *parallel, FileIO *fileIO, Geometry *geo) 
: Fields(setup, grid, parallel,fileIO, geo)
{
      
  // Move to Master PETSc class (create one ...) 
  
  PetscInitialize(&setup->argc, &setup->argv, (char *) 0,  FieldsHermite_help);
  
  //PetscPopSignalHandler();
  
  GyroMatrix = Array4<Matrix *>(grid->RkyLD, grid->RzLD, grid->RmLD, grid->RsLD); 
  MatrixPoissonSolverLHS = Array2<MatrixSolver *>(grid->RkyLD, grid->RzLD);
       
  interpolationOrder = setup->get("FieldsHermite.InterpolationOrder", 5);

  // Create gyro-averaging matrix Matrix
  for(int s = NsLlD; s <= NsLuD; s++) {  for(int   m = NmLlD ; m   <= NmLuD ; m++  ) {
  for(int z = NzLlD; z <= NzLuD; z++) {  for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { 
        
    GyroMatrix(y_k, z, m, s) =  getGyroAveragingMatrix(M[m], y_k, z, s);
          
  } } } }

    
  setDoubleGyroAverageMatrix(setup);
        
  // Create Vectors
  VecCreateMPI(parallel->Comm[DIR_X], NxLD, grid->NxGD, &GyroVector_X);
  VecCreateMPI(parallel->Comm[DIR_X], NxLD, grid->NxGD, &GyroVector_XAvrg);
  VecAssemblyBegin(GyroVector_X    );    VecAssemblyEnd  (GyroVector_X    );
  VecAssemblyBegin(GyroVector_XAvrg);    VecAssemblyEnd  (GyroVector_XAvrg);

}


void FieldsHermite::gyroAverage(const CComplex In [Nq][NzLD][NkyLD][NxLD], 
                                      CComplex Out[Nq][NzLD][NkyLD][NxLD],
                                const int m, const int s, const bool forward, const bool stack)
{
  
  CComplex *vec_X, *vec_XAvrg;
  
  // Create Matrix
  // Is it possible to fill array at once ?! otherwise parallelization efficiency will be bad
  for(int q = 0; q < Nq; q++) {
        
  for(int  z = NzLlD; z <= NzLuD; z++) {  for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { 
        
    VecGetArray    (GyroVector_X, (PetscScalar **) &vec_X);
    for(int x = NxLlD, n = 0; x <= NxLuD; x++)  vec_X[n++] = In[q][z][y_k][x];
    VecRestoreArray(GyroVector_X, (PetscScalar **) &vec_X);
       
    MatMult(GyroMatrix(y_k,z,m,s)->getMat(), GyroVector_X, GyroVector_XAvrg);
           
    VecGetArrayRead    (GyroVector_XAvrg, (const PetscScalar **) &vec_XAvrg);
    for(int x=NxLlD, n = 0; x <= NxLuD; x++)  Out[q][z][y_k][x] = vec_XAvrg[n++]; 
    VecRestoreArrayRead(GyroVector_XAvrg, (const PetscScalar **) &vec_XAvrg);
        
    
  } } } 

  return;

};


CComplex FieldsHermite::getElements(const int x, const int n, const double r, int y_k, const int z) 
{
  
  const double ky = (2. * M_PI/Ly) * y_k;
  // Assign the lambda expression that adds two numbers to an auto variable.
  auto Integrand = [=] (double alpha) -> Complex { 
           
  
    /////////////////////////////////////////////////////////////////////////////////////
    // Linearize the metric, take care of the errors (Lapillone, 2009, PhD, Thesis).
    //  (Do we use covariant or contra-variant metric ?!)
    //
    //  g^2     = g^{xx} g^{yy} - (g_{xy})^2
    // 
    //  rho_x   = sqrt{g_{xx}} \rho \cos{(\alpha)}
    //  rho_y   = frac{g_{xy}}{\sqrt{g_{xx}}} rho \cos{(\alpha}} 
    //              + \frac{g}{\sqrt{g_{xx}} \rho sin{(\alpha)}
    //
    const double g = sqrt( geo->g_xx(x,z) * geo->g_yy(x,z) - pow2(geo->g_xy(x,z)) );
  
    const double r_x = sqrt(geo->g_xx(x, z)) * r * cos(alpha);

    const double r_y = (geo->g_xy(x, z) / geo->g_xx(x,z)) * r * cos(alpha) +
                               g / geo->g_xx(x,z) * r * sin(alpha);
           
    // Is this correct ? (Matrix is still not singular ?)
    //     if( ( i <= (NxGlD+1)) || ( i  >= (NxGuD-1))) return 0.;
    //     if( ( n <= (NxGlD+1)) || ( n  >= (NxGuD-1))) return 0.;
    if( ((X[x] - r_x) <= X[NxGlD]) || ((X[x] - r_x) >= X[NxGuD])) return 0.;

    return Lambda(X[x] - r_x, n) * exp(Complex(0.,1.) * ky * r_y);

    
  };
        
  
  // Integration is very sensitif on order ( why ? some bug ?? )
  check(-1, DMESG("Fix _Complex double -> std::complex<double> cast"));
  
  return 1./(2.*M_PI) * Integrate::GaussLegendre(Integrand, 0., 2.*M_PI, 128);
        

};

FieldsHermite::~FieldsHermite()
{
  // shutdown PETSc !?
}



//Array4C FieldsHermite::solveFieldEquations(Array4C Qn, Timing timing) 
void FieldsHermite::solveFieldEquations(const CComplex Q     [Nq][NxLD][NkyLD][Nz],
                                              CComplex Field0[Nq][NxLD][NkyLD][Nz])
{
   
  // for 3 fields phi and B_perp are coupled and need to be solved differently
  if( solveEq & Field::Iphi) solvePoissonEquation  (Q, Field0);
  if( solveEq & Field::IAp ) check(-1, DMESG("Not Implemented"));
  if( solveEq & Field::IBp ) check(-1, DMESG("Not Implemented"));
    
  return;

}


// Take care of X[n] overflow
double FieldsHermite::Lambda(const double x, const int n) 
{
        
  // Scale input to [ 0 , 1 ]
  auto N = [=] (const double x, const int n) -> double { return (x-X[n])/dx ; };
        
  
  // Hermite Interpolation has to base functions corresponding to 1 and 0 and their derivatives 
  // We use Central Difference n derivative m order stencil (CD_n-m)

  // 1st order interpolation
  if      (interpolationOrder == 1) {
            check(-1, DMESG("Interpolation order not implemented"));
  }
  // 3rd order interpolation
  else if(interpolationOrder == 3) {

    const double H3_0_n   = HermiteInterpolation::H3_00(N(x,n))   + HermiteInterpolation::H3_10(N(x,n));

    const double H3_1_np1 = HermiteInterpolation::H3_01(N(x,n+1)) + HermiteInterpolation::H3_11(N(x,n+1));
    const double H3_1_nm1 = HermiteInterpolation::H3_01(N(x,n-1)) + HermiteInterpolation::H3_11(N(x,n-1));
    const double H3_1_np2 = HermiteInterpolation::H3_01(N(x,n+2)) + HermiteInterpolation::H3_11(N(x,n+2));
    const double H3_1_nm2 = HermiteInterpolation::H3_01(N(x,n-2)) + HermiteInterpolation::H3_11(N(x,n-2));

    
    // D_0 H3_0 + D_1 H3_1
    
    return H3_0_n  
           + ( (H3_1_nm2 - H3_1_np2)  -  8. * ( H3_1_nm1 - H3_1_np1 )) / ( 12.* dx);
  
  } 
  // 5th order interpolation
  else if(interpolationOrder == 5) {

  
    const double H5_0_n   = HermiteInterpolation::H5_00(N(x,n))   + HermiteInterpolation::H5_10(N(x,n));

    const double H5_1_np1 = HermiteInterpolation::H5_01(N(x,n+1)) + HermiteInterpolation::H5_11(N(x,n+1));
    const double H5_1_nm1 = HermiteInterpolation::H5_01(N(x,n-1)) + HermiteInterpolation::H5_11(N(x,n-1));
    const double H5_1_np2 = HermiteInterpolation::H5_01(N(x,n+2)) + HermiteInterpolation::H5_11(N(x,n+2));
    const double H5_1_nm2 = HermiteInterpolation::H5_01(N(x,n-2)) + HermiteInterpolation::H5_11(N(x,n-2));

    const double H5_2_n   = HermiteInterpolation::H5_02(N(x,n  )) + HermiteInterpolation::H5_12(N(x,n  ));
    const double H5_2_np1 = HermiteInterpolation::H5_02(N(x,n+1)) + HermiteInterpolation::H5_12(N(x,n+1));
    const double H5_2_nm1 = HermiteInterpolation::H5_02(N(x,n-1)) + HermiteInterpolation::H5_12(N(x,n-1));
    const double H5_2_np2 = HermiteInterpolation::H5_02(N(x,n+2)) + HermiteInterpolation::H5_12(N(x,n+2));
    const double H5_2_nm2 = HermiteInterpolation::H5_02(N(x,n-2)) + HermiteInterpolation::H5_12(N(x,n-2));
            
    // D_0 H3_0 + D_1 (CD_1-4) H3_1 + D_2 (CD_2-4) H3_1
    return H5_0_n  
           + ( (H5_1_nm2 - H5_1_np2)  -  8. * ( H5_1_nm1 - H5_1_np1 )) / ( 12.* dx)
           + (-(H5_2_nm2 + H5_2_np2 ) + 16. * ( H5_2_nm1 + H5_2_np1) - 30. * H5_2_n ) / (12. * pow2(dx));
  }
  // 7th order interpolation
  else if      (interpolationOrder == 7) {
            check(-1, DMESG("Interpolation order not implemented"));
  }
  else check(-1, DMESG("Interpolation order not implemented"));
        
  return 0.;      
}


void FieldsHermite::printOn(std::ostream &output) const 
{
  output   << "Fields    |  Hermite      Order : "  << interpolationOrder << std::endl;

  return;
}
 

void FieldsHermite::solvePoissonEquation(const CComplex Q     [Nq][NxLD][NkyLD][Nz],
                                               CComplex Field0[Nq][NxLD][NkyLD][Nz])
{

  CComplex *vec_Q, *vec_Field;
     
  for(int z = NzLlD; z <= NzLuD; z++) {  for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) {

    VecGetArray    (GyroVector_X, (PetscScalar **) &vec_Q);
    for(int x = NxLlD, n = 0; x <= NxLuD; x++)  vec_Q[n++] = Q[Field::phi][z][y_k][x];
    VecRestoreArray(GyroVector_X, (PetscScalar **) &vec_Q);
       
    MatrixPoissonSolverLHS(y_k, z)->solve(GyroVector_X, GyroVector_XAvrg);
           
    VecGetArrayRead    (GyroVector_XAvrg, (const PetscScalar **) &vec_Field);
    for(int x=NxLlD, n = 0; x <= NxLuD; x++) Field0[Field::phi][z][y_k][x] = vec_Field[n++]; 
    VecRestoreArrayRead(GyroVector_XAvrg, (const PetscScalar **) &vec_Field);

  } }

  return;

};
           
         
Matrix* FieldsHermite::getGyroAveragingMatrix(const double mu, const int y_k, const int z, const int s)
{

  Matrix *M = new Matrix(NxLD, grid->NxGD, DIR_X, parallel);         
           
  // fill martrix
  for(int x = NxLlD; x <= NxLuD; x++) { for(int n = NxLlD; n <= NxLuD; n++) {
        
    const double rho_t2  = plasma->species[s].T0 * plasma->species[s].m / (pow2(plasma->species[s].q) * plasma->B0); 
    const double lambda2 = 2. * mu * rho_t2;
                
    const CComplex value  = getElements(x, n, sqrt(lambda2), y_k, z);

    // Note that our domain starts at NxGlD, while PETSc col/row @ 0
    M->setValue(x - NxGlD, n - NxGlD, value);
        
  } }

  // MatSetValues (PETSc) caches values thus MatAssembly needs to be called
  M->assemble();

  return M;
}
  
// not implemented
void FieldsHermite::getFieldEnergy(double& phiEnergy, double& ApEnergy, double& BpEnergy) 
{
   phiEnergy = 0.; ApEnergy = 0.; BpEnergy = 0.;

   return ;
};
        

void FieldsHermite::setDoubleGyroAverageMatrix(Setup *setup)
{

  ///////////////////////     Create Double GyroAverageMatrix /////////////////////
  // Note : int_0^infty J cdot J e^{mu} is sensitiv, thus we use higher order matrix
  //        to setup this matrix (linear matrix) Increase by factor 8
  Matrix Matrix_Gamma0    (NxLD, grid->NxGD, DIR_X, parallel);
  Matrix Matrix_PoissonLHS(NxLD, grid->NxGD, DIR_X, parallel);

  int integrationOrder = setup->get("Fields.Hermite.DoubleGyroIntegrationOrder", 32);
        
  Integrate GRQuad("Gauss-Laguerre", integrationOrder);

  parallel->print("Intializing Double Gyro-Averaging Matrix");
        
  for(int z = NzLlD; z <= NzLuD; z++) {  for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { 
              
    // Is out matrix hermitian (propably epends on the geometry) ?!
    Matrix_PoissonLHS.setZero();
           
    // Create 1-\Gamma for each species s
    for(int s = NsLlD; s <= NsLuD; s++) { 
       
      Matrix_Gamma0.setZero();
              
      //  Perform loop for higher order integration of Double Gyro Average Matrix
      //
      // NOT YET Supported by intel: for(int n : Tools::ParallelRange(integrationOrder, parallel, DIR_M)) { 
      auto v = Tools::ParallelRange(integrationOrder, parallel, DIR_M);
      
      for(auto it = v.begin(); it != v.end(); ++it) {
        
        int n = *it;
        const double mu = GRQuad.x(n);
                   
        // ToDO : what if BT != 1 ??!!
        const double BT = plasma->B0 / plasma->species[s].T0; 

        // no exp in case of Laguerre intergration const Complex w  = BT * exp(- BT * mu) * GRQuad.w(n); 
        const CComplex w  = BT * GRQuad.w(n); 
        
        Matrix *GyroM = getGyroAveragingMatrix(mu/BT, y_k, z, s);
                    
        // what about imaginary values
        Matrix_Gamma0 += - w * ( (*GyroM) * (*GyroM) );

        delete GyroM;
        
        // some feedback to used that program is not stalled
        std::cout << "*";
  
      }
      
      Matrix_Gamma0.reduce(DIR_M);
      Matrix_Gamma0.addDiagonal(1.);

      const double qqnT  = plasma->species[s].n0 * pow2(plasma->species[s].q)/plasma->species[s].T0;
      Matrix_PoissonLHS += qqnT * Matrix_Gamma0;
           
    } // s
    
    Matrix_PoissonLHS.reduce(DIR_S); 
           
    // add adibatic response, solve for correction ?!
    const double adiab = plasma->species[0].n0 * pow2(plasma->species[0].q)/plasma->species[0].T0;
    Matrix_PoissonLHS.addDiagonal(adiab);
           
    MatrixPoissonSolverLHS(y_k, z) = new MatrixSolver(parallel, Matrix_PoissonLHS.getMat(), DIR_X, true, "SuperLU");

    // Poisson Matrix (lambda_D2 ) $ g^{xx} D_x^2  + 2 i g^{xy} k_y D_x - g^{yy} k_y^2$
    // ignore now, only do quasi neutrality.
    
  } } // y_k, z

          
}
