/*
 * =====================================================================================
 *
 *       Filename: FieldsHypre.cpp
 *
 *    Description: Interface to Hypre Matrix Solver routine
 *
 *         Author: Paul P. Hilscher (2011), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include "FieldsHypre.h"

FieldsHypre::FieldsHypre(Setup *setup, Grid *grid, Parallel *parallel, FileIO *fileIO, Geometry<HELIOS_GEOMETRY> *geo, FFTSolver *fft) : FieldsFFT(setup, grid, parallel,fileIO, geo, fft)
{
  
   // setup grid 
   NXky_LlD[0] =  NxLlD; NXky_LlD[1] =  NkyLlD;
   NXky_LuD[0] =  NxLuD; NXky_LuD[1] =  NkyLuD;

     
   HYPRE_StructGridCreate(parallel->Comm[DIR_XYZ], 2, &poisson_grid);
   HYPRE_StructGridSetExtents(poisson_grid, NXky_LlD, NXky_LuD);


   //int BD_periodic[] = { NxGuD, 0 };
   //HYPRE_StructGridSetPeriodic (poisson_grid, BD_periodic);


   HYPRE_StructGridAssemble(poisson_grid);
   
   // set Matrixes
   setupMatrixPoissonEquation(setup->get("Hypre.PoissonEquationOrder", 2));


   // Setup solution vector and b
   HYPRE_StructVectorCreate(parallel->Comm[DIR_XYZ], poisson_grid, &vec_b);
   HYPRE_StructVectorCreate(parallel->Comm[DIR_XYZ], poisson_grid, &vec_x);

   HYPRE_StructVectorInitialize(vec_b);
   HYPRE_StructVectorInitialize(vec_x);
     
   
   // Initialize solver
   tolerance = setup->get("Hypre.Tolerance", 1.0e-12);
   
   HYPRE_StructPCGCreate(parallel->Comm[DIR_XYZ], &solver);
   //HYPRE_StructGMRESCreate(parallel->Comm[DIR_XYZ], &solver);
   HYPRE_StructPCGSetTol(solver, tolerance);
   HYPRE_StructPCGSetPrintLevel(solver, 1);
   HYPRE_StructPCGSetup(solver, A, vec_b, vec_x);
   //HYPRE_StructGMRESSetup(solver, A, vec_b, vec_x);

   HYPRE_StructDiagScaleSetup(solver, A,vec_b, vec_x);
   // initialize Arrays
   rho_star.resize(RxLD, RkyLD, RzLD);
}

FieldsHypre::~FieldsHypre()
{
   HYPRE_StructPCGDestroy(solver);
   HYPRE_StructVectorDestroy(vec_b);
   HYPRE_StructVectorDestroy(vec_x);
   HYPRE_StructMatrixDestroy(A);
   HYPRE_StructGridDestroy(poisson_grid);
}




Array3z FieldsHypre::solvePoissonEquation(Array3z rho, Timing timing)
{
   double values_I[NxLD * NkyLD]; 
   double values_R[NxLD * NkyLD]; 
    
   HYPRE_StructVectorAssemble(vec_x);

   for(int z = NzLlD; z <= NzLuD; z++) {

        calcRhoStar(Q);
        // We need to solve real/imag values separately
        for(int x = NxLlD, k = 0; x <= NxLuD; x++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++, k++) { 
            values_R[k] = real(rho_star(x, y_k, z));
            values_I[k] = imag(rho_star(x, y_k, z));
        }}

        // solve Real
        HYPRE_StructVectorSetBoxValues(vec_b, NXky_LlD, NXky_LuD, values_R);
        HYPRE_StructVectorAssemble(vec_b);
    
        HYPRE_StructDiagScaleSetup(solver, A,vec_b, vec_x);
        HYPRE_StructPCGSetup(solver, A, vec_b, vec_x);
        HYPRE_StructPCGSolve(solver, A, vec_b, vec_x);
        HYPRE_StructVectorGetBoxValues(vec_x, NXky_LlD, NXky_LuD, values_R);
    
        // solve Imag
        HYPRE_StructVectorSetBoxValues(vec_b, NXky_LlD, NXky_LuD, values_I);
        HYPRE_StructVectorAssemble(vec_b);
   
        HYPRE_StructDiagScaleSetup(solver, A,vec_b, vec_x);
        HYPRE_StructPCGSetup(solver, A, vec_b, vec_x);
        HYPRE_StructPCGSolve(solver, A, vec_b, vec_x);
        HYPRE_StructVectorGetBoxValues(vec_x, NXky_LlD, NXky_LuD, values_I);
        
        // back to real imaginary
        for(int x = NxLlD, k = 0; x <= NxLuD; x++)  { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++, k++) {
            Field0(x, y_k, z, Field::phi) = cmplxd(values_R[k],values_I[k]);
        } }
  
   }
    
   return rho;

}

   
void FieldsHypre::setupMatrixPoissonEquation(int n_points) {


        HYPRE_StructStencil  stencil;
        HYPRE_StructStencilCreate(2, n_points, &stencil);
        int offsets[4][2] = {{0,0}, {1,0}, {2,0}, {3,0}};

        // Assign each of the 5 stencil entries 
        for (int entry = 0; entry < n_points; entry++) HYPRE_StructStencilSetElement(stencil, entry, offsets[entry]);
          
        HYPRE_StructMatrixCreate(parallel->Comm[DIR_XYZ], poisson_grid, stencil, &A);
        HYPRE_StructMatrixSetSymmetric(A, 1);
        HYPRE_StructMatrixInitialize(A);

        int stencil_indices[] = {0,1,2,3,4,5}; 
        int nvalues           = n_points * NxLD * NkyLD;
 
        double values[nvalues];

        /// create stencil set
        //  only works in certrain geometries (2D Sheared Slab, Shearless Slab)!
        const double norm   = 1.;//plasma->species(1).n0 * pow2(plasma->species(1).q)/plasma->species(1).T0;
        const double rho_t2 = 1.;//plasma->species(1).T0 * plasma->species(1).m / (pow2(plasma->species(1).q) * plasma->B0);
        const double adiab  = 1.;//plasma->species(0).n0 * pow2(plasma->species(0).q)/plasma->species(0).T0;
     
        // Setup the Matrix
        for(int x = NxLlD, k = 0; x <= NxLuD; x++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) {

            // set x-zero BC (if stencil reaches to the boundary, we need to set it zero)
            const double ky2 = pow2(fft->ky(y_k));
           
            // 2nd order discretization
            if(n_points == 2) {
                const double h = dx*dx;
                //values[k++] = ((x+1) >= NxGlD) ? 1. + (2.+plasma->debye2) * (  -2./h + ky2) : 0.;
                //values[k++] = ((x+1) <= NxGuD) ?      (2.+plasma->debye2)/h                  : 0.;
                values[k++] = (1) ? 1. + (2.+plasma->debye2) * ( 2./h - ky2) : 0.;
                values[k++] = (1) ?      -(2.+plasma->debye2)/h                  : 0.;
            }
            // 4th order discretization
            else if(n_points == 3) {
                const double h = 12.*dx*dx;
                values[k++] = (plasma->debye2 + 2.) + 2./h + ky2 ;
                values[k++] = ((x+1) <= NxGuD) ? -16./h : 0.;
                values[k++] = ((x+2) <= NxGuD) ?   1./h : 0.;
            }
            // 6th order discretization
            else if(n_points == 4) {
                const double h = 12.*dx*dx;
                values[k++] = plasma->debye2 + 2. + 2./h + ky2 ;
                values[k++] = ((x+1) != NxGuD) ? -1./h : 0.;
                values[k++] = ((x+2) != NxGuD) ? -1./h : 0.;
                values[k++] = ((x+3) != NxGuD) ? -1./h : 0.;
            }
            else check(-1, DMESG("No such order"));

        }}

        HYPRE_StructMatrixSetBoxValues(A, NXky_LlD, NXky_LuD, n_points, stencil_indices, values);
        HYPRE_StructMatrixAssemble(A);
        HYPRE_StructStencilDestroy(stencil);
 }
   
Array4z FieldsHypre::solveFieldEquations(Array4z Q, Timing timing) 
{
   
  
   // for 3 fields phi and B_perp are coupled and need to be solved differently
   if( solveEq & Field::phi) solvePoissonEquation  (Q(RxLD, RkyLD, RzLD, Field::phi), timing);
   if( solveEq & Field::Ap ) check(-1, DMESG("Not Implemented"));
   if( solveEq & Field::Bpp) check(-1, DMESG("Not Impleemnted"));
    
   return Field0(RxLD, RkyLD, RzLD, RFields);

}



void FieldsHypre::calcRhoStar(Array4z Q) 
{
    
    //if(plasma->species(0).doGyro)  phi_yz = calcFluxSurfAvrg(fft->kXOut);

    
//    rho_star(RxLD, RkyLD, RzLD) = Q(RxLD, RkyLD, RzLD, Field::phi);
//    return;

    fft->rXIn(RxLD, RkyLD, RzLD, RFields) = Q(RxLD, RkyLD, RzLD, RFields);
    fft->solve(FFT_X, FFT_FORWARD, FFT_FIELDS);

    const double rho_t2 = plasma->species(1).T0  * plasma->species(1).m / pow2(plasma->species(1).q * plasma->B0);
    
    #pragma omp parallel for
    for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int z=NzLlD; z<=NzLuD;z++) { for(int x_k=fft->K1xLlD; x_k<= fft->K1xLuD;x_k++) {

          if((x_k == 0) && (y_k == 0)) { fft->kXIn(x_k,y_k,z, Field::phi) = (cmplxd) 0.e0 ; continue; }
          
          const double b = rho_t2 * fft->k2_p(x_k,y_k,z);
          
          //const cmplxd rhs  = (fft->kXOut(x_k,y_k,z, Q::rho) + adiab * phi_yz(x_k))/fft->Norm_X; 
          const cmplxd rhs  = fft->kXOut(x_k,y_k,z, Q::rho)/fft->Norm_X; 
          
          fft->kXIn(x_k,y_k,z, Field::phi) = (1. + b) * rhs; 
         
    } } }
   
    fft->solve(FFT_X, FFT_BACKWARD, NkyLD * NzLD * plasma->nfields);
   
    rho_star(RxLD, RkyLD, RzLD) = fft->rXOut(RxLD, RkyLD, RzLD, Field::phi);

    return;
}

/* 
void FieldsHypre::calcRhoStar(Array4z Q) 
{


      for(int z = NzLlD; z <= NzLuD; z++) {
      for(int x = NxLlD; x <= NxLuD; x++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) {
                
                const cmplxd d2phi_dx = (16. * (Q(x+1,y_k,z,Field::phi) + Q(x-1,y_k,z,Field::phi)) - (Q(x+2,y_k,z,Field::phi) + Q(x-2,y_k,z,Field::phi)) - 30. * Q(x,y_k,z,Field::phi)) / (12. * pow2(dx));
                const cmplxd d2phi_dy = pow2(fft->ky(y_k)) * Q(x,y_k,z,Field::phi);

                rho_star(x,y_k,z) = Q(x,y_k,z,Field::phi) + (d2phi_dx + d2phi_dy);
      }}

    }

}
 * */
