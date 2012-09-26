/*
 * =====================================================================================
 *
 *       Filename:  FFTSolver_fftw3.cpp
 *    Description:  Fourier Solver implementation for fftw3-mpi 
 *                  see www.fftw.org
 *         Author:  Paul P. Hilscher (2009, 2011), 
 *
 * =====================================================================================
 */


#include "FFTSolver/FFTSolver_fftw3.h"
#include "Plasma.h"

#include <fftw3-mpi.h>


// we have to place plans here to avoid namespace errors 
fftw_plan plan_YForward_Field, plan_YBackward_Field, 
          plan_YForward_PSF, plan_YBackward_PSF,
          plan_YForward_NL , plan_YBackward_NL;
fftw_plan plan_XForward_Fields, plan_XBackward_Fields;
fftw_plan plan_FieldTranspose_1, plan_FieldTranspose_2;
fftw_plan plan_AA_YForward, plan_AA_YBackward;


// from http://agentzlerich.blogspot.jp/2010/01/using-fftw-for-in-place-matrix.html
// any license issues ? 
fftw_plan plan_transpose(char storage_type, int rows, int cols, double *in, double *out);



FFTSolver_fftw3::FFTSolver_fftw3(Setup *setup, Parallel *parallel, Geometry *geo) : FFTSolver(setup, parallel, geo, Nx*(2*Nky-2)*Nz, Nx*(2*Nky-2), Nx,  (2*Nky-2)) {

  
  if(parallel->Coord[DIR_V] == 0) {          // Fourier solver required in velocity space

  
    // Setup plans
    int perf_flag = FFTW_ESTIMATE;
    
    plan   = setup->get("FFTW3.Plan", "");
    if      (plan == "Estimate"   ) perf_flag = FFTW_ESTIMATE;
    else if (plan == "Measure"    ) perf_flag = FFTW_MEASURE;
    else if (plan == "Exhaustive" ) perf_flag = FFTW_EXHAUSTIVE;

    // Setup wisedom
    wisdom = setup->get("FFTW3.Wisdom", "");
    if     (wisdom == "System") fftw_import_system_wisdom();
    else if(wisdom != ""      ) fftw_import_wisdom_from_filename(wisdom.c_str());

    const int nfields = plasma->nfields;

#ifdef PARALLEL_OPENMP
    fftw_init_threads();
    fftw_plan_with_nthreads(parallel->numThreads);
#endif

#ifdef GKC_PARALLEL_MPI
    fftw_mpi_init();
#endif
    
    // needs for poissons equation
    
         
    
      // set and check bounds 
      long X_NxLD, X_NxLlD, X_NkxL, X_NkxLlD, X_numElements, X_Nx = Nx; 
      
      X_numElements = fftw_mpi_local_size_1d(Nx, parallel->Comm[DIR_X], FFTW_FORWARD, 0, &X_NxLD, &X_NxLlD, &X_NkxL, &X_NkxLlD);
      
      FFTSolver::X_NkxL = X_NkxL;
      
      // Prefactor of 2 for safety (is required otherwise we get crash, but why ?)
      int numAlloc = 2 * X_numElements * NkyLD * NzLD * nfields;

      // allocate arrays 
      data_kXIn       = (CComplex *) fftw_alloc_complex(numAlloc);
      data_kXOut      = (CComplex *) fftw_alloc_complex(numAlloc);
      data_X_rOut     = (CComplex *) fftw_alloc_complex(numAlloc);
      data_X_rIn      = (CComplex *) fftw_alloc_complex(numAlloc);
      data_X_Transp_1 = (CComplex *) fftw_alloc_complex(numAlloc);
      data_X_Transp_2 = (CComplex *) fftw_alloc_complex(numAlloc);
      
      check((NxLD != X_NxLD) ? -1 : 0, DMESG("Bounds to not align")); 
         
      // set and check bounds 
      K1xLlD = X_NkxLlD;       K1xLuD = X_NkxLlD + X_NkxL - 1;
      Rk1xL.setRange(K1xLlD, K1xLuD);

      // used only to calculate offset
      nct::allocate Array_kX = nct::allocate(nct::Range(1,plasma->nfields), nct::Range(NzLlD, NzLD), nct::Range(NkyLlD,NkyLD), nct::Range(X_NkxLlD, X_NkxL));
      
      // Calculated shifted pointer for CEAN
      kXIn  = Array_kX.zero(data_kXIn );
      kXOut = Array_kX.zero(data_kXOut);

      // Note : We should use transformed out to improve parallelization
      // fftw_plan fftw_mpi_plan_many_dft(int rnk, const ptrdiff_t *n, 
      //           ptrdiff_t howmany, ptrdiff_t block, ptrdiff_t tblock, fftw_complex *in, fftw_complex *out,
      //           MPI_Comm comm, int sign, unsigned flags);
      
      long numTrans = NkyLD * NzLD * nfields;
      
      plan_XForward_Fields  = fftw_mpi_plan_many_dft(1, &X_Nx, numTrans, NxLD, X_NkxL, (fftw_complex *) data_X_rIn, (fftw_complex *) data_kXOut , parallel->Comm[DIR_X], FFTW_FORWARD , perf_flag);
      plan_XBackward_Fields = fftw_mpi_plan_many_dft(1, &X_Nx, numTrans, NxLD, X_NkxL, (fftw_complex *) data_kXIn , (fftw_complex *) data_X_rOut, parallel->Comm[DIR_X], FFTW_BACKWARD, perf_flag);

      // Fields have to be continuous in howmanyfields, and thus we have to transpose the array (use in-place) 
      // add factor of 2 because we deal with complex numbers not real numbers
      // plan_FieldTranspose = plan_transpose('R', 2* NxLD, 2 * NkyLD * NzLD * nfields, (double *) kXOut.data(), (double *) kXOut.data());
      // plan_FieldTranspose_1 = plan_transpose('R', 2 * NkyLD * NzLD * nfields, 2*   NxLD, (double *) data_X_Transp_1, (double *) data_X_Transp_2);
      // plan_FieldTranspose_2 = plan_transpose('C', 2 * NkyLD * NzLD * nfields, 2 * X_NkxL, (double *) data_X_Transp_1, (double *) data_X_Transp_2);
      // check(((plan_FieldTranspose_1 == NULL) || (plan_FieldTranspose_2 == NULL)) ? -1 : 0, DMESG("Transpose planner null"));

           
            
   
    // Needed to calculate non-linearity in real space
    
    ////////////////////////////////// Define Fourier transforms for Y ///////////////// 
    
      //                                                  howmany                                  stride distance
      //  leave space for boundary conditions
      // 
      //   kYIn[Nky][Nx]
      //
      //   [ y0x0 y0x1 y0x2 ... ] [y1x0 y1x1 y1x2 ...] [y2x0 y2x1 ... ] ...
      //
      //   thus our FFT routine has stride of NxLD, however our output is same
      //
      //   [ y0x0 y0x1 y0x2 ... ] [y1x0 y1x1 y1x2 ...] [y2x0 y2x1 ... ] ...
      //
      //   
      //    from fftw3-doc :
      //
      //       fftw_plan fftw_plan_many_dft(int rank, const int *n, int howmany,
      //                                    fftw_complex *in, const int *inembed,
      //                                    int istride, int idist,
      //                                    fftw_complex *out, const int *onembed,
      //                                    int ostride, int odist,
      //                                    int sign, unsigned flags);
      //
      //      location of input : in + k * idist
      //      The stride parameters indicate that the j-th element of the input or output arrays is
      //      located at j*istride 
      //
      //
      //   Thus :
      //           X is on fast dimension, thus our fourier mode is the stride NxLD+BD
      //           Distance for the next mode is thus 1          
      //
      
      perf_flag |= FFTW_UNALIGNED;
      
      // Orginal
      
      double   rY_BD4[NyLD ][NxLB+4]; CComplex kY_BD4[NkyLD][NxLB+4];
      plan_YBackward_Field = fftw_plan_many_dft_c2r(1, &NyLD, NxLB+4, (fftw_complex *) kY_BD4, NULL, NxLB+4, 1, (double *) rY_BD4, NULL, NxLB+4, 1, perf_flag);
      

      
      double   rY_BD2[NyLD ][NxLB]; CComplex kY_BD2[NkyLD][NxLB];
      plan_YBackward_PSF   = fftw_plan_many_dft_c2r(1, &NyLD, NxLB  , (fftw_complex *) kY_BD2, NULL, NxLB  , 1, (double *) rY_BD2, NULL, NxLB  , 1, perf_flag);
      

      
      double   rY_BD0[NyLD ][NxLD]; CComplex kY_BD0[NkyLD][NxLD];
      plan_YForward_NL     = fftw_plan_many_dft_r2c(1, &NyLD, NxLD, (double *     ) rY_BD0, NULL, NxLD, 1, (fftw_complex *) kY_BD0, NULL, NxLD  , 1, perf_flag);
      plan_YBackward_NL    = fftw_plan_many_dft_c2r(1, &NyLD, NxLD, (fftw_complex*) kY_BD0, NULL, NxLD, 1, (double       *) rY_BD0, NULL, NxLD , 1, perf_flag);
               
      
      ////////////////////////   Define Anti-Aliased Arrays /////////////////////////////////////
      
      int AA_NkyLD  = 3 * Nky   / 2       , 
          AA_NkyLlD = NkyLlD              , 
          AA_NkyLuD = NkyLlD + AA_NkyLD -1,
          AA_NyLD  = 2 * AA_NkyLD - 2     ;
                
          
      double   rY_AA[AA_NyLD][NxLD]; CComplex kY_AA[AA_NkyLD][NxLD];
      

      
      plan_AA_YForward  = fftw_plan_many_dft_r2c(1, &AA_NyLD, NxLD, (double      *) rY_AA, NULL, 1, AA_NyLD ,  (fftw_complex*) kY_AA, NULL, 1, AA_NkyLD, perf_flag);
      plan_AA_YBackward = fftw_plan_many_dft_c2r(1, &AA_NyLD, NxLD, (fftw_complex*) kY_AA, NULL, 1, AA_NkyLD,  (double      *) rY_AA, NULL, 1, AA_NyLD , perf_flag);


    setNormalizationConstants();
   
    // export it again (only by root job ?!)
   
    if(wisdom != "") fftw_export_wisdom_to_filename(wisdom.c_str());
   
  }

}


/// too much crap here ..... :(
void FFTSolver_fftw3::solve(const FFT_Type type, const FFT_Sign direction, void *in, void *out) 
{

   if(type == FFT_Type::X_FIELDS) {
             
   
     if(in == nullptr)  check(-1, DMESG("Need Pointer to array"));

     
     if     (direction == FFT_Sign::Forward )  {
    
       // fftw3-mpi many transform requires specific input (thus we have to transpose our data)
       transpose(NxLD, NkyLD, NzLD, plasma->nfields      , (A4zz) ((CComplex *) in)             , (A4zz) ((CComplex *) data_X_Transp_1));                
       fftw_mpi_execute_dft(plan_XForward_Fields , (fftw_complex *) data_X_Transp_1 , (fftw_complex *) data_X_Transp_2); 
       transpose_rev(X_NkxL, NkyLD, NzLD, plasma->nfields, (A4zz) ((CComplex *) data_X_Transp_2), (A4zz) ((CComplex *) data_kXOut));               

     }
    
     else if(direction == FFT_Sign::Backward) {
                  
       // fftw3-mpi many transform requires specific input (thus we have to transpose our data and backtransform)
       transpose(X_NkxL, NkyLD, NzLD, plasma->nfields, (A4zz) ((CComplex *) data_kXIn), (A4zz) ((CComplex *) data_X_Transp_1));                
       fftw_mpi_execute_dft(plan_XBackward_Fields, (fftw_complex *) data_X_Transp_1, (fftw_complex *) data_X_Transp_2 ); 
       transpose_rev(NxLD, NkyLD, NzLD, plasma->nfields, (A4zz) ((CComplex *) data_X_Transp_2), (A4zz) ((CComplex *) in));                
     }
     
     else   check(-1, DMESG("No such FFT direction"));

     
   }  
  
   // These are speed critical (move above x-transformation)
   else if(type == FFT_Type::Y_FIELDS ) {
            
             // Need to cast between bit comparible complex types
             //if     (direction == FFT_Sign::Forward )  fftw_execute_dft_r2c(plan_YForward_Field , (double   *) in, (fftw_complex *) out); 
             if(direction == FFT_Sign::Backward)  fftw_execute_dft_c2r(plan_YBackward_Field, (fftw_complex *) in, (double   *) out); 
             else   check(-1, DMESG("No such FFT direction"));
   }
        
   
   else if(type == FFT_Type::Y_PSF ) {
            
             // Need to cast between bit comparible complex types
             //if     (direction == FFT_Sign::Forward )  fftw_execute_dft_r2c(plan_YForward_PSF , (double   *) in, (fftw_complex *) out); 
             if(direction == FFT_Sign::Backward)  fftw_execute_dft_c2r(plan_YBackward_PSF, (fftw_complex *) in, (double   *) out); 
             else   check(-1, DMESG("No such FFT direction"));
   
   }
   
   else if(type == FFT_Type::Y_NL  ) {
            
             // Need to cast between bit comparible complex types
             if     (direction == FFT_Sign::Forward )  fftw_execute_dft_r2c(plan_YForward_NL , (double   *) in, (fftw_complex *) out); 
             else if(direction == FFT_Sign::Backward)  fftw_execute_dft_c2r(plan_YBackward_NL, (fftw_complex *) in, (double   *) out); 
             else   check(-1, DMESG("No such FFT direction"));
   
   }

   else  check(-1, DMESG("Unknown FFT type or not supported"));
   
   return;
}


std::string FFTSolver_fftw3::getLibraryName() 
{
        return std::string(fftw_version) + std::string("-mpi");
}
  

FFTSolver_fftw3::~FFTSolver_fftw3()
{

    //  release fftw-3  
        fftw_destroy_plan(plan_XForward_Fields);
       fftw_destroy_plan(plan_XBackward_Fields);
       
       fftw_destroy_plan(plan_YForward_Field);
       fftw_destroy_plan(plan_YBackward_Field);
        
       fftw_destroy_plan(plan_AA_YForward);
       fftw_destroy_plan(plan_AA_YBackward);

#ifdef PARALLEL_OPENMP
    fftw_cleanup_threads();
#endif
    fftw_free(data_X_rOut);
    fftw_free(data_X_rIn );
    fftw_free(data_kXOut );
    fftw_free(data_kXIn  );
    

    fftw_destroy_plan(plan_FieldTranspose_1);
    fftw_destroy_plan(plan_FieldTranspose_2);
}


// Nice, no X-parallelization required !!
void FFTSolver_fftw3::multiply(const CComplex A[NkyLD][NxLD], const CComplex B[NkyLD][NxLD],
                               CComplex R[NkyLD][NxLD])
{
   

   // Create antialiased arrays and copy and transform to real space
   CComplex AA_A[AA_NkyLD][NxLD];
   double   RS_A[AA_NyLD ][NxLD], RS_B[AA_NyLD ][NxLD];
   
 
   // Copy to larger Anti-Aliased array and transform to real space
   AA_A[NkyLlD:NkyLD][:] = A[NkyLlD:NkyLD][:];
   fftw_execute_dft_c2r(plan_AA_YBackward, (fftw_complex *) AA_A, (double *) RS_A);

   AA_A[NkyLlD:NkyLD][:] = B[NkyLlD:NkyLD][:];
   fftw_execute_dft_c2r(plan_AA_YBackward, (fftw_complex *) AA_A, (double *) RS_B);

   //////////////////////// Real Space (multiply values) /////////////////// 
   
   const double _kw_fft_Norm = 1./(1.5 * Norm_Y_Forward * pow2(Norm_Y_Backward));
   
   RS_A[:][:] *= RS_B[:][:] * _kw_fft_Norm;
   
   //////////////////////// End Real Space (multiply values) /////////////////// 
  
   fftw_execute_dft_r2c(plan_AA_YForward, (double *) RS_A, (fftw_complex *) R);
   
   return;

};

/*  
void somelogic()
{
    char storage_type = 'R';
    int rows          = 3;
    int cols          = 3;
    double *in        = p;
    double *out       = p2;

    // Plan the transpose once; transpose is in-place if in == out 
    fftw_plan transpose = plan_transpose(storage_type, rows, cols, in, out);
    assert(transpose);

    // Execute the plan potentially many times 
    fftw_execute(transpose);

    // FFTW New-array Execute functions should be callable, too 
    // Beware of mixing in-place and out-of-place planning and usage 
    double *another_in  = a;
    double *another_out = b;
    fftw_execute_r2r(transpose, another_in, another_out);

    // Destroy the plan when completely done 
    fftw_destroy_plan(transpose);
}
 *  */


fftw_plan plan_transpose(char storage_type, int rows, int cols, double *in, double *out)
{
    const unsigned flags = FFTW_MEASURE; /* Do not destroy input */

    fftw_iodim howmany_dims[2];
    switch (toupper(storage_type)) {
        case 'R':
            howmany_dims[0].n  = rows;
            howmany_dims[0].is = cols;
            howmany_dims[0].os = 1;
            howmany_dims[1].n  = cols;
            howmany_dims[1].is = 1;
            howmany_dims[1].os = rows;
            break;
        case 'C':
            howmany_dims[0].n  = rows;
            howmany_dims[0].is = 1;
            howmany_dims[0].os = cols;
            howmany_dims[1].n  = cols;
            howmany_dims[1].is = rows;
            howmany_dims[1].os = 1;
            break;
        default:
            return NULL;
    }
    const int howmany_rank = sizeof(howmany_dims)/sizeof(howmany_dims[0]);

    return fftw_plan_guru_r2r(0, NULL, 2, howmany_dims, in, out, NULL, flags);
    //return fftw_plan_guru_r2r(/*rank*/0, /*dims*/NULL, howmany_rank, howmany_dims, in, out, /*kind*/NULL, flags);
}

    
void FFTSolver_fftw3::printOn(std::ostream &output) const {

         output   << "FFTSolver  |  using fftw-3 interface for (" << std::string(fftw_version) << ")" << std::endl;
         output   << "           |  Plan : " << plan << " Wisdom : " << ((wisdom=="") ? "None" : wisdom) << std::endl;
         
}



// restrict pointers
// transpose from C-ordering to Fortran ordering (required by fftw3-mpi) 
void FFTSolver_fftw3::transpose(int Nx, int Ny, int Nz, int Nq, CComplex In[Nq][Nz][Ny][Nx], CComplex OutT[Nx][Ny][Nz][Nq])
{
     #pragma ivdep
     for(int x=0; x < Nx; x++ ) { for(int y=0; y < Ny; y++ ) {  
     for(int z=0; z < Nz; z++ ) { for(int q=0; q < Nq; q++ ) {  

        OutT[x][y][z][q] = In[q][z][y][x];
   
     } } } }
  
}

// transpose from Fortran ordering to C-ordering (required by fftw3-mpi) 
void FFTSolver_fftw3::transpose_rev(int Nx, int Ny, int Nz, int Nq, CComplex In[Nx][Ny][Nz][Nq], CComplex OutT[Nq][Nz][Ny][Nx])
{
     #pragma ivdep
     for(int x=0; x < Nx; x++ ) { for(int y=0; y < Ny; y++ ) {  
     for(int z=0; z < Nz; z++ ) { for(int q=0; q < Nq; q++ ) {  

        OutT[q][z][y][x] = In[x][y][z][q];
   
       } } } }
  
}

