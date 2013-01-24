/*
 * =====================================================================================
 *
 *       Filename: Vlasov.cpp
 *
 *    Description: Vlasov Solver Interface
 *
 *         Author: Paul P. Hilscher (2009-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include "Vlasov.h"


Vlasov::Vlasov(Grid *_grid, Parallel *_parallel, Setup *_setup, FileIO *fileIO, Geometry *_geo, FFTSolver *_fft, Benchmark *_bench, Collisions *_coll)

: fft(_fft), bench(_bench), parallel(_parallel), grid(_grid), setup(_setup), geo(_geo), coll(_coll)

, _kw_12_dx_dx( 1./(12.*dx*dx)   )
, _kw_12_dv   ( 1./(12.*dv)      )
, _kw_12_dx   ( 1./(12.*dx)      )
, _kw_12_dz   ( 1./(12.*dz)      )
, _kw_12_dv_dv( 1./(12.*dv*dv)   )
, _kw_16_dx4  ( 1./(16.*pow4(dx)))


{

   ArrayPhase = nct::allocate(grid->RsLD, grid->RmLD, grid->RzLB, grid->RkyLD, grid->RxLB, grid->RvLB);
   ArrayPhase(&f0, &f, &fss, &fs, &f1, &ft, &Coll);
   
   ArrayXi = nct::allocate(grid->RzLB , grid->RkyLD, grid->RxLB4, grid->RvLB)(&Xi);
   ArrayG  = nct::allocate(grid->RzLB , grid->RkyLD, grid->RxLB , grid->RvLB)(&G );
   ArrayNL = nct::allocate(grid->RkyLD, grid->RxLD , grid->RvLD)(&nonLinearTerm);
   
   // allocate boundary (mpi) buffers
   int BoundX_num =    2 * Nky * NzLD * NvLD * NmLD * NsLD;
   int BoundZ_num = NxLD * Nky *    2 * NvLD * NmLD * NsLD;
   int BoundV_num = NxLD * Nky * NzLD *    2 * NmLD * NsLD;

   ArrayBoundX = nct::allocate(nct::Range(0 , BoundX_num ))(&SendXu, &SendXl, &RecvXl, &RecvXu);
   ArrayBoundZ = nct::allocate(nct::Range(0 , BoundZ_num ))(&SendZu, &SendZl, &RecvZl, &RecvZu);
   ArrayBoundV = nct::allocate(nct::Range(0 , BoundV_num ))(&SendVu, &SendVl, &RecvVl, &RecvVu);
  
   equation_type       = setup->get("Vlasov.Equation" , "ES");        
   doNonLinear         = setup->get("Vlasov.doNonLinear", 0      );
   doNonLinearParallel = setup->get("Vlasov.doNonLinearParallel", 0);

   const std::string dir_string[] = { "X", "Y", "Z", "V", "M", "S" };
   for(int dir = DIR_X; dir <= DIR_S; dir++) hyp_visc[dir] = setup->get("Vlasov.HyperViscosity." + dir_string[dir], 0.0);
   
   dataOutputF1      = Timing( setup->get("DataOutput.Vlasov.Step", -1),
                               setup->get("DataOutput.Vlasov.Time", -1.));

   Xi_max[:] = 0.;

   initData(setup, fileIO);
}



Vlasov::~Vlasov() 
{
   closeData();
};


void Vlasov::solve(Fields *fields, CComplex  *_fs, CComplex  *_fss, 
                   double dt, int rk_step, const double rk[3], bool useNonBlockingBoundary)
{
   static CComplex *f_boundary = nullptr;

   useNonBlockingBoundary =  false;

   // Need boundary_isclean to avoid deadlock at first iteration
   #pragma omp single nowait
   {
   if((f_boundary != nullptr) && useNonBlockingBoundary) setBoundary(f_boundary, Boundary::RECV);
   }

   Xi_max[:] = 0.; // Needed to calculate CFL time step 
  
   // Calculate the collision operator
   // BUG : how to deal with velocity space decomposition ?!
   coll->solve(fields, _fs, f0, Coll, dt, rk_step);
       
   // Calculate the Vlasov equation
   #pragma omp barrier
   solve(equation_type, fields, _fs, _fss, dt, rk_step, rk);


   // Note : we have non-blocking boundaries as Poisson solver does not require ghosts
   // Set nowait, as field solver does not require boundaries & analysis too ... ?! 
   #pragma omp single nowait
   {
   (useNonBlockingBoundary) ? setBoundary(_fss, Boundary::SEND) :  setBoundary(_fss, Boundary::SENDRECV); 
   }

   f_boundary = _fss; 
  
   return;
}

void Vlasov::setBoundary(CComplex *f) 
{ 
  setBoundary(f, Boundary::SENDRECV); 
};


void Vlasov::setBoundary(CComplex *f, Boundary boundary_type)
{
    [=] (
         CComplex g     [NsLD][NmLD][NzLB][Nky][NxLB][NvLB],
         CComplex SendXl[NsLD][NmLD][NzLD][Nky][GC2 ][NvLD], CComplex SendXu[NsLD][NmLD][NzLD][Nky][GC2 ][NvLD], 
         CComplex RecvXl[NsLD][NmLD][NzLD][Nky][GC2 ][NvLD], CComplex RecvXu[NsLD][NmLD][NzLD][Nky][GC2 ][NvLD], 
         CComplex SendZl[NsLD][NmLD][GC2 ][Nky][NxLD][NvLD], CComplex SendZu[NsLD][NmLD][GC2 ][Nky][NxLD][NvLD],
         CComplex RecvZl[NsLD][NmLD][GC2 ][Nky][NxLD][NvLD], CComplex RecvZu[NsLD][NmLD][GC2 ][Nky][NxLD][NvLD],
         CComplex SendVl[NsLD][NmLD][NzLD][Nky][NxLD][GC2 ], CComplex SendVu[NsLD][NmLD][NvLD][Nky][NxLD][GC2 ],
         CComplex RecvVl[NsLD][NmLD][NzLD][Nky][NxLD][GC2 ], CComplex RecvVu[NsLD][NmLD][NvLD][Nky][NxLD][GC2 ])
  {

  /////////////////////////// Send Boundaries //////////////////////////////
  if(boundary_type & Boundary::SEND) {
   
    // X-Boundary (Note, we may have different boundaries for global simulations)
    SendXl[:][:][:][:][:][:] = g[NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][:][NxLlD  :2][NvLlD:NvLD];
    SendXu[:][:][:][:][:][:] = g[NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][:][NxLuD-1:2][NvLlD:NvLD];
    parallel->updateBoundaryVlasov(Vlasov::SendXu, Vlasov::SendXl, Vlasov::RecvXu, Vlasov::RecvXl, ArrayBoundX.getNum(), DIR_X);
   
    // We do not domain decompose poloidal (y) fourier modes, thus boundaries not required
  
    // Z-Boundary (skip exchange in case we use only 2D (Nz == 1) simulations 
    if(Nz > 1) {
    for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { for(int x = NxLlD; x <= NxLuD; x++) { 

      const CComplex shift_l = (NzLlD == NzGlD) ? cexp(-_imag * 2.* M_PI * geo->nu(x)) : 1.;
      const CComplex shift_u = (NzLuD == NzGuD) ? cexp(+_imag * 2.* M_PI * geo->nu(x)) : 1.;
            
      // NzLlD == NzGlD -> Connect only physcial boundaries after mode made one loop 
      SendZl[:][:][:][y_k][x-NxLlD][:] = g[NsLlD:NsLD][NmLlD:NmLD][NzLlD  :2][y_k][x][NvLlD:NvLD] * shift_l;
      SendZu[:][:][:][y_k][x-NxLlD][:] = g[NsLlD:NsLD][NmLlD:NmLD][NzLuD-1:2][y_k][x][NvLlD:NvLD] * shift_u;
               
    } }
    parallel->updateBoundaryVlasov(Vlasov::SendZu, Vlasov::SendZl, Vlasov::RecvZu, Vlasov::RecvZl, ArrayBoundZ.getNum(), DIR_Z);
    } // (Nz > 1)

    // We do not need to communicate for M and S as we do not have boundary cells (yet)
  
    // Decomposition in velocity is rather unlikely, thus give non-decomposed version too
    // We set endpoint in velocity space to zero
    if(parallel->decomposition[DIR_V] > 1) {
      SendVl[:][:][:][:][:][:] = g[NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][:][NxLlD:NxLD][NvLlD  :2]; 
      SendVu[:][:][:][:][:][:] = g[NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][:][NxLlD:NxLD][NvLuD-1:2]; 
      parallel->updateBoundaryVlasov(Vlasov::SendVu, Vlasov::SendVl, Vlasov::RecvVu, Vlasov::RecvVl, ArrayBoundV.getNum(), DIR_V);
    } else {
      g[NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][:][NxLlD:NxLD][NvLlB  :2] = 0.;
      g[NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][:][NxLlD:NxLD][NvLuD+1:2] = 0.;
    }

  }

  /////////////////////////// Receive Boundaries //////////////////////////////
  if(boundary_type & Boundary::RECV) {
     
    // Wait until boundaries are communicated
    parallel->updateBoundaryVlasovBarrier();
   
    // Set boundary in X (take care of Neumann boundary ?!) 
    g[NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][:][NxLlB  :2][NvLlD:NvLD] = RecvXl[:][:][:][:][:][:]; 
    g[NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][:][NxLuD+1:2][NvLlD:NvLD] = RecvXu[:][:][:][:][:][:]; 
    
    // Set boundary in Z
    if(Nz > 1) {
      g[NsLlD:NsLD][NmLlD:NmLD][NzLlB  :2][:][NxLlD:NxLD][NvLlD:NvLD] = RecvZl[:][:][:][:][:][:]; 
      g[NsLlD:NsLD][NmLlD:NmLD][NzLuD+1:2][:][NxLlD:NxLD][NvLlD:NvLD] = RecvZu[:][:][:][:][:][:]; 
    }

    // Set boundary in V
    if(parallel->decomposition[DIR_V] > 1) {
      g[NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][:][NxLlD:NxLD][NvLlB  :2] = RecvVl[:][:][:][:][:][:]; 
      g[NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][:][NxLlD:NxLD][NvLuD+1:2] = RecvVu[:][:][:][:][:][:]; 
    }
   
  }

  }  ((A6zz) f, 
      (A6zz) SendXl,  (A6zz) SendXu,  (A6zz) RecvXl,  (A6zz) RecvXu,
      (A6zz) SendZl,  (A6zz) SendZu,  (A6zz) RecvZl,  (A6zz) RecvZu,
      (A6zz) SendVl,  (A6zz) SendVu,  (A6zz) RecvVl,  (A6zz) RecvVu);

}



double Vlasov::getMaxNLTimeStep(const double maxCFL) 
{
  double dt_NL = 0.;
 
  #pragma omp single copyprivate(dt_NL)
  {
    // from non-linear ExB Term
    const double NL_ExB_v    = std::max(Xi_max[DIR_X]/dy , Xi_max[DIR_Y]/dx);

    // from non-linear Landau damping (particle trapping)
    const double NL_Landau_v = 0.; // Not included yet

    // if non-linear terms are very small, there is no restriction on
    // time step, add 1.e-10 to avoid division by 0 
    double NL = std::max(NL_ExB_v, NL_Landau_v) + 1.e-10; 

    // get global maximum timestep time step 
    dt_NL = maxCFL / parallel->reduce(NL, Op::max);
  }

  return dt_NL;
}


///////////////////////////////////////////   File I/O    ///////////////////////////////////////


void Vlasov::initData(Setup *setup, FileIO *fileIO) 
{
   
  //// Phasespace Group 
  hid_t psfGroup = fileIO->newGroup("/Vlasov", fileIO->getFileID());
  //set object attributes "Stores the phase space function");

  check(H5LTset_attribute_double(psfGroup, ".", "HyperViscosity", hyp_visc, 6), DMESG("Attribute"));
  
  
  // Phase space dimensions
  hsize_t psf_dim[]       = { Ns , Nm, Nz, Nky, Nx, Nv,             1 };
  hsize_t psf_maxdim[]    = { Ns , Nm, Nz, Nky, Nx, Nv, H5S_UNLIMITED };
  hsize_t psf_chunkBdim[] = { NsLB      , NmLB      , NzLB      , Nky , NxLB      , NvLB      ,             1 };
  hsize_t psf_chunkdim[]  = { NsLD      , NmLD      , NzLD      , Nky , NxLD      , NvLD      ,             1 };
  hsize_t psf_offset[]    = { NsLlB-1   , NmLlB-1   , NzLlB-1   , NkyLlB, NxLlB-1   , NvLlB-1   ,             0 };
  hsize_t psf_moffset[]   = { 0         , 0         , 2         , 0     , 2         , 2         , 0             };
     
  // ignore, as HDF-5 automatically (?) allocated data for it
/*
  FA_f0       = new FileAttr("f0", psfGroup, fileIO->file, 7, psf_dim, psf_maxdim, psf_chunkdim, psf_moffset,  psf_chunkBdim, psf_offset, true, fileIO->complex_tid);
  FA_f1       = new FileAttr("f1", psfGroup, fileIO->file, 7, psf_dim, psf_maxdim, psf_chunkdim, psf_moffset,  psf_chunkBdim, psf_offset, true, fileIO->complex_tid);
  FA_psfTime  = fileIO->newTiming(psfGroup);
  // call additional routines
  //initData(vlasovGroup, fileIO);  
 * */

  H5Gclose(psfGroup);

  if(fileIO->resumeFile == true) {

    // Currently we handle it here (should find better place)
    hid_t file_in = check(H5Fopen (fileIO->inputFileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), DMESG("H5Fopen"));

    // have to read from new file 
    FileAttr *FA_in_f0 = new FileAttr("/Vlasov/f0", file_in, fileIO->file, 7, psf_dim, psf_maxdim, psf_chunkdim, psf_moffset,  psf_chunkBdim, psf_offset, true, fileIO->complex_tid, false);
    FileAttr *FA_in_f1 = new FileAttr("/Vlasov/f1", file_in, fileIO->file, 7, psf_dim, psf_maxdim, psf_chunkdim, psf_moffset,  psf_chunkBdim, psf_offset, true, fileIO->complex_tid, false);

    FA_in_f0->read(ArrayPhase.data(f0));
    FA_in_f1->read(ArrayPhase.data(f ));
         
    delete FA_in_f0;
    delete FA_in_f1;

    H5Fclose(file_in); 

  }

  return;
}

void Vlasov::closeData() 
{
 // delete FA_f0;
 // delete FA_f1;
 // delete FA_psfTime;
}


void Vlasov::writeData(const Timing &timing, const double dt) 
{
  
  if (timing.check(dataOutputF1, dt)       )   {
      
    FA_f0->write(ArrayPhase.data(f0));
    FA_f1->write(ArrayPhase.data(f ));
    FA_psfTime->write(&timing);
      
    parallel->print("Wrote phase-space data ... "); 
  }
}


void Vlasov::loadData(FileIO *fileIO)
{
  // what is exactly the difference between reading and writting HDF-5 ? Setup should be same.
  // Phase space dimensions
  /*
  hsize_t psf_dim[]       = { grid->NsGD, grid->NmGD, grid->NvGD, grid->NzGD, grid->NkyGD, grid->NxGD,  1 };
  hsize_t psf_maxdim[]    = { grid->NsGD, grid->NmGD, grid->NvGD, grid->NzGD, grid->NkyGD, grid->NxGD, H5S_UNLIMITED};
  hsize_t psf_moffset[]   = { 0, 0, 2, 0, 2, 2, 0 };
  hsize_t psf_chunkBdim[] = { NsLB, NmLB, NvLB, NzLB, Nky, NxLB, 1};
  hsize_t psf_chunkdim[]  = { NsLD, NmLD, NvLD, NzLD, Nky, NxLD, 1};
  */
  return;
}

void Vlasov::printOn(std::ostream &output) const
{
  output << "Vlasov     | Type : " << equation_type <<  " Non-Linear : " << (doNonLinear ? "yes" : "no") << std::endl ;
  output << "Vlasov     | Hyperviscosity [ " ;
  for(int dir = DIR_X; dir <= DIR_S; dir++) output << hyp_visc[dir] << " ";
  output << " ] " << std::endl;
}

