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


Vlasov::Vlasov(Grid *_grid, Parallel *_parallel, Setup *_setup, FileIO *fileIO, Geometry *_geo, FFTSolver *_fft, Benchmark *_bench)
: fft(_fft), bench(_bench),
boundary_isclean(true),   parallel(_parallel), grid(_grid), setup(_setup), geo(_geo),
f0(GKCStorage), f(GKCStorage), fs(GKCStorage), fss(GKCStorage),
ft(GKCStorage), G(GKCStorage4), Xi(GKCStorage4),    f1(GKCStorage), nonLinearTerms(GKCStorage3)
{

   allocate(RxLB, RkyLD, RzLB, RvLB, RmLD, RsLD, f0, f, fss, fs, f1, ft);
   
   // for electro-static simulations we don"t need this but overhead
   allocate(RxLB , RkyLD, RzLB, RvLB, G);
   allocate(RxLB4, RkyLD, RzLB, RvLB, Xi);
   allocate(RxLD , RkyLD, RvLD, nonLinearTerms);
   
   // allocate boundary (mpi) buffers
   allocate(RB  , RkyLD, RzLD, RvLD, RmLD, RsLD, SendXu, SendXl, RecvXu, RecvXl);
   allocate(RxLD, RB   , RzLD, RvLD, RmLD, RsLD, SendYu, SendYl, RecvYu, RecvYl);
   allocate(RxLD, RkyLD, RB  , RvLD, RmLD, RsLD, SendZu, SendZl, RecvZu, RecvZl);
   allocate(RxLD, RkyLD, RzLD, RB  , RmLD, RsLD, SendVu, SendVl, RecvVu, RecvVl);

   equation_type       = setup->get("Vlasov.Equation" , "2D_ES");        
   nonLinear           = setup->get("Vlasov.NonLinear", 0      );

   std::string dir_string[] = { "X", "Y", "Z", "V", "M", "S" };
   

   for(int dir = DIR_X ; dir <= DIR_S ; dir++) hyper_visc[dir] = setup->get("Vlasov.HyperViscosity." + dir_string[dir]   ,  0.0);
       
   
   collisions = new Collisions(grid, parallel, setup, fileIO, geo, fft);
    
}



Vlasov::~Vlasov() 
{

};


int Vlasov::solve(Fields *fields, Array6C  _fs, Array6C  _fss, double dt, int rk_step, const double rk[3], bool useNonBlockingBoundary)
{

   // use static function
   useNonBlockingBoundary = false; 
   // Need boundary_isclean to avoid deadlock at first iteration
   if((boundary_isclean == false) && useNonBlockingBoundary) setBoundary(f_boundary, BOUNDARY_RECV);
  
   Xi_max[:] = 0.; // Needed to calculate CFL time step 
   
   solve(equation_type, fields, _fs, _fss, dt, rk_step, rk);


   // Note : we have non-blocking boundaries as Poisson solver does not require ghosts
   (useNonBlockingBoundary) ? setBoundary(_fss, BOUNDARY_SEND) :  setBoundary(_fss, BOUNDARY_SENDRECV); 

   f_boundary.reference(_fss); boundary_isclean = false;
  
   return 1;
}

void Vlasov::setBoundary(Array6C A) 
{ 
  setBoundary(A, BOUNDARY_SENDRECV); 
};


void Vlasov::setBoundary(Array6C A, int boundary_type)
{


    [=] (
         CComplex g     [NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB],
         CComplex SendXl[NsLD][NmLD][NzLD][NkyLD][GC2 ][NvLD], CComplex SendXu[NsLD][NmLD][NzLD][NkyLD][GC2 ][NvLD], 
         CComplex RecvXl[NsLD][NmLD][NzLD][NkyLD][GC2 ][NvLD], CComplex RecvXu[NsLD][NmLD][NzLD][NkyLD][GC2 ][NvLD], 
         CComplex SendZl[NsLD][NmLD][GC2 ][NkyLD][NxLD][NvLD], CComplex SendZu[NsLD][NmLD][GC2 ][NkyLD][NxLD][NvLD],
         CComplex RecvZl[NsLD][NmLD][GC2 ][NkyLD][NxLD][NvLD], CComplex RecvZu[NsLD][NmLD][GC2 ][NkyLD][NxLD][NvLD],
         CComplex SendVl[NsLD][NmLD][NzLD][NkyLD][NxLD][GC2 ], CComplex SendVu[NsLD][NmLD][NvLD][NkyLD][NxLD][GC2 ],
         CComplex RecvVl[NsLD][NmLD][NzLD][NkyLD][NxLD][GC2 ], CComplex RecvVu[NsLD][NmLD][NvLD][NkyLD][NxLD][GC2 ],
         int boundary_type)
{

  /////////////////////////// Send Boundaries //////////////////////////////
  if(!(boundary_type & BOUNDARY_SEND)) {


   // X-Boundary (Note, we may have different boundaries for global simulations)
   SendXl[:][:][:][:][:][:] = g[NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][NkyLlD:NkyLD][NxLlD  :2][NvLlD:NvLD];
   SendXu[:][:][:][:][:][:] = g[NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][NkyLlD:NkyLD][NxLuD-1:2][NvLlD:NvLD];
   parallel->updateNeighbours(Vlasov::SendZu, Vlasov::SendZl, Vlasov::RecvZu,Vlasov::RecvZl, DIR_Z);
  
   // We do not domain decompose poloidal (y) fourier modes, thus boundaries not required
  
   // Z-Boundary 
   if(Nz > 1) for(int x=NxLlD; x<= NxLuD;x++) { omp_for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) {
           
            //const Complex a = newComplex(0., 2. * M_PI * (2.* M_PI/Ly) * y_k);
            
            SendZl[:][:][:][:][y_k][x-3] = g[NsLlD:NsLD][NmLlD:NmLD][NzLlD  :2][y_k][x][NvLlD:NvLD]; //* exp( ((NzLlD == NzGlD) ? a : 0. ) *geo->nu(x));
            SendZu[:][:][:][:][y_k][x-3] = g[NsLlD:NsLD][NmLlD:NmLD][NzLuD-1:2][y_k][x][NvLlD:NvLD]; //* (exp(-((Complex) ((NzLlD == NzGlD) ? a : 0.) * geo->nu(x)))));
               
  } }
  parallel->updateNeighbours(Vlasov::SendZu, Vlasov::SendZl, Vlasov::RecvZu, Vlasov::RecvZl, DIR_Z);
  
  // We do not need to communicate for M and S as we do not have boundary cells (yet)
  
  // Decomposition in velocity is rather unlikely, thus give non-decomposed version too
  // We set endpoint in velocity space to zero
  if(parallel->decomposition[DIR_V] > 1) {
       SendVl[:][:][:][:][:][:] = g[NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][NkyLlD:NkyLD][NxLlD:NxLD][NvLlD  :2]; 
       SendVu[:][:][:][:][:][:] = g[NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][NkyLlD:NkyLD][NxLlD:NxLD][NvLuD-1:2]; 
       parallel->updateNeighbours(Vlasov::SendVu, Vlasov::SendVl, Vlasov::RecvVu, Vlasov::RecvVl, DIR_V);
  } else {
       g[NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][NkyLlD:NkyLD][NxLlD:NxLD][NvLlD  :2] = 0.;
       g[NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][NkyLlD:NkyLD][NxLlD:NxLD][NvLuD-1:2] = 0.;
  }


  }

  /////////////////////////// Receive Boundaries //////////////////////////////
  else if(!(boundary_type & BOUNDARY_RECV))
  {
 
    // Get boundaries
     parallel->updateNeighboursBarrier();
   
   // Set boundary in X (take care of Neumann boundary ?!) 
   g[NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][NkyLlD:NkyLD][NxLlB  :2][NvLlD:NvLD] = RecvXl[:][:][:][:][:][:]; 
   g[NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][NkyLlD:NkyLD][NxLuD+1:2][NvLlD:NvLD] = RecvXu[:][:][:][:][:][:]; 
    
   // Set boundary in Z
   if(Nz > 1) {
   
       g[NsLlD:NsLD][NmLlD:NmLD][NzLlB  :2][NkyLlD:NkyLD][NxLlD:NxLD][NvLlD:NvLD] = RecvZl[:][:][:][:][:][:]; 
       g[NsLlD:NsLD][NmLlD:NmLD][NzLuD+1:2][NkyLlD:NkyLD][NxLlD:NxLD][NvLlD:NvLD] = RecvZu[:][:][:][:][:][:]; 
   }

   // Set boundary in V
   if(parallel->decomposition[DIR_V] > 1) {
       
       g[NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][NkyLlD:NkyLD][NxLlD:NxLD][NvLlB  :2] = RecvVl[:][:][:][:][:][:]; 
       g[NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][NkyLlD:NkyLD][NxLlD:NxLD][NvLuD+1:2] = RecvVu[:][:][:][:][:][:]; 

   }
   
   boundary_isclean = true;
  
  }

}  ((A6zz) A.dataZero(),  (A6zz) SendXl.data(),  (A6zz) SendXu.data(),  (A6zz) RecvXl.data(),  (A6zz) RecvXu.data(),
                   (A6zz) SendZl.data(),  (A6zz) SendZu.data(),  (A6zz) RecvZl.data(),  (A6zz) RecvZu.data(),
                   (A6zz) SendVl.data(),  (A6zz) SendVu.data(),  (A6zz) RecvVl.data(),  (A6zz) RecvVu.data(),
     boundary_type);



}



double Vlasov::getMaxTimeStep(int dir, const double maxCFL) 
{
 // simplify !
  double v_scale = 0., dt=0.;
  
  for(int s=NsGlD; s<=NsGuD; s++) v_scale = max(v_scale, plasma->species(s).scale_v);
   
  if     (dir == DIR_X  ) dt =  maxCFL / max(1.e-99, parallel->collect(Xi_max[DIR_X]/dx, OP_MAX));
  else if(dir == DIR_Y  ) dt =  maxCFL / max(1.e-99, parallel->collect(Xi_max[DIR_Y]/dy, OP_MAX));
  if     (dir == DIR_XY ) dt =  maxCFL / max(1.e-99, parallel->collect(Xi_max[DIR_X]/dy + Xi_max[DIR_Y]/dx, OP_MAX));
  else if(dir == DIR_Z  ) dt =  maxCFL / max(1.e-99, parallel->collect(Xi_max[DIR_Z]/dz, OP_MAX));
  else if(dir == DIR_V  ) dt =  maxCFL / (v_scale * Lv/(sqrt(geo->eps_hat) * dz));
  else if(dir == DIR_ALL) dt =  maxCFL / parallel->collect(Xi_max[DIR_Y]/dx + Xi_max[DIR_X]/dy +  
                                              v_scale*Lv*Xi_max[DIR_Z]/(sqrt(geo->eps_hat)*dz) +  v_scale * Lv/(sqrt(geo->eps_hat)*dz), OP_MAX);
  return dt;
}


///////////////////////////////////////////   File I/O    ///////////////////////////////////////


// BUG : Not working (HDF-5 intialized the whole size even if we don't write
void Vlasov::initDataOutput(FileIO *fileIO) {

   
  //// Phasespace Group 
   
  hid_t psfGroup = check(H5Gcreate(fileIO->getFileID(), "/Vlasov",H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), DMESG("Error creating group file for Phasespace : H5Gcreate"));
  
  /*
   
   hsize_t psf_offset[7] =  { NsLlB-1, NmLlB-1, NvLlB-1, NzLlB-1, NyLlB-1, NxLlB-1, 0 };
   
   // FA_psf->create(psfGroup, "Data", psf_offset);
   // FA_psfTime->create(psfGroup, "Timing", offset0);
        
   hsize_t psf_dim[]       = { grid->NsGD, grid->NmGD, grid->NvGD, grid->NzGD, grid->NyGD, grid->NxGD,  1 };
   hsize_t psf_maxdim[]    = { grid->NsGD, grid->NmGD, grid->NvGD, grid->NzGD, grid->NyGD, grid->NxGD, H5S_UNLIMITED};
   hsize_t psf_moffset[]   = { 0, 0, 2, 2, 2, 2, 0 };
   hsize_t psf_chunkBdim[] = { grid->NsGD, grid->NmGD, grid->NvLB, grid->NzLB,  grid->NyLB, grid->NxLB, 1};
   hsize_t psf_chunkdim[]  = {NsLD, NmLD, NvLD, NzLD, NyLD, NxLD, 1};
     
   // FA_psf      = new FileAttr("Unnamed",7, psf_dim, psf_maxdim, psf_chunkdim, psf_moffset,  psf_chunkBdim, true);
   // FA_psfTime  = new FileAttr("Unnamed",1, time_dim, timing_maxdim, timing_chunkdim, offset0,  timing_chunkdim, parallel->myRank == 0, timing_tid);
  */
          
  H5Gclose(psfGroup);
}

/*
// Fix it ! 
 // update only from non-linear term (df1_dv)
void Vlasov::updateCFL(const Complex dphi_dx, const Complex dphi_dy, const Complex dphi_dz)
{
  Xi_max(DIR_X) = max(Xi_max(DIR_X), abs(dphi_dx));
  Xi_max(DIR_Y) = max(Xi_max(DIR_Y), abs(dphi_dy));
  Xi_max(DIR_Z) = max(Xi_max(DIR_Z), abs(dphi_dz));
};
 
 * */







void Vlasov::printOn(ostream &output) const
{
   output << "Vlasov     | Type : " << equation_type <<  " Non-Linear : " << (nonLinear ? "yes" : "no") << std::endl ;
   output << "Vlasov     | Hyperviscosity [ " ;
   for(int dir = DIR_X ; dir <= DIR_S ; dir++) output << hyper_visc[dir] << " ";
   output << " ] " << std::endl;
};


