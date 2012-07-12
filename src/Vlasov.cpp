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


Vlasov::Vlasov(Grid *_grid, Parallel *_parallel, Setup *_setup, FileIO *fileIO, Geometry<HELIOS_GEOMETRY> *_geo, FFTSolver *(_fft))    : fft(_fft),
boundary_isclean(true),   parallel(_parallel), grid(_grid), setup(_setup), geo(_geo),
f0(HeliosStorage), f(HeliosStorage), fs(HeliosStorage), fss(HeliosStorage),
ft(HeliosStorage), G(HeliosStorage4), Xi(HeliosStorage4),    f1(HeliosStorage)
{

   allocate(RxLB, RkyLD, RzLB, RvLB, RmLD, RsLD, f0, f, fss, fs, f1, ft);
   
   // for electro-static simulations we don"t need this but overhead
   allocate(RxLB, RkyLD, RzLB, RvLB, G, Xi);
   
   // allocate boundary (mpi) buffers
   
   allocate(RB  , RkyLD, RzLD, RvLD, RmLD, RsLD, SendXu, SendXl, RecvXu, RecvXl);
   allocate(RxLD, RB   , RzLD, RvLD, RmLD, RsLD, SendYu, SendYl, RecvYu, RecvYl);
   allocate(RxLD, RkyLD, RB  , RvLD, RmLD, RsLD, SendZu, SendZl, RecvZu, RecvZl);
   allocate(RxLD, RkyLD, RzLD, RB  , RmLD, RsLD, SendVu, SendVl, RecvVu, RecvVl);

   
   equation_type       = setup->get("Vlasov.Equation", "2D_ES");        
   calculate_nonLinear = setup->get("Vlasov.NonLinear", 0);
   useAntiAliasing     = setup->get("Vlasov.useAA", 0);

   
   std::string dir_string[] = { "X", "Y", "Z", "V", "M", "S" };
   
   for(int dir = DIR_X ; dir <= DIR_S ; dir++) hyper_visc[dir] = setup->get("Vlasov.HyperViscosity." + dir_string[dir]   ,  0.0);
       
   
   // needed for CFL condition
   Xi_max.resize(Range(DIR_X, DIR_Z)); Xi_max = 0.;
   collisions = new Collisions(grid, parallel, setup, fileIO, geo, fft);
    
   useBoundary     = setup->get("Vlasov.Boundary"    , "Periodic");
   krook_nu = setup->get("Vlasov.Krook.Nu", 0.);

}



Vlasov::~Vlasov() 
{

};


int Vlasov::solve(Fields *fields, Array6z  _fs, Array6z  _fss, double dt, int rk_step, int user_boundary_type)
{
   if(boundary_isclean == false) cleanBoundary(f_boundary);
   solve(equation_type, fields, _fs, _fss, dt, rk_step);

   // collision->addCollision(equation_type, fields, _fs, _fss, dt, rk_step);

   #pragma omp single
   setBoundary(_fss, user_boundary_type);       
   if(user_boundary_type == BOUNDARY_DIRTY) f_boundary.reference(_fss);
  
   return 1;
}

int Vlasov::cleanBoundary(Array6z A)
{
#ifdef GKC_PARALLEL_MPI
   // now wait until MPI communication finished and arrays are updated
   parallel->updateNeighboursBarrier();
    
   if(parallel->decomposition(DIR_X) > 1) {
    	A(Range(NxLlB  , NxLlD-1), RkyLD, RzLD, RvLD, RmLD, RsLD) = RecvXl(RB, RkyLD, RzLD, RvLD, RmLD, RsLD);
    	A(Range(NxLuD+1, NxLuB  ), RkyLD, RzLD, RvLD, RmLD, RsLD) = RecvXu(RB, RkyLD, RzLD, RvLD, RmLD, RsLD);
   }
     
   if((parallel->decomposition(DIR_Z) > 1) && (Nz > 1)) {
    	A(RxLD, RkyLD, Range(NzLlB  , NzLlB+1), RvLD, RmLD, RsLD) = RecvZl(RxLD, RkyLD, RB, RvLD, RmLD, RsLD);
    	A(RxLD, RkyLD, Range(NzLuD+1, NzLuB  ), RvLD, RmLD, RsLD) = RecvZu(RxLD, RkyLD, RB, RvLD, RmLD, RsLD);
   }

   if(parallel->decomposition(DIR_V) > 1) {
    	A(RxLD, RkyLD, RzLD, Range(NvLlB  , NvLlB+1), RmLD, RsLD) = RecvVl(RxLD, RkyLD, RzLD, RB, RmLD, RsLD);
    	A(RxLD, RkyLD, RzLD, Range(NvLuD+1, NvLuB  ), RmLD, RsLD) = RecvVu(RxLD, RkyLD, RzLD, RB, RmLD, RsLD); 
   }
#endif 
   boundary_isclean = true;

   return HELIOS_SUCCESS;
}

int Vlasov::setBoundary(Array6z  A , int boundary_type) {

   // still use periodicity but rescale f1/ use Krook instead
   /*
    if(useBoundary == "NoSlip") {
      const double sigma = 32.;
       for(int x = NxLlD; x <= NxLuD; x++) {
         A(x, RkyLD, RzLD, RvLD, RmLD, RsLD) *= 2./(1.+exp(- pow4((X(x)-Lx/2.)/sigma))+ exp(-pow4((X(x)+Lx/2.)/sigma)))-1.;
         
         A(Range(NxGlB, NxLlD), RkyLD, RzLD, RvLD, RmLD, RsLD) = 0.;
         A(Range(NxLuD, NxGuB), RkyLD, RzLD, RvLD, RmLD, RsLD) = 0.;
       }

    };
*/

   // X-Boundary
   if(parallel->decomposition(DIR_X) > 1) {
#ifdef GKC_PARALLEL_MPI
   // We do not need to communicate for M and S as we do not have boundary cells
   SendXl(RB, RkyLD, RzLD, RvLD, RmLD, RsLD) = A(Range(NxLlD  , NxLlD+1), RkyLD, RzLD, RvLD, RmLD, RsLD);
   SendXu(RB, RkyLD, RzLD, RvLD, RmLD, RsLD) = A(Range(NxLuD-1, NxLuD  ), RkyLD, RzLD, RvLD, RmLD, RsLD);
   parallel->updateNeighbours(SendXu, SendXl, RecvXu, RecvXl, DIR_X);
#else
   check(-1, DMESG("decompositionn in X, but compiled without MPI support"));
#endif 
  } else {
   A(Range(NxLuD+1, NxLuB)  , RkyLD, RzLD, RvLD, RmLD, RsLD) = A(Range(NxLlD  , NxLlD+1), RkyLD, RzLD, RvLD, RmLD, RsLD);
   A(Range(NxLlB  , NxLlD-1), RkyLD, RzLD, RvLD, RmLD, RsLD) = A(Range(NxLuD-1, NxLuD)  , RkyLD, RzLD, RvLD, RmLD, RsLD);
  }

   // Z-Boundary 
   if((parallel->decomposition(DIR_Z) > 1) && (Nz > 1)) {
#ifdef GKC_PARALLEL_MPI
	
   check(-1, DMESG("Boundary not implemented"));

   /*
        for(int x=NxLlD; x<= NxLuD;x++) { for(int y=NyLlD; y<= NyLuD;y++) {
            ShearB b = geo->getYPos(x,y);
      // original
           // SendZl(x, y, RB, RvLD, RmLD, RsLD) = A(x, b.ly, Range(NzLlD  , NzLlD+1), RvLD, RmLD, RsLD);
           // SendZu(x, y, RB, RvLD, RmLD, RsLD) = A(x, b.uy, Range(NzLuD-1, NzLuD  ), RvLD, RmLD, RsLD);
            
            SendZl(x, b.ly, RB, RvLD, RmLD, RsLD) = A(x, y, Range(NzLlD  , NzLlD+1), RvLD, RmLD, RsLD);
            SendZu(x, y, RB, RvLD, RmLD, RsLD) = A(x, b.ly, Range(NzLuD-1, NzLuD  ), RvLD, RmLD, RsLD);

            check(-1, DMESG("NOT FIXED"));
        }}
   */

   parallel->updateNeighbours(SendZu, SendZl, RecvZu, RecvZl, DIR_Z);
#else
   check(-1, DMESG("decompositionn in Y, but compiled without MPI support"));
#endif
   } else if (Nz > 1){ 
	check(-1, DMESG("Boundary not implemented"));
   /*
     // For z-we need to connect the magnetic field lines
     for(int x=NxLlD; x<= NxLuD;x++) { for(int y=NyLlD; y<= NyLuD;y++) {
         ShearB b = geo->getYPos(x,y);
         // orginal
         //A(x, y, Range(NzLuD+1, NzLuB  ), RvLD, RmLD, RsLD) = A(x, b.ly, Range(NzLlD  , NzLlD+1), RvLD, RmLD, RsLD);
         //A(x, y, Range(NzLlB,   NzLlB+1), RvLD, RmLD, RsLD) = A(x, b.uy, Range(NzLuD-1, NzLuD  ), RvLD, RmLD, RsLD);
        
         // test
         A(x, b.ly, Range(NzLuD+1, NzLuB  ), RvLD, RmLD, RsLD) = A(x, y   , Range(NzLlD  , NzLlD+1), RvLD, RmLD, RsLD);
         A(x,    y, Range(NzLlB,   NzLlB+1), RvLD, RmLD, RsLD) = A(x, b.ly, Range(NzLuD-1, NzLuD  ), RvLD, RmLD, RsLD);
     }}
  */
  } else ; // for Nz = 1 (2D simulations) no z boundary is required

  /*
        A(x, y,    Range(NzLlB, NzLlB+1), RvLD, RmLD, Range::all())    = b.y0 * A(x, b.ypos, Range(NzLuD-1, NzLuD), RvLD, RmLD, Range::all()) + (1. - b.y0) *  A(x, b.ypos-1, Range(NzLuD-1, NzLuD), RvLD, RmLD, Range::all());
        A(x, b.ypos, Range( NzLuD+1, NzLuB), RvLD, RmLD, Range::all()) = b.y0 * A(x, y, Range(NzLlD, NzLlD+1), RvLD, RmLD, Range::all()) + (1. - b.y0) * A(x, y+1, Range(NzLlD, NzLlD+1), RvLD, RmLD, Range::all());
*/  

   if(parallel->decomposition(DIR_V) > 1) {
#ifdef GKC_PARALLEL_MPI
      SendVl(RxLD, RkyLD,  RzLD, RB, RmLD, RsLD) = A(RxLD, RkyLD, RzLD, Range(NvLlD  , NvLlD+1), RmLD, RsLD);
      SendVu(RxLD, RkyLD,  RzLD, RB, RmLD, RsLD) = A(RxLD, RkyLD, RzLD, Range(NvLuD-1, NvLuD  ), RmLD, RsLD);
      parallel->updateNeighbours(SendVu, SendVl, RecvVu, RecvVl, DIR_V);
#else
      check(-1, DMESG("decompositionn in Y, but compiled without MPI support"));
#endif 
   } else {
      A(RxLD, RkyLD, RzLD, Range(NvGlB  , NvGlB+1), RmLD, RsLD) = 0.e0;
      A(RxLD, RkyLD, RzLD, Range(NvGuD+1, NvGuB  ), RmLD, RsLD) = 0.e0;
   }

   // For the Vlasov equation non-blocking IO is possible (Poisson eq. does not need boundary values)
   if      (boundary_type == BOUNDARY_CLEAN) cleanBoundary(A);
   else if (boundary_type == BOUNDARY_DIRTY) boundary_isclean = false;
   else    check(-1, DMESG("Vlasov : Only clean/dirty values are allowed"));

   return HELIOS_SUCCESS;

}
        
void Vlasov::printOn(ostream &output) const
{
   output << "Vlasov     | Type : " << equation_type <<  " Non-Linear : " << (calculate_nonLinear ? "yes" : "no") << std::endl ;
   output << "Vlasov     | Hyperviscoisty [ " ;
   for(int dir = DIR_X ; dir <= DIR_S ; dir++) output << hyper_visc[dir] << " ";
   output << " ] " << std::endl;
};

double Vlasov::getMaxTimeStep(int dir, const double maxCFL) 
{

  double v_scale = 0., dt=0.;
  
  for(int s=NsGlD; s<=NsGuD; s++) v_scale = max(v_scale, plasma->species(s).scale_v);
	
  if     (dir == DIR_X  ) dt =  maxCFL / max(1.e-99, parallel->collect(Xi_max(DIR_X)/dx, OP_MAX));
  else if(dir == DIR_Y  ) dt =  maxCFL / max(1.e-99, parallel->collect(Xi_max(DIR_Y)/dy, OP_MAX));
  if     (dir == DIR_XY ) dt =  maxCFL / max(1.e-99, parallel->collect(Xi_max(DIR_X)/dy + Xi_max(DIR_Y)/dx, OP_MAX));
  else if(dir == DIR_Z  ) dt =  maxCFL / max(1.e-99, parallel->collect(Xi_max(DIR_Z)/dz, OP_MAX));
  else if(dir == DIR_V  ) dt =  maxCFL / (v_scale * Lv/(sqrt(geo->eps_hat) * dz));
  else if(dir == DIR_ALL) dt =  maxCFL / parallel->collect(Xi_max(DIR_Y)/dx + Xi_max(DIR_X)/dy +  
                                              v_scale*Lv*Xi_max(DIR_Z)/(sqrt(geo->eps_hat)*dz) +  v_scale * Lv/(sqrt(geo->eps_hat)*dz), OP_MAX);
  return dt;
}


///////////////////////////////////   File I/O ////////////////////


// BUG : Not working (HDF-5 intialized the whole size even if we don't write
void Vlasov::initDataOutput(FileIO *fileIO) {

   
  //////////////////////////////////////////////////////////////// Phasespace Group ////////////////////////////////////////////////////////
   
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
        
void Vlasov::updateCFL(const cmplxd dphi_dx, const cmplxd dphi_dy, const cmplxd dphi_dz)
{

  Xi_max(DIR_X) = max(Xi_max(DIR_X), abs(dphi_dx));
  Xi_max(DIR_Y) = max(Xi_max(DIR_Y), abs(dphi_dy));
  Xi_max(DIR_Z) = max(Xi_max(DIR_Z), abs(dphi_dz));

};
 
