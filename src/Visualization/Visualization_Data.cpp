/*
 * =====================================================================================
 *
 *       Filename: Visualization_Data.h
 *
 *    Description: Writes only slides of potential to disk.
 *
 *         Author: Paul P. Hilscher (2011-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#include "Visualization/Visualization_Data.h"


Visualization_Data::Visualization_Data(Grid *grid, Parallel *_parallel, Setup *setup, FileIO *fileIO, Vlasov *_vlasov, Fields *_fields) 
: Visualization(_vlasov, _fields, grid, setup, _parallel) 
{

  // Define some common use variables
  hsize_t offset0[] = { 0, 0, 0, 0, 0,  0, 0};
    
  hsize_t timing_chunkdim[] = {1}; hsize_t 
  timing_maxdim          [] = {H5S_UNLIMITED};
  hsize_t time_dim       [] = { 1 }; 
    
  hid_t visGroup = fileIO->newGroup("Visualization");
    
  ////////////////////////// Visualization Data for potentials ///////////////////////////
  {   
    // We only save the first Z(z=0) slide, especially useful in 2D simulations
    int numZSlide = 1;
    hsize_t dim  [] = { numZSlide, Nky, Nx     ,             1 };
    hsize_t mdim [] = { numZSlide, Nky, Nx     , H5S_UNLIMITED };
    hsize_t cBdim[] = { numZSlide, Nky, NxLD   ,             1 };
    hsize_t cdim [] = { numZSlide, Nky, NxLD   ,             1 };
    hsize_t moff [] = { 0        , 0  , 0      ,             0 };
    hsize_t off  [] = { NzLlB-1  , 0  , NxLlB-1,             0 }; 
     
    bool write = (parallel->Coord[DIR_VMS] == 0) && (parallel->Coord[DIR_Z] == 0);
    
    FA_slphi  = new FileAttr("Phi", visGroup, fileIO->file, 4, dim, mdim, cdim, moff, cBdim, off, write && Nq >= 1, fileIO->complex_tid);
    FA_slAp   = new FileAttr("Ap" , visGroup, fileIO->file, 4, dim, mdim, cdim, moff, cBdim, off, write && Nq >= 2, fileIO->complex_tid);
    FA_slBp   = new FileAttr("Bp" , visGroup, fileIO->file, 4, dim, mdim, cdim, moff, cBdim, off, write && Nq >= 3, fileIO->complex_tid);
    FA_slTime = fileIO->newTiming(visGroup);

  }
  /////////////////////// Velocity space (X(0)-V) ///////////////////////////////////////
  
  visXV = setup->get("Visualization.XV", 0);
  
  if(visXV == true) {
  
    hsize_t dim  [] = { Nx  , Nv  , Ns  ,             1 };
    hsize_t mdim [] = { Nx  , Nv  , Ns  , H5S_UNLIMITED };
    hsize_t cBdim[] = { Nx  , NvLD, NsLD,             1 };
    hsize_t cdim [] = { NxLD, NvLD, NsLD,             1 };
    hsize_t moff [] = { 0   , 0   , 0   ,             0 };
     
    FA_XV  = new FileAttr("XV",  visGroup, fileIO->file, 4, dim, mdim, cdim, moff,  cBdim, offset0, true, fileIO->complex_tid);
  
  } 
     
  H5Gclose(visGroup);

}   
    
Visualization_Data::~Visualization_Data() 
{

  delete FA_slphi;
  delete FA_slAp;
  delete FA_slBp;
  delete FA_slTime;
    
  if(visXV == true) delete FA_XV;

}

void Visualization_Data::writeData(const Timing &timing, const double dt, const bool force) 
{
    
  if (timing.check(dataOutputVisual,dt) || force) {

    [=](CComplex Field0[Nq][NzLD][Nky][NxLD]) {
    
      if(Nq >= 1) FA_slphi->write(&Field0[Field::phi][NzLlD][0][NxLlD]);
      if(Nq >= 2) FA_slAp ->write(&Field0[Field::Ap ][NzLlD][0][NxLlD]);
      if(Nq >= 3) FA_slBp ->write(&Field0[Field::Bp ][NzLlD][0][NxLlD]);
      
    }((A4zz) fields->Field0);
    
    FA_slTime->write(&timing);

    if(visXV == true) {
      
      //ArrXV(RxLD, RvLD, RsLD) = vlasov->f(RxLD, 0, NzLlD, RvLD, NmLlD, RsLD);
      //if(plasma->global == false) ArrXV(RxLD, RvLD, RsLD) += vlasov->f0(RxLD, 0, NzLlD, RvLD, NmLlD, RsLD);
      //FA_XV->write(ArrXV);
    
    }

    parallel->print("Wrote Visual  data ... "); 
    
  }
}
     
