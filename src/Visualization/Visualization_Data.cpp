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
    
    
     ////////////////////////// Visualization Data for potentials ///////////////////////////
     
     // We only save the first Z(z=0) slide, especially useful in 2D simulations
     int numZSlide = 1;
     hsize_t Fields_dim      [] = { numZSlide, NkyLD, Nx  ,             1 };
     hsize_t Fields_maxdim   [] = { numZSlide, NkyLD, Nx  , H5S_UNLIMITED };
     hsize_t Fields_chunkBdim[] = { numZSlide, NkyLD, NxLD,             1 };
     hsize_t Fields_chunkdim [] = { numZSlide, NkyLD, NxLD,             1 };
     hsize_t Fields_moffset  [] = { 0        , 0    ,   0 , 0             };
     hsize_t Fields_offset   [] = {NzLlB-1, NkyLlB , NxLlB-1, 0 }; 
     
     bool phiWrite = (parallel->Coord[DIR_VMS] == 0) && (parallel->Coord[DIR_Z] == 0);
    
     hid_t visualGroup = fileIO->newGroup("Visualization");

     FA_slphi      = new FileAttr("Phi", visualGroup, fileIO->file, 4, Fields_dim, Fields_maxdim, Fields_chunkdim, Fields_moffset,  Fields_chunkBdim, Fields_offset, phiWrite && plasma->nfields >= 1,fileIO->complex_tid );
     FA_slAp       = new FileAttr("Ap",  visualGroup, fileIO->file, 4, Fields_dim, Fields_maxdim, Fields_chunkdim, Fields_moffset,  Fields_chunkBdim, Fields_offset, phiWrite && plasma->nfields >= 2, fileIO->complex_tid);
     FA_slBp       = new FileAttr("Bp",  visualGroup, fileIO->file, 4, Fields_dim, Fields_maxdim, Fields_chunkdim, Fields_moffset,  Fields_chunkBdim, Fields_offset, phiWrite && plasma->nfields >= 3, fileIO->complex_tid);
     FA_slphiTime  = fileIO->newTiming(visualGroup);


     /////////////////////// Velocity space (X(0)-V) ///////////////////////////////////////
      visXV = setup->get("Visualization.XV", 0);
      if(visXV == true) {
        hsize_t XV_dim[]      =  { Nx, Nv,  Ns, 1  };
        hsize_t XV_maxdim[]   =  { Nx, Nv, Ns, H5S_UNLIMITED };
        hsize_t XV_chunkBdim[] = { Nx, NvLD      , NsLD, 1 };
        hsize_t XV_chunkdim[] =  { NxLD      , NvLD      , NsLD, 1 };
        hsize_t XV_moffset[]   = { 0, 0, 0, 0, 0, 0 };
     
        FA_XV  = new FileAttr("XV",  visualGroup, fileIO->file, 4, XV_dim, XV_maxdim, XV_chunkdim, XV_moffset,  XV_chunkBdim, offset0, true, fileIO->complex_tid);
        //ArrXV.resize(RxLD, RvLD, RsLD);
     
      } 
     

     H5Gclose(visualGroup);
     

}   
    
Visualization_Data::~Visualization_Data() 
{

  delete FA_slphi;
  delete FA_slAp;
  delete FA_slBp;
  delete FA_slphiTime;
    
  if(visXV == true) delete FA_XV;

};

void Visualization_Data::writeData(const Timing &timing, const double dt, const bool force) 
{
    
  if (timing.check(dataOutputVisual,dt) || force) {

    [=] (CComplex Field0[Nq][NzLD][NkyLD][NxLD]) {
    
      if(Nq >= 1) FA_slphi->write(&Field0[Field::phi][NzLlD][0][NxLlD]);
      if(Nq >= 2) FA_slAp ->write(&Field0[Field::Ap ][NzLlD][0][NxLlD]);
      if(Nq >= 3) FA_slBp ->write(&Field0[Field::Bp ][NzLlD][0][NxLlD]);
      
    }((A4zz) fields->Field0);
    
    FA_slphiTime->write(&timing);

    if(visXV == true) {
      
      //ArrXV(RxLD, RvLD, RsLD) = vlasov->f(RxLD, 0, NzLlD, RvLD, NmLlD, RsLD);
      //if(plasma->global == false) ArrXV(RxLD, RvLD, RsLD) += vlasov->f0(RxLD, 0, NzLlD, RvLD, NmLlD, RsLD);
      //FA_XV->write(ArrXV);
    
    }

    parallel->print("Wrote Visual  data ... "); 
    
  }
}
     
