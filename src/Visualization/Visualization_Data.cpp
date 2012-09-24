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
     hsize_t Fields_dim      [] = { numZSlide  , NkyLD, grid->NxGD,             1 };
     hsize_t Fields_maxdim   [] = { numZSlide  , NkyLD, grid->NxGD, H5S_UNLIMITED };
     hsize_t Fields_chunkBdim[] = { numZSlide  , NkyLD, NxLD      ,             1 };
     //hsize_t Fields_chunkBdim[] = { NzLB,  NkyLD, NxLB+4,        1};
     hsize_t Fields_chunkdim [] = { numZSlide  , NkyLD, NxLD      ,             1 };
     hsize_t Fields_moffset  [] = { 0, 0, 0, 0 };
     hsize_t Fields_offset   [] = {NzLlB-1, NkyLlB , NxLlB-1, 0 }; 
     // hsize_t Fields_moffset[]   = { 2, 0, 4, 0 };
     
     bool phiWrite = (parallel->Coord[DIR_VMS] == 0) && (Z[NzLlD] == 0.);
    
     hid_t visualGroup = check(H5Gcreate(fileIO->getFileID(), "/Visualization",H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), DMESG("Error creating group file for Phi : H5Gcreate"));
     FA_slphi      = new FileAttr("Phi", visualGroup, 4, Fields_dim, Fields_maxdim, Fields_chunkdim, Fields_moffset,  Fields_chunkBdim, Fields_offset, phiWrite && plasma->nfields >= 1,fileIO->complex_tid );
     FA_slAp       = new FileAttr("Ap",  visualGroup, 4, Fields_dim, Fields_maxdim, Fields_chunkdim, Fields_moffset,  Fields_chunkBdim, Fields_offset, phiWrite && plasma->nfields >= 2, fileIO->complex_tid);
     FA_slBp       = new FileAttr("Bp",  visualGroup, 4, Fields_dim, Fields_maxdim, Fields_chunkdim, Fields_moffset,  Fields_chunkBdim, Fields_offset, phiWrite && plasma->nfields >= 3 , fileIO->complex_tid);
     FA_slphiTime  = fileIO->newTiming(visualGroup);



      ////////////////////////// Other stuff (please cleanup) /////////////////////////////
//     FA_sln        = new FileAttr("Density",  visualGroup,     5, A4_dim, A4_maxdim, A4_chunkdim, A4_moffset, A4_chunkBdim, A4_offset, phiWrite);
//     FA_slT        = new FileAttr("Temperature",visualGroup,   5, A4_dim, A4_maxdim, A4_chunkdim, A4_moffset, A4_chunkBdim, A4_offset, phiWrite);
     
//     FA_slphiTime  = new FileAttr("Time", visualGroup,   1, time_dim, timing_maxdim, timing_chunkdim, offset0,  offset0, timing_chunkdim, parallel->myRank == 0, fileIO->timing_tid);
       

     // Visualization of Velocity space
     //hsize_t B4_offset[] = {NsLlB-1, 0, NzLlB-1, NyLlB-1, NxLlB-1, NvLlB-1, 0 }; 
     //FA_slF1       = new FileAttr("F1",  visualGroup,      7, B4_dim, B4_maxdim, B4_chunkdim, B4_moffset,  B4_chunkBdim, B4_offset, true);
     //hsize_t B4_offset[] = {NsLlB-1, 0, NzLlB-1, NyLlB-1, NxLlB-1, NvLlB-1, 0 }; 
     //FA_slF1       = new FileAttr("F1",  visualGroup,      4, V3_dim, V3_maxdim, V3_chunkdim, V3_moffset,  V3_chunkBdim, offset0, true, fileIO->complex_tid);
     //V3.resize(RxLD, RvLD, Range(NkyLlD, NkyLlD+nkyslides));

 //hsize_t Fields_offset[] = {NzLlB-1, NyLlB-1, NxLlB-1, 0 }; 
        
//     FA_slphiTime->create(visualGroup, "Time", offset0);
     /*  4D Stuff  XYZ + Species */
     int nslides=1;
     hsize_t A4_dim      [] =  { grid->NsGD, nslides, NkyLD, grid->NxGD,             1};
     hsize_t A4_maxdim   [] =  { grid->NsGD, nslides, NkyLD, grid->NxGD, H5S_UNLIMITED};
     hsize_t A4_chunkBdim[] = {       NsLD, NzLD   ,       NkyLD,        NxLD,        1};
     hsize_t A4_chunkdim [] =  { NsLD, nslides, NkyLD, NxLD, 1};
     hsize_t A4_moffset  [] = { 0, 0, 0, 0, 0 };
     hsize_t A4_offset[] = {NsLlB-1, NzLlB-1, NkyLlB, NxLlB-1, 0 }; 
     
     /////////////////////// Velocity space (X(0)-V) ///////////////////////////////////////
      visXV = setup->get("Visualization.XV", 0);
      if(visXV == true) {
        hsize_t XV_dim[]      =  { grid->NxGD, grid->NvGD,  grid->NsGD, 1  };
        hsize_t XV_maxdim[]   =  { grid->NxGD, grid->NvGD, grid->NsGD, H5S_UNLIMITED };
        hsize_t XV_chunkBdim[] = { grid->NxGD, NvLD      , NsLD, 1 };
        hsize_t XV_chunkdim[] =  { NxLD      , NvLD      , NsLD, 1 };
        hsize_t XV_moffset[]   = { 0, 0, 0, 0, 0, 0 };
     
        FA_XV  = new FileAttr("XV",  visualGroup,      4, XV_dim, XV_maxdim, XV_chunkdim, XV_moffset,  XV_chunkBdim, offset0, true, fileIO->complex_tid);
        //ArrXV.resize(RxLD, RvLD, RsLD);
     
      } 
     

     H5Gclose(visualGroup);
     

}   
    
Visualization_Data::~Visualization_Data() {

    delete FA_slphi;
    delete FA_slAp;
    delete FA_slBp;
    delete FA_slphiTime;
    
    if(visXV == true) delete FA_XV;

};

void Visualization_Data::writeData(Timing timing, double dt, bool force) 
{
        if (timing.check(dataOutputVisual,dt) || (timing.step == 1) || force) {
            FA_slphi->write(fields->ArrayField0.data(fields->Field0));
            FA_slphiTime->write(&timing);

      if(visXV == true) {
            //ArrXV(RxLD, RvLD, RsLD) = vlasov->f(RxLD, 0, NzLlD, RvLD, NmLlD, RsLD);
            //if(plasma->global == false) ArrXV(RxLD, RvLD, RsLD) += vlasov->f0(RxLD, 0, NzLlD, RvLD, NmLlD, RsLD);
            //FA_XV->write(ArrXV);

      }


      parallel->print("Wrote Visual  data ... "); 
      }
}
     
