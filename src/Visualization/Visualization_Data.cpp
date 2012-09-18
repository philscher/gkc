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

// ToDO : Use Fiedls instead of phi , Ap, Bp


#include "Visualization/Visualization_Data.h"

    

Visualization_Data::Visualization_Data(Grid *grid, Parallel *_parallel, Setup *setup, FileIO *fileIO, Vlasov *_vlasov, Fields *_fields) 
: Visualization(_vlasov, _fields, grid, setup, _parallel) 
{

   // we wanna save only slides
     int nslides = 1;
     /*  3D Stuff */
     hsize_t A3_dim[]      =  { nslides, NkyLD, grid->NxGD,             1};
     hsize_t A3_maxdim[]   =  { nslides, NkyLD, grid->NxGD, H5S_UNLIMITED};
     hsize_t A3_chunkBdim[] = { NzLB,  NkyLD, NxLB+4,        1};
     hsize_t A3_chunkdim[] =  { nslides, NkyLD, NxLD, 1};
     hsize_t A3_moffset[]   = { 2, 0, 4, 0 };
     
     /*  4D Stuff  XYZ + Species */
     hsize_t A4_dim[]      =  { grid->NsGD, nslides, NkyLD, grid->NxGD,             1};
     hsize_t A4_maxdim[]   =  { grid->NsGD, nslides, NkyLD, grid->NxGD, H5S_UNLIMITED};
     hsize_t A4_chunkBdim[] = {       NsLD, NzLD   ,       NkyLD,        NxLD,        1};
     hsize_t A4_chunkdim[] =  { NsLD, nslides, NkyLD, NxLD, 1};
     hsize_t A4_moffset[]   = { 0, 0, 0, 0, 0 };
     
     /*  4D Stuff  XkyZV + Species */
     //hsize_t B4_dim[]      =  { grid->NsGD, 1, nslides, grid->NyGD, grid->NxGD, grid->NvGD ,           1};
     //hsize_t B4_maxdim[]   =  { grid->NsGD, 1, nslides, grid->NyGD, grid->NxGD, grid->NvGD , H5S_UNLIMITED};
     //hsize_t B4_chunkBdim[] = {       NsLD, 1, NzLD   ,       NyLD,        NxLD, NvLD,       1};
     //hsize_t B4_chunkdim[] =  { NsLD, 1, NyLD, NxLD, NvLD, 1};
     // hsize_t B4_moffset[]   = { 0, 0, 0, 0, 0, 0 };
     int nkyslides =1; 
     hsize_t V3_dim[]      =  { grid->NxGD, grid->NvGD, nkyslides, 1  };
     hsize_t V3_maxdim[]   =  { grid->NxGD, grid->NvGD, nkyslides , H5S_UNLIMITED };
     hsize_t V3_chunkBdim[] = { grid->NxGD, NvLD      , nkyslides , 1 };
     hsize_t V3_chunkdim[] =  { NxLD      , NvLD      , nkyslides , 1 };
     hsize_t V3_moffset[]   = { 0, 0, 0, 0, 0, 0 };
     
     
     
     
     
     
     bool phiWrite = (parallel->Coord[DIR_VMS] == 0) && (Z[NzLlD] == 0.);

     hsize_t offset0[] = { 0, 0, 0, 0, 0,  0, 0};
    
     hsize_t timing_chunkdim[] = {1}; hsize_t timing_maxdim[] = {H5S_UNLIMITED};
     hsize_t time_dim[1] = { 1 }; 
    
     hsize_t A3_offset[] = {NzLlB-1, NkyLlB, NxLlB-1, 0 }; 
     hsize_t A4_offset[] = {NsLlB-1, NzLlB-1, NkyLlB, NxLlB-1, 0 }; 
    
     hid_t visualGroup = check(H5Gcreate(fileIO->getFileID(), "/Visualization",H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), DMESG("Error creating group file for Phi : H5Gcreate"));
     FA_slphi      = new FileAttr("Phi", visualGroup,       4, A3_dim, A3_maxdim, A3_chunkdim, A3_moffset,  A3_chunkBdim, A3_offset, phiWrite && plasma->nfields >= 1,fileIO->complex_tid );
     FA_slAp       = new FileAttr("Ap",  visualGroup,      4, A3_dim, A3_maxdim, A3_chunkdim, A3_moffset,  A3_chunkBdim, A3_offset, phiWrite && plasma->nfields >= 2, fileIO->complex_tid);
     FA_slBp       = new FileAttr("Bp",  visualGroup,      4, A3_dim, A3_maxdim, A3_chunkdim, A3_moffset,  A3_chunkBdim, A3_offset, phiWrite && plasma->nfields >= 3 , fileIO->complex_tid);

//     FA_sln        = new FileAttr("Density",  visualGroup,     5, A4_dim, A4_maxdim, A4_chunkdim, A4_moffset, A4_chunkBdim, A4_offset, phiWrite);
//     FA_slT        = new FileAttr("Temperature",visualGroup,   5, A4_dim, A4_maxdim, A4_chunkdim, A4_moffset, A4_chunkBdim, A4_offset, phiWrite);
     
//     FA_slphiTime  = new FileAttr("Time", visualGroup,   1, time_dim, timing_maxdim, timing_chunkdim, offset0,  offset0, timing_chunkdim, parallel->myRank == 0, fileIO->timing_tid);
     FA_slphiTime  = fileIO->newTiming(visualGroup);
       

     // Visualization of Velocity space
     //hsize_t B4_offset[] = {NsLlB-1, 0, NzLlB-1, NyLlB-1, NxLlB-1, NvLlB-1, 0 }; 
     //FA_slF1       = new FileAttr("F1",  visualGroup,      7, B4_dim, B4_maxdim, B4_chunkdim, B4_moffset,  B4_chunkBdim, B4_offset, true);
     //hsize_t B4_offset[] = {NsLlB-1, 0, NzLlB-1, NyLlB-1, NxLlB-1, NvLlB-1, 0 }; 
     //FA_slF1       = new FileAttr("F1",  visualGroup,      4, V3_dim, V3_maxdim, V3_chunkdim, V3_moffset,  V3_chunkBdim, offset0, true, fileIO->complex_tid);
     //V3.resize(RxLD, RvLD, Range(NkyLlD, NkyLlD+nkyslides));

 //hsize_t A3_offset[] = {NzLlB-1, NyLlB-1, NxLlB-1, 0 }; 
        
//     FA_slphi->create(visualGroup,"phi", A3_offset);
//     FA_slAp->create(visualGroup, "Ap", A3_offset);
//     FA_slBp->create(visualGroup, "Bp", A3_offset);

//     FA_sln->create(visualGroup, "n", A4_offset);
     //FA_sldn->create(visualGroup, "dn", A3_offset);
//     FA_slT->create(visualGroup, "T", A4_offset);

//     FA_slphiTime->create(visualGroup, "Time", offset0);
   
        // Visualize X-V
      visXV = setup->get("Visualization.XV", 0);
      if(visXV == true) {
        hsize_t XV_dim[]      =  { grid->NxGD, grid->NvGD,  grid->NsGD, 1  };
        hsize_t XV_maxdim[]   =  { grid->NxGD, grid->NvGD, grid->NsGD, H5S_UNLIMITED };
        hsize_t XV_chunkBdim[] = { grid->NxGD, NvLD      , NsLD, 1 };
        hsize_t XV_chunkdim[] =  { NxLD      , NvLD      , NsLD, 1 };
        hsize_t XV_moffset[]   = { 0, 0, 0, 0, 0, 0 };
     
        FA_XV  = new FileAttr("XV",  visualGroup,      4, XV_dim, XV_maxdim, XV_chunkdim, XV_moffset,  XV_chunkBdim, offset0, true, fileIO->complex_tid);
        ArrXV.resize(RxLD, RvLD, RsLD);
     
      } 
     

     H5Gclose(visualGroup);
     

}   
    
Visualization_Data::~Visualization_Data() {
    delete FA_slphi;
    delete FA_slAp;
    delete FA_slBp;
    //    delete FA_slF1;
    //    delete FA_sln;
    //      delete FA_sldn;
    //    delete FA_slT;
    if(visXV == true) delete FA_XV;
    delete FA_slphiTime;

};

void Visualization_Data::writeData(Timing timing, double dt, bool force) 
{
        if (timing.check(dataOutputVisual,dt) || (timing.step == 1) || force) {
            FA_slphi->write(fields->phi);
            FA_slAp ->write(fields->Ap);
            FA_slBp ->write(fields->Bp);
            FA_slphiTime->write(&timing);

/* 
            // get slides for velocity plot
            for(int y_k = NkyLlD; y_k < NkyLlD + 1 ; y_k++) V3(RxLD, RvLD, y_k) = vlasov->f0(RxLD, y_k, NzLlD, RvLD, NmLlD, NsLlD);// + vlasov->f(RxLD, y_k, NzLlD, RvLD, NmLlD, NsLlD);
            FA_slF1 ->write(V3);


 * */
      if(visXV == true) {
            ArrXV(RxLD, RvLD, RsLD) = vlasov->f(RxLD, 0, NzLlD, RvLD, NmLlD, RsLD);
            if(plasma->global == false) ArrXV(RxLD, RvLD, RsLD) += vlasov->f0(RxLD, 0, NzLlD, RvLD, NmLlD, RsLD);
            FA_XV->write(ArrXV);

      }


      parallel->print("Wrote Visual  data ... "); 
      }
}
     
