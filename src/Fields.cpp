/*
 * =====================================================================================
 *
 *       Filename: Fields.cpp
 *
 *    Description: Main Interface for solving fields equations. Implements
 *                 source terms calculations.
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include "Fields.h"


Range RFields;

Fields::Fields(Setup *setup, Grid *_grid, Parallel *_parallel, FileIO *fileIO, Geometry<HELIOS_GEOMETRY> *_geo)  : 
  grid(_grid),    parallel(_parallel), geo(_geo), 
     Q(FortranArray<4>()),     Field0(FortranArray<4>()),      Field(FortranArray<6>()), solveEq(0)
{

    RFields.setRange(1, plasma->nfields);
    // Note : Fixed fields are initialized in Init.cpp
    solveEq = 0;
    solveEq |=  ((plasma->nfields >= 1) && (setup->get("Init.FixedPhi", ".0") == ".0")) ? Field::phi : 0;
    solveEq |=  ((plasma->nfields >= 2) && (setup->get("Init.FixedAp" , ".0") == ".0")) ? Field::Ap  : 0;
    solveEq |=  ((plasma->nfields >= 3) && (setup->get("Init.FixedBp" , ".0") == ".0")) ? Field::Bpp : 0;



    // for phi terms
    Field.resize(RxLB4,RkyLD,RzLB,RmLB, RsLB, RFields); Field  = 0.; 
    allocate(RxLD,RkyLD,RzLD, RFields, Q, Field0);
      

    if(plasma->nfields >= Field::phi) { phi.reference(Field(RxLB4,RkyLD,RzLB,RmLB, RsLB, Field::phi)); phi = 0.; }
    if(plasma->nfields >= Field::Ap ) { Ap.reference (Field(RxLB4,RkyLD,RzLB,RmLB, RsLB, Field::Ap )); Ap  = 0.; }
    if(plasma->nfields >= Field::Bp ) { Bp.reference (Field(RxLB4,RkyLD,RzLB,RmLB, RsLB, Field::Bp )); Bp  = 0.; }
    
    //  brackets should be 1/2 but due to numerical error we should include calucalte ourselves, see Dannert[2] 
    Yeb = (1./sqrt(M_PI) * sum(pow2(V(RvLD)) * exp(-pow2(V(RvLD)))) * dv) * geo->eps_hat * plasma->beta; 


    // Allocate boundary conditions, allocate Send/Recv buffers, note we have 4 ghost cells for X & Y
    allocate(RB4 , RkyLD, RzLD, RmLD, RsLD, RFields, SendXu, SendXl, RecvXu, RecvXl);
    allocate(RxLD, RB4  , RzLD, RmLD, RsLD, RFields, SendYu, SendYl, RecvYu, RecvYl);
    allocate(RxLD, RkyLD, RB  , RmLD, RsLD, RFields, SendZu, SendZl, RecvZu, RecvZl);


    gyroAverageModel = setup->get("Fields.GyroAvrgModel", (Nm == 1) ? "Gyro-1" : "Gyro");


    initDataOutput(setup, fileIO);
} 


   // * 1. Integrate over v and back-transform (pull-back transform) over m to get the charge number density.
   // only DIR_V = 0 have density unequal zero
   // calculare sources terms
int Fields::solve(Array6z f0, Array6z  f, Timing timing, int rk_step)
{
   
  // calculate source terms
   Q = 0.;
   for(int s = NsLlD; s <= NsLuD; s++) { for(int m = NmLlD; m <= NmLuD; m++) {

        if(solveEq & Field::phi) calculateChargeDensity               (f0, f, m, s);
        if(solveEq & Field::Ap ) calculateParallelCurrentDensity      (f0, f, m, s);
        if(solveEq & Field::Bpp) calculatePerpendicularCurrentDensity (f0, f, m, s);
     
        parallel->collect(Field0, OP_SUM, DIR_V, true); //parallel->collect(rho, OP_SUM, DIR_V, false);
        // perform gyroAverage (note, gyorAvearge gives input array if in drift kinetic mode)
        // first test, skip electrons
        if (parallel->Coord(DIR_V) == 0) Q(RxLD, RkyLD, RzLD, RFields) += gyroAverage(Field0(RxLD, RkyLD, RzLD, RFields), m,s, Q::rho);
    
   } } // for m, for s
    


   if(parallel->Coord(DIR_V) == 0) {

     // integrate over mu-space and over species
     parallel->collect(Q, OP_SUM, DIR_MS);
     if((parallel->Coord(DIR_MS) == 0)) Field0(RxLD, RkyLD, RzLD, RFields) = solveFieldEquations(Q, timing);
     parallel->send(Field0, DIR_MS);    

	    
     // gyro-averaging procedure for each species and magnetic moment
     for(int s = NsLlD; s <= NsLuD; s++) { for(int m = NmLlD; m <= NmLuD; m++) {

            Field(RxLD,RkyLD,RzLD,m, s, RFields ) = gyroAverage(Field0(RxLD, RkyLD, RzLD, RFields ), m, s, Field::phi, true);
        
	 } }
   
     setBoundary(Field);
  }   
    
  parallel->send(Field, DIR_V);
        
   
  return HELIOS_SUCCESS;

}


// integrate Phasespace function to get number density / charge density!
Array3z Fields::calculateChargeDensity(Array6z f0, Array6z f, const int m, const int s) {

   
     // re-normalize with \f[ \hat{v}^3
     const double pqnB_d6Z = M_PI * plasma->species(s).n0 * plasma->species(s).q   * plasma->B0 * dv * grid->Ipol_M->d(m) ;
    
     // integrate over v
     for(int z=NzLlD; z<= NzLuD;z++) { 
       #pragma omp parallel for
       for(int y_k=NkyLlD; y_k<= NkyLuD; y_k++) { for(int x=NxLlD; x<= NxLuD;x++) {

              Field0(x, y_k, z, Q::rho) =  ( sum(f(x, y_k, z, RvLD, m, s) - (plasma->global ? sum(f0(x,y_k,z,RvLD,m,s)) : 0. ) ) ) * pqnB_d6Z;
     
     } } }
    

   return Q(RxLD, RkyLD, RzLD, Q::rho); 
}

// integrate Phasespace function to get number density / charge density!
Array3z Fields::calculateParallelCurrentDensity(Array6z f0, Array6z f, const int m, const int s) {
  
     // re-normalize with \f[ \hat{v}^3
     const double qa_d6Z = plasma->species(s).alpha * plasma->species(s).q   * plasma->B0 * M_PI * dv * grid->Ipol_M->d(m) ;
   
     // int_v  v_\parallel * g_j dv
     for(int z=NzLlD; z<= NzLuD;z++) { for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int x=NxLlD; x<= NxLuD;x++){

                Field0(x, y_k, z, Q::jp) = - sum(V(RvLD) * f(x, y_k, z, RvLD, m, s)) * qa_d6Z;

     } } }

    return Q(RxLD, RkyLD, RzLD, Q::jp); 
}

   
Array3z Fields::calculatePerpendicularCurrentDensity(Array6z f0, Array6z f, const int m, const int s) {
   
     // re-normalize with \f[ \hat{v}^3
     const double qn_d6Z = -plasma->species(s).alpha * plasma->species(s).q   * plasma->B0 * M_PI * dv * grid->Ipol_M->d(m) ;
   
     // int_v  v_\parallel * g_j dv
     for(int z=NzLlD; z<= NzLuD;z++) { for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int x=NxLlD; x<= NxLuD;x++){

                Field0(x, y_k, z, Q::jo) = M(m) * sum(f(x, y_k, z, RvLD, m, s)) * qn_d6Z;
            
     } } }

     return Q(RxLD, RkyLD, RzLD, Q::jo); 
}

   
   
   
int Fields::setBoundary(Array6z  A) {

#ifdef GKC_PARALLEL_MPI

        SendXl(RB4  , RkyLD, RzLD, RmLD, RsLD, Range::all()) = A(Range(NxLlD  , NxLlD+3), RkyLD, RzLD, RmLD, RsLD, Range::all());
        SendXu(RB4  , RkyLD, RzLD, RmLD, RsLD, Range::all()) = A(Range(NxLuD-3, NxLuD  ), RkyLD, RzLD, RmLD, RsLD, Range::all());

//        SendYl(RxLD, RB4  , RzLD, RmLD, RsLD, Range::all()) = A(RxLD, Range(NyLlD  , NyLlD+3), RzLD, RmLD, RsLD, Range::all());
//        SendYu(RxLD, RB4  , RzLD, RmLD, RsLD, Range::all()) = A(RxLD, Range(NyLuD-3, NyLuD  ), RzLD, RmLD, RsLD, Range::all());
        
        // For z-we need to connect the magnetic field lines
        for(int x=NxLlD; x<= NxLuD;x++) { for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) {
            const ShearB b = geo->getYPos(x,y_k);

            double res;
            // original
//            SendZl(x, y, RB, RmLD, RsLD, Range::all()) = A(x, b.ly, Range(NzLlD  , NzLlD+1), RmLD, RsLD, Range::all());
//            SendZu(x, y, RB, RmLD, RsLD, Range::all()) = A(x, b.uy, Range(NzLuD-1, NzLuD  ), RmLD, RsLD, Range::all());
            
            // test
            SendZl(x, b.ly, RB, RmLD, RsLD, Range::all()) = A(x, y_k, Range(NzLlD  , NzLlD+1), RmLD, RsLD, Range::all());
            SendZu(x, y_k, RB, RmLD, RsLD, Range::all())    = A(x, b.ly, Range(NzLuD-1, NzLuD  ), RmLD, RsLD, Range::all());

        }}

        parallel->updateNeighbours(  SendXl,   SendXu,   SendYl,   SendYu,  SendZl,  SendZu, 
                                     RecvXl,   RecvXu,   RecvYl,   RecvYu,  RecvZl,  RecvZu); 
        
        A(Range(NxLlB-2, NxLlD-1), RkyLD, RzLD, RmLD, RsLD, Range::all()) = RecvXl(RB4, RkyLD, RzLD, RmLD, RsLD, Range::all());
        A(Range(NxLuD+1, NxLuB+2), RkyLD, RzLD, RmLD, RsLD, Range::all()) = RecvXu(RB4, RkyLD, RzLD, RmLD, RsLD, Range::all());
 
//        A(RxLD, Range(NyLlB-2, NyLlD-1), RzLD, RmLD, RsLD, Range::all()) = RecvYl(RxLD, RB4, RzLD, RmLD, RsLD, Range::all());
//        A(RxLD, Range(NyLuD+1, NyLuB+2), RzLD, RmLD, RsLD, Range::all()) = RecvYu(RxLD, RB4, RzLD, RmLD, RsLD, Range::all());
        
        A(RxLD, RkyLD, Range(NzLlB  , NzLlB+1), RmLD, RsLD, Range::all()) = RecvZl(RxLD, RkyLD, RB, RmLD, RsLD, Range::all());
        A(RxLD, RkyLD, Range(NzLuD+1, NzLuB  ), RmLD, RsLD, Range::all()) = RecvZu(RxLD, RkyLD, RB, RmLD, RsLD, Range::all());

#else

        A(Range(NxLlB-2, NxLlB+1), RkyLD, RzLD, RmLD, RsLD, Range::all()) = A(Range(NxLuD-3, NxLuD),   RkyLD, RzLD, RmLD, RsLD, Range::all());
        A(Range(NxLuD+1, NxLuB+2), RkyLD, RzLD, RmLD, RsLD, Range::all()) = A(Range(NxLlD, NxLlD+3),   RkyLD, RzLD, RmLD, RsLD, Range::all());
    
//        A(RxLD, Range(NyLlB-2, NyLlB+1), RzLD, RmLD, RsLD, Range::all()) = A(RxLD, Range(NyLuD-3, NyLuD), RzLD, RmLD, RsLD, Range::all());
//        A(RxLD, Range(NyLuD+1, NyLuB+2), RzLD, RmLD, RsLD, Range::all()) = A(RxLD, Range(NyLlD, NyLlD+3), RzLD, RmLD, RsLD, Range::all());

        // For z-we need to connect the magnetic field lines
        for(int x=NxLlD; x<= NxLuD;x++) { for(int y=NyLlD; y<= NyLuD;y++) {
           ShearB b = geo->getYPos(x,y);
           A(x, y, Range(NzLlB, NzLlB+1), RmLD, RsLD, Range::all())   = A(x, b.ly, Range(NzLuD-1, NzLuD), RmLD, RsLD, Range::all());
           A(x, y, Range(NzLuD+1, NzLuB), RmLD, RsLD, Range::all())   = A(x, b.uy, Range(NzLlD, NzLlD+1), RmLD, RsLD, Range::all());
        }} 




     //   A(x, y, Range(NzLlB, NzLlB+1), RmLD, RsLD, Range::all())  = b.y0 * A(x, b.ypos, Range(NzLuD-1, NzLuD), RmLD, RsLD, Range::all()) + (1. - b.y0) *  A(x, b.ypos-1, Range(NzLuD-1, NzLuD), RmLD, RsLD, Range::all());
     //   A(x, y, Range( NzLuD+1, NzLuB), RmLD, RsLD, Range::all()) = b.y0 * A(x, y, Range(NzLlD, NzLlD+1), RmLD, RsLD, Range::all())      + (1. - b.y0) *  A(x, y+1     , Range(NzLlD, NzLlD+1), RmLD, RsLD, Range::all());
       // A(x, y,    Range(NzLlB, NzLlB+1), RmLD, RsLD, Range::all())    = b.y0 * A(x, b.ypos, Range(NzLuD-1, NzLuD), RmLD, RsLD, Range::all()) + (1. - b.y0) *  A(x, b.ypos-1, Range(NzLuD-1, NzLuD), RmLD, RsLD, Range::all());
       // A(x, b.ypos, Range( NzLuD+1, NzLuB), RmLD, RsLD, Range::all()) = b.y0 * A(x, y, Range(NzLlD, NzLlD+1), RmLD, RsLD, Range::all())      + (1. - b.y0) *  A(x, y+1     , Range(NzLlD, NzLlD+1), RmLD, RsLD, Range::all());
//        A(x, y,    Range(NzLlB, NzLlB+1), RmLD, RsLD, Range::all())    = A(x, b.ypos, Range(NzLuD-1, NzLuD), RmLD, RsLD, Range::all()) ;
//        A(x, b.ypos, Range( NzLuD+1, NzLuB), RmLD, RsLD, Range::all()) = A(x, y, Range(NzLlD, NzLlD+1), RmLD, RsLD, Range::all());
#endif


   return HELIOS_SUCCESS; 



  };


Fields::~Fields() {
      // ------------------ Phi Potential ----------------------//
    closeData() ;

}

  void Fields::initDataOutput(Setup *setup, FileIO *fileIO) {
     
     hsize_t field_dim[]      =  { grid->NzGD, grid->NkyGD, grid->NxGD  ,             1};
     hsize_t field_chunkBdim[] = { NzLB      , grid->NkyGD, NxLB+4      ,             1};
     hsize_t field_chunkdim[] =  { NzLD      , NkyLD      , NxLD        ,             1};
     hsize_t field_maxdim[]   =  { grid->NzGD, grid->NkyGD, grid->NxGD  , H5S_UNLIMITED};
     hsize_t field_moffset[]   = { 2, 0, 4, 0 };
     hsize_t field_offset[] = {NzLlB-1, NkyLlB, NxLlB-1, 0 }; 
     
     bool phiWrite = (parallel->Coord(DIR_VMS) == 0);
     
     hid_t fieldsGroup = check(H5Gcreate(fileIO->getFileID(), "/Potential",H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), DMESG("Error creating group file for Phi : H5Gcreate"));
     
     FA_phi      = new FileAttr("Phi" , fieldsGroup, 4, field_dim , field_maxdim   , field_chunkdim   , field_moffset    ,  field_chunkBdim  , field_offset, phiWrite && plasma->nfields >= 1);
     FA_Ap       = new FileAttr("Ap"  , fieldsGroup, 4, field_dim , field_maxdim   , field_chunkdim   , field_moffset    ,  field_chunkBdim  , field_offset, phiWrite && plasma->nfields >= 2);
     FA_Bp       = new FileAttr("Bp"  , fieldsGroup, 4, field_dim , field_maxdim   , field_chunkdim   , field_moffset    ,  field_chunkBdim  , field_offset, phiWrite && plasma->nfields >= 3);
     FA_phiTime  = fileIO->newTiming(fieldsGroup);
        
     check(H5LTset_attribute_string(fieldsGroup, ".", "GyroAverageModel", gyroAverageModel.c_str()), DMESG("H5LTset_attribute"));
        
     H5Gclose(fieldsGroup);
      
     
     dataOutputFields        = Timing(setup->get("DataOutput.Phi.Step", -1)       , setup->get("DataOutput.Phi.Time", -1.));

  }   
     
  void Fields::writeData(Timing timing, double dt) 
{
      if (timing.check(dataOutputFields, dt)       )   {
         FA_phi->write(phi);
         FA_Ap->write(Ap);
         FA_Bp->write(Bp);
         FA_phiTime->write(&timing);
      
        writeMessage("Wrote Potential data ... "); 
      }
  } 
      
  void Fields::closeData() {
      delete FA_phi;
      delete FA_Ap;
      delete FA_Bp;
      delete FA_phiTime;
  }
