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

Fields::Fields(Setup *setup, Grid *_grid, Parallel *_parallel, FileIO *fileIO, Geometry<GKC_GEOMETRY> *_geo)  : 
grid(_grid),    parallel(_parallel), geo(_geo), 
Q(FortranArray<4>()), Field0(FortranArray<4>()), Field(FortranArray<6>()), solveEq(0)
{

   RFields.setRange(1, plasma->nfields);
   // Note : Fixed fields are initialized in Init.cpp
   solveEq |=  ((plasma->nfields >= 1) && (setup->get("Init.FixedPhi", ".0") == ".0")) ? Field::phi : 0;
   solveEq |=  ((plasma->nfields >= 2) && (setup->get("Init.FixedAp" , ".0") == ".0")) ? Field::Ap  : 0;
   solveEq |=  ((plasma->nfields >= 3) && (setup->get("Init.FixedBp" , ".0") == ".0")) ? Field::Bpp : 0;


   // for phi terms
   allocate(RxLB4,RkyLD,RzLB, RmLB, RsLB, RFields, Field);
   allocate(RxLD ,RkyLD,RzLD, RFields, Q, Field0);
      

   if(plasma->nfields >= Field::phi) { phi.reference(Field(RxLB4,RkyLD,RzLB,RmLB, RsLB, Field::phi)); phi = 0.; }
   if(plasma->nfields >= Field::Ap ) { Ap.reference (Field(RxLB4,RkyLD,RzLB,RmLB, RsLB, Field::Ap )); Ap  = 0.; }
   if(plasma->nfields >= Field::Bp ) { Bp.reference (Field(RxLB4,RkyLD,RzLB,RmLB, RsLB, Field::Bp )); Bp  = 0.; }
    
   //  brackets should be 1/2 but due to numerical errors, we should calculate it ourselves, see Dannert[2] 
   Yeb = (1./sqrt(M_PI) * sum(pow2(V(RvLD)) * exp(-pow2(V(RvLD)))) * dv) * geo->eps_hat * plasma->beta; 


   // Allocate boundary conditions, allocate Send/Recv buffers, note we have 4 ghost cells for X, 0 for Y
   allocate(RB4 , RkyLD, RzLD, RmLD, RsLD, RFields, SendXu, SendXl, RecvXu, RecvXl);
   allocate(RxLD, RB4  , RzLD, RmLD, RsLD, RFields, SendYu, SendYl, RecvYu, RecvYl);
   allocate(RxLD, RkyLD, RB  , RmLD, RsLD, RFields, SendZu, SendZl, RecvZu, RecvZl);

   initDataOutput(setup, fileIO);
} 


int Fields::solve(Array6z f0, Array6z  f, Timing timing, int rk_step)
{
   
  ////////////// / calculate source terms  (do we have to set it to zero ??) ///////////
   
   Q = 0.;
   for(int s = NsLlD; s <= NsLuD; s++) { for(int m = NmLlD; m <= NmLuD; m++) {

      if(solveEq & Field::phi) calculateChargeDensity               (f0, f, m, s);
      if(solveEq & Field::Ap ) calculateParallelCurrentDensity      (f0, f, m, s);
      if(solveEq & Field::Bpp) calculatePerpendicularCurrentDensity (f0, f, m, s);
     
      // thus uses AllReduce, Reduce is more effective (with false flag...)
      parallel->collect(Field0, OP_SUM, DIR_V, true); 
        
	  // gyro-average : back-transformation 
	  if (parallel->Coord(DIR_V) == 0) Q(RxLD, RkyLD, RzLD, RFields) += gyroAverage(Field0(RxLD, RkyLD, RzLD, RFields), m,s, Q::rho);
    
     } 
   //if(plasma->species(s).gyroModel != "Gyro") break; 

   } // for m, for s
    

   /////////////////////////////// Solve for the corresponding fields //////////////
  
   // Note :  Fields are only solved by root nodes in real space (X=0, V=0, S=0)
   //        Gyro-averaging is done for (V=0)
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
        
   return GKC_SUCCESS;

}


Array3z Fields::calculateChargeDensity(Array6z f0, Array6z f, const int m, const int s) 
{
  
   // In case of a full-f simulation the Maxwellian is subtracted

   // re-normalize with \f[ \hat{v}^3
   const double pqnB_dV = M_PI * plasma->species(s).n0 * plasma->species(s).q   * plasma->B0 * dv * grid->dm(m) ;
    
   for(int z=NzLlD; z<= NzLuD;z++) {  omp_for(int y_k=NkyLlD; y_k<= NkyLuD; y_k++) { for(int x=NxLlD; x<= NxLuD;x++) {

              Field0(x, y_k, z, Q::rho) =  ( sum(f(x, y_k, z, RvLD, m, s) - (plasma->global ? sum(f0(x,y_k,z,RvLD,m,s)) : 0. ) ) ) * pqnB_dV;
     
   } } }
    
   return Q(RxLD, RkyLD, RzLD, Q::rho); 
}

Array3z Fields::calculateParallelCurrentDensity(Array6z f0, Array6z f, const int m, const int s) 
{
  
   const double qa_dV = plasma->species(s).alpha * plasma->species(s).q   * plasma->B0 * M_PI * dv * grid->dm(m) ;
   
   for(int z=NzLlD; z<= NzLuD;z++) {  omp_for(int y_k=NkyLlD; y_k<= NkyLuD; y_k++) { for(int x=NxLlD; x<= NxLuD;x++) {

                Field0(x, y_k, z, Q::jp) = - sum(V(RvLD) * f(x, y_k, z, RvLD, m, s)) * qa_dV;

   } } }

   return Q(RxLD, RkyLD, RzLD, Q::jp); 

}

  
// Note : Did not checked correctness of this function
Array3z Fields::calculatePerpendicularCurrentDensity(Array6z f0, Array6z f, const int m, const int s) 
{
   
   const double qan_dV = - plasma->species(s).q * plasma->species(s).alpha   * plasma->B0 * M_PI * dv * grid->dm(m) ;
   
   for(int z=NzLlD; z<= NzLuD;z++) { omp_for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int x=NxLlD; x<= NxLuD;x++){

      Field0(x, y_k, z, Q::jo) = M(m) * sum(f(x, y_k, z, RvLD, m, s)) * qan_dV;
            
   } } }

   return Q(RxLD, RkyLD, RzLD, Q::jo); 
}

   
// @todo this needs some cleanup work  and retesting for 3-dimensional case.
// only support parallelized version.
int Fields::setBoundary(Array6z  A) {

#ifdef GKC_PARALLEL_MPI

   SendXl(RB4  , RkyLD, RzLD, RmLD, RsLD, RFields) = A(Range(NxLlD  , NxLlD+3), RkyLD, RzLD, RmLD, RsLD, RFields);
   SendXu(RB4  , RkyLD, RzLD, RmLD, RsLD, RFields) = A(Range(NxLuD-3, NxLuD  ), RkyLD, RzLD, RmLD, RsLD, RFields);

   // Not required
   // SendYl(RxLD, RB4  , RzLD, RmLD, RsLD, RFields) = A(RxLD, Range(NyLlD  , NyLlD+3), RzLD, RmLD, RsLD, RFields);
   // SendYu(RxLD, RB4  , RzLD, RmLD, RsLD, RFields) = A(RxLD, Range(NyLuD-3, NyLuD  ), RzLD, RmLD, RsLD, RFields);
        
   // For z-we need to connect the magnetic field lines
   for(int x=NxLlD; x<= NxLuD;x++) { for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) {

      const ShearB b = geo->getYPos(x,y_k);
      double res;
      // original
      // SendZl(x, y, RB, RmLD, RsLD, RFields) = A(x, b.ly, Range(NzLlD  , NzLlD+1), RmLD, RsLD, RFields);
      // SendZu(x, y, RB, RmLD, RsLD, RFields) = A(x, b.uy, Range(NzLuD-1, NzLuD  ), RmLD, RsLD, RFields);
            
      // test (but what kind of ? not documented, not remembering :(  )
      SendZl(x, b.ly, RB, RmLD, RsLD, RFields) = A(x, y_k , Range(NzLlD  , NzLlD+1), RmLD, RsLD, RFields);
      SendZu(x, y_k , RB, RmLD, RsLD, RFields) = A(x, b.ly, Range(NzLuD-1, NzLuD  ), RmLD, RsLD, RFields);

   }  }

   // Parallized version is using SendRecv call
   parallel->updateNeighbours(  SendXl,   SendXu,  SendYl, SendYu, SendZl,  SendZu, 
                                RecvXl,   RecvXu,  RecvYl, RecvYu, RecvZl,  RecvZu); 
        
   A(Range(NxLlB-2, NxLlD-1), RkyLD, RzLD, RmLD, RsLD, RFields) = RecvXl(RB4, RkyLD, RzLD, RmLD, RsLD, RFields);
   A(Range(NxLuD+1, NxLuB+2), RkyLD, RzLD, RmLD, RsLD, RFields) = RecvXu(RB4, RkyLD, RzLD, RmLD, RsLD, RFields);
 
   // A(RxLD, Range(NyLlB-2, NyLlD-1), RzLD, RmLD, RsLD, RFields) = RecvYl(RxLD, RB4, RzLD, RmLD, RsLD, RFields);
   // A(RxLD, Range(NyLuD+1, NyLuB+2), RzLD, RmLD, RsLD, RFields) = RecvYu(RxLD, RB4, RzLD, RmLD, RsLD, RFields);
        
   A(RxLD, RkyLD, Range(NzLlB  , NzLlB+1), RmLD, RsLD, RFields) = RecvZl(RxLD, RkyLD, RB, RmLD, RsLD, RFields);
   A(RxLD, RkyLD, Range(NzLuD+1, NzLuB  ), RmLD, RsLD, RFields) = RecvZu(RxLD, RkyLD, RB, RmLD, RsLD, RFields);

#else

   A(Range(NxLlB-2, NxLlB+1), RkyLD, RzLD, RmLD, RsLD, RFields) = A(Range(NxLuD-3, NxLuD),   RkyLD, RzLD, RmLD, RsLD, RFields);
   A(Range(NxLuD+1, NxLuB+2), RkyLD, RzLD, RmLD, RsLD, RFields) = A(Range(NxLlD, NxLlD+3),   RkyLD, RzLD, RmLD, RsLD, RFields);
    
   // A(RxLD, Range(NyLlB-2, NyLlB+1), RzLD, RmLD, RsLD, RFields) = A(RxLD, Range(NyLuD-3, NyLuD), RzLD, RmLD, RsLD, RFields);
   // A(RxLD, Range(NyLuD+1, NyLuB+2), RzLD, RmLD, RsLD, RFields) = A(RxLD, Range(NyLlD, NyLlD+3), RzLD, RmLD, RsLD, RFields);

   // For z-we need to connect the magnetic field lines
   for(int x=NxLlD; x<= NxLuD;x++) { omp_for(int y=NyLlD; y<= NyLuD;y++) {

      ShearB b = geo->getYPos(x,y);
      A(x, y, Range(NzLlB, NzLlB+1), RmLD, RsLD, RFields)   = A(x, b.ly, Range(NzLuD-1, NzLuD), RmLD, RsLD, RFields);
      A(x, y, Range(NzLuD+1, NzLuB), RmLD, RsLD, RFields)   = A(x, b.uy, Range(NzLlD, NzLlD+1), RmLD, RsLD, RFields);

   }   } 

   // Interpolation of connecting field lines if they do not math exactly
   // A(x, y, Range(NzLlB, NzLlB+1), RmLD, RsLD, RFields)  = b.y0 * A(x, b.ypos, Range(NzLuD-1, NzLuD), RmLD, RsLD, RFields) + (1. - b.y0) *  A(x, b.ypos-1, Range(NzLuD-1, NzLuD), RmLD, RsLD, RFields);
   // A(x, y, Range( NzLuD+1, NzLuB), RmLD, RsLD, RFields) = b.y0 * A(x, y, Range(NzLlD, NzLlD+1), RmLD, RsLD, RFields)      + (1. - b.y0) *  A(x, y+1     , Range(NzLlD, NzLlD+1), RmLD, RsLD, RFields);
   // A(x, y,    Range(NzLlB, NzLlB+1), RmLD, RsLD, RFields)    = b.y0 * A(x, b.ypos, Range(NzLuD-1, NzLuD), RmLD, RsLD, RFields) + (1. - b.y0) *  A(x, b.ypos-1, Range(NzLuD-1, NzLuD), RmLD, RsLD, RFields);
   // A(x, b.ypos, Range( NzLuD+1, NzLuB), RmLD, RsLD, RFields) = b.y0 * A(x, y, Range(NzLlD, NzLlD+1), RmLD, RsLD, RFields)      + (1. - b.y0) *  A(x, y+1     , Range(NzLlD, NzLlD+1), RmLD, RsLD, RFields);
   // A(x, y,    Range(NzLlB, NzLlB+1), RmLD, RsLD, RFields)    = A(x, b.ypos, Range(NzLuD-1, NzLuD), RmLD, RsLD, RFields) ;
   // A(x, b.ypos, Range( NzLuD+1, NzLuB), RmLD, RsLD, RFields) = A(x, y, Range(NzLlD, NzLlD+1), RmLD, RsLD, RFields);
     
#endif  // GKC_PARALLEL_MPI

   return GKC_SUCCESS; 

};


void Fields::initDataOutput(Setup *setup, FileIO *fileIO) {
    
   // Set sizes : Note, we use fortran ordering for field variables 
   hsize_t field_dim[]       = { grid->NzGD, grid->NkyGD, grid->NxGD  ,             1};
   hsize_t field_chunkBdim[] = { NzLB      , grid->NkyGD, NxLB+4      ,             1};
   hsize_t field_chunkdim[]  = { NzLD      , NkyLD      , NxLD        ,             1};
   hsize_t field_maxdim[]    = { grid->NzGD, grid->NkyGD, grid->NxGD  , H5S_UNLIMITED};
   hsize_t field_moffset[]   = { 2, 0, 4, 0 };
   hsize_t field_offset[]    = {NzLlB-1, NkyLlB, NxLlB-1, 0 }; 
     
   bool phiWrite = (parallel->Coord(DIR_VMS) == 0);
     
   hid_t fieldsGroup = check(H5Gcreate(fileIO->getFileID(), "/Potential",H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), DMESG("Error creating group file for Phi : H5Gcreate"));
     
   FA_phi      = new FileAttr("Phi" , fieldsGroup, 4, field_dim , field_maxdim   , field_chunkdim   , field_moffset    ,  field_chunkBdim  , field_offset, phiWrite && plasma->nfields >= 1, fileIO->complex_tid);
   FA_Ap       = new FileAttr("Ap"  , fieldsGroup, 4, field_dim , field_maxdim   , field_chunkdim   , field_moffset    ,  field_chunkBdim  , field_offset, phiWrite && plasma->nfields >= 2, fileIO->complex_tid);
   FA_Bp       = new FileAttr("Bp"  , fieldsGroup, 4, field_dim , field_maxdim   , field_chunkdim   , field_moffset    ,  field_chunkBdim  , field_offset, phiWrite && plasma->nfields >= 3, fileIO->complex_tid);
   FA_phiTime  = fileIO->newTiming(fieldsGroup);
        
   H5Gclose(fieldsGroup);
      
     
   dataOutputFields        = Timing(setup->get("DataOutput.Phi.Step", -1)       , setup->get("DataOutput.Phi.Time", -1.));

}   



// to do : Improve writeMessage (timestep, offset, etc.)
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


// @todo : how to interact with FieldFFT, FieldHermite ?
void Fields::printOn(ostream &output) const {

         output   << "Poisson    |  " << "Base class" << std::endl;
         output   << "Ampere     |  " << ((plasma->nfields >= 2) ? "beta :  " + Num2String(plasma->beta) : " --- no electromagnetic effects ---") << std::endl;
         output   << "B_parallel |  " << ((plasma->nfields >= 3) ? "beta :  " + Num2String(plasma->beta) : " --- no electromagnetic effects ---") << std::endl;
         output   << "           |  Pert. : " << ((!(solveEq & Field::Ap) && (plasma->nfields >= 2))     ? ApPerturbationStr  : " ") << std::endl;
}


Fields::~Fields() {

   closeData() ;

}

