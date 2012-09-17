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

int GC2, GC4, Nq;

Fields::Fields(Setup *setup, Grid *_grid, Parallel *_parallel, FileIO *fileIO, Geometry *_geo)  : 
grid(_grid), parallel(_parallel), geo(_geo), 
Q  (FortranArray<4>()),  Qm(FortranArray<4>()),  Field0(FortranArray<4>()), Field(FortranArray<6>()), solveEq(0),

SendXu(FortranArray<6>()), SendXl(FortranArray<6>()), SendZu(FortranArray<6>()), SendZl(FortranArray<6>()), 
RecvXu(FortranArray<6>()), RecvXl(FortranArray<6>()), RecvZu(FortranArray<6>()), RecvZl(FortranArray<6>()) 

{
   GC2 = 2; GC4 = 4, Nq=plasma->nfields;
 
   RFields.setRange(1, plasma->nfields);
   // Note : Fixed fields are initialized in Init.cpp
   solveEq |=  ((plasma->nfields >= 1) && (setup->get("Init.FixedPhi", ".0") == ".0")) ? Field::phi : 0;
   solveEq |=  ((plasma->nfields >= 2) && (setup->get("Init.FixedAp" , ".0") == ".0")) ? Field::Ap  : 0;
   solveEq |=  ((plasma->nfields >= 3) && (setup->get("Init.FixedBp" , ".0") == ".0")) ? Field::Bpp : 0;


   // for phi terms
   allocate(RxLB4,RkyLD,RzLB, RmLB, RsLB, RFields, Field);
   allocate(RxLD ,RkyLD,RzLD, RFields, Q, Qm, Field0);
      

   if(plasma->nfields >= Field::phi) { phi.reference(Field(RxLB4,RkyLD,RzLB,RmLB, RsLB, Field::phi)); phi = 0.; }
   if(plasma->nfields >= Field::Ap ) { Ap.reference (Field(RxLB4,RkyLD,RzLB,RmLB, RsLB, Field::Ap )); Ap  = 0.; }
   if(plasma->nfields >= Field::Bp ) { Bp.reference (Field(RxLB4,RkyLD,RzLB,RmLB, RsLB, Field::Bp )); Bp  = 0.; }
    
   //  brackets should be 1/2 but due to numerical errors, we should calculate it ourselves, see Dannert[2] 
   Yeb = (1./sqrt(M_PI) * sum(pow2(V(RvLD)) * exp(-pow2(V(RvLD)))) * dv) * geo->eps_hat * plasma->beta; 

   // Allocate boundary conditions, allocate Send/Recv buffers, note we have 4 ghost cells for X, 0 for Y
   allocate(RB4 , RkyLD, RzLD, RmLD, RsLD, RFields, SendXu, SendXl, RecvXu, RecvXl);
   //allocate(RxLD, RB4  , RzLD, RmLD, RsLD, RFields, SendYu, SendYl, RecvYu, RecvYl);
   allocate(RxLD, RkyLD, RB  , RmLD, RsLD, RFields, SendZu, SendZl, RecvZu, RecvZl);

   initDataOutput(setup, fileIO);
} 

Fields::~Fields() {

   closeData() ;

}


void Fields::solve(Array6C f0, Array6C  f, Timing timing)
{
  
  // calculate source terms  Q (Q is overwritten in the first iteration )
  for(int s = NsLlD, loop=0; s <= NsLuD; s++) { for(int m = NmLlD; m <= NmLuD; m++, loop++) {

      if(solveEq & Field::phi) calculateChargeDensity               ((A6zz) f0.dataZero(), (A6zz) f.dataZero(), (A4zz) Field0.dataZero(),               m, s);
      if(solveEq & Field::Ap ) calculateParallelCurrentDensity      ((A6zz) f0.dataZero(), (A6zz) f.dataZero(), (A4zz) Field0.dataZero(), V.dataZero(), m, s);
      if(solveEq & Field::Bpp) calculatePerpendicularCurrentDensity ((A6zz) f0.dataZero(), (A6zz) f.dataZero(), (A4zz) Field0.dataZero(), M.dataZero(), m, s);
  
      // thus uses AllReduce, Reduce is more effective (with false flag...)
      parallel->collect(Field0, OP_SUM, DIR_V, true); 
 
      // OPTIM : Normally we would decompose in m&s, thus no need for Qm                         
      // backward-transformation from gyro-center to drift-center 
      if (parallel->Coord[DIR_V] == 0) {

            gyroAverage(Field0, Qm, m, s, true);            
            // Lambda function to integrate over m ( for source terms), If loop=0, we overwrite value of Q
           if(loop==0) [=] (CComplex Q[Nq][NzLD][NkyLD][NxLD], CComplex Qm[Nq][NzLD][NkyLD][NxLD]) { Q[:][:][:][:]  = Qm[:][:][:][:]; } ((A4zz) Q.data(), (A4zz) Qm.data()); 
           else        [=] (CComplex Q[Nq][NzLD][NkyLD][NxLD], CComplex Qm[Nq][NzLD][NkyLD][NxLD]) { Q[:][:][:][:] += Qm[:][:][:][:]; } ((A4zz) Q.data(), (A4zz) Qm.data()); 
      
      }
      
   }  } // for m, for s

   /////////////////////////////// Solve for the corresponding fields //////////////
   // Note :  Fields are only solved by root nodes  (X=0, V=0, S=0), Gyro-averaging is done for (V=0)
   if(parallel->Coord[DIR_V] == 0) {

      // integrate over mu-space and over species
      parallel->collect(Q, OP_SUM, DIR_MS);

      // Solve field equation in drift coordinates
      // This routine is solved only on a part of notes when decomposed in m thus efficiency is crucial
      if(parallel->Coord[DIR_MS] == 0) solveFieldEquations((A4zz) Q.dataZero(), (A4zz) Field0.dataZero());

      parallel->send(Field0, DIR_MS);    

      // Gyro-averaging procedure for each species and magnetic moment ( drift-coord -> gyro-coord )
      // OPTIM : We can skip foward transform after first call
      for(int s = NsLlD; s <= NsLuD; s++) { for(int m = NmLlD; m <= NmLuD; m++) {

           // Field has complicated stride, thus cannot be easily used in FFTSolver.
           // Use temporary Qm, and copy to Fields afterwards.

           // forward-transformation from drift-center -> gyro-center 
           gyroAverage(Field0, Qm, m, s, true);

           // Field is first index, is last index not better ? e.g. write as vector ?
           [=] (CComplex Qm[Nq][NzLD][NkyLD][NxLD], CComplex Field[Nq][NsLD][NmLD][NzLB][NkyLD][NxLB+4]) 
           { 
                  Field[1:Nq][s][m][NzLlD:NzLD][NkyLlD:NkyLD][NxLlD:NxLD] = Qm[1:Nq][NzLlD:NzLD][NkyLlD:NkyLD][NxLlD:NxLD]   ;
           } ((A4zz) Qm.dataZero(), (A6zz) Field.dataZero());
        
      } }
   
      updateBoundary();
   
   }   
  
   parallel->send(Field, DIR_V);
        
   return;

}


void Fields::calculateChargeDensity(const CComplex f0         [NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB],
                                    const CComplex f          [NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB],
                                    CComplex Field0     [Nq]        [NzLD][NkyLD][NxLD]      ,
                                    const int m, const int s) 
{
 
   // In case of a full-f simulation the Maxwellian is subtracted

   // re-normalize with \f[ \hat{v}^3
   //const double pqnB_dvdm = M_PI * plasma->species[s].n0 * plasma->species[s].q   * plasma->B0 * dv * grid->dm(m) ;
   const double pqnB_dvdm = M_PI * plasma->species[s].q   * plasma->B0 * dv * grid->dm(m) ;
    
   omp_for_C2(int z=NzLlD; z<= NzLuD;z++) {  for(int y_k=NkyLlD; y_k<= NkyLuD; y_k++) { for(int x=NxLlD; x<= NxLuD;x++) {

              Field0[Q::rho][z][y_k][x] = ( __sec_reduce_add(f [s][m][z][y_k][x][NvLlD:NvLD]) 
                        - (plasma->global ? __sec_reduce_add(f0[s][m][z][y_k][x][NvLlD:NvLD]) : 0)) * pqnB_dvdm;
     
   } } } // z, y_k, x
   
   return; 
}



void Fields::calculateParallelCurrentDensity(const CComplex f0[NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB],
                                             const CComplex f [NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB],
                                             CComplex Field0[Nq][NzLD][NkyLD][NxLD],
                                             const double V[NvGB], const int m, const int s          ) 
{
  
   const double qa_dvdm = plasma->species[s].q * plasma->species[s].alpha  * plasma->B0 * M_PI * dv * grid->dm(m) ;
   
   omp_for_C2(int z=NzLlD; z<= NzLuD;z++) {  for(int y_k=NkyLlD; y_k<= NkyLuD; y_k++) { for(int x=NxLlD; x<= NxLuD;x++) {

                Field0[Q::jp][z][y_k][x] = -__sec_reduce_add(V[NvLlD:NvLD] * f[s][m][z][y_k][x][NvLlD:NvLD]) * qa_dvdm;

   } } } // z, y_k, x

   return ; 

}

  
// Note : Did not checked correctness of this function
void Fields::calculatePerpendicularCurrentDensity(const CComplex f0 [NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB],
                                                  const CComplex f  [NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB],
                                                  CComplex Field0[Nq][NzLD][NkyLD][NxLD],
                                                  const double M[NmGB], const int m, const int s          ) 
{
   
   const double qan_dvdm = - plasma->species[s].q * plasma->species[s].alpha   * plasma->B0 * M_PI * dv * grid->dm(m) ;
   
   for(int z=NzLlD; z<= NzLuD;z++) { omp_for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int x=NxLlD; x<= NxLuD;x++){

      Field0[Q::jo][z][y_k][x] =  M[m] * __sec_reduce_add(f[s][m][z][y_k][x][NvLlD:NvLD]) * qan_dvdm;
            
   } } } // z, y_k, x

   return;
}

   
// @todo this needs some cleanup work  and retesting for 3-dimensional case.
// only support parallelized version (assuming we run always parallized).
void Fields::updateBoundary()
{
   updateBoundary((A6zz ) Field.dataZero(), 
                  (A6zz ) SendXl.data()   , (A6zz) SendXu.data(), (A6zz) RecvXl.data(), (A6zz) RecvXu.data(),
                  (A6zz ) SendZl.data()   , (A6zz) SendZu.data(), (A6zz) RecvZl.data(), (A6zz) RecvZu.data());
   return;
}



void Fields::updateBoundary(
         CComplex Field [Nq][NsLD][NmLD][NzLB][NkyLD][NxLB+4], 
         CComplex SendXl[Nq][NsLD][NmLD][NzLD][NkyLD][GC4   ], CComplex SendXu[Nq][NsLD][NmLD][NzLD][NkyLD][GC4 ],
         CComplex RecvXl[Nq][NsLD][NmLD][NzLD][NkyLD][GC4   ], CComplex RecvXu[Nq][NsLD][NmLD][NzLD][NkyLD][GC4 ],
         CComplex SendZl[Nq][NsLD][NmLD][GC2 ][NkyLD][NxLD  ], CComplex SendZu[Nq][NsLD][NmLD][GC2 ][NkyLD][NxLD],
         CComplex RecvZl[Nq][NsLD][NmLD][GC2 ][NkyLD][NxLD  ], CComplex RecvZu[Nq][NsLD][NmLD][GC2 ][NkyLD][NxLD])
{
   
   // X-Boundary (we have extended BC  - 4 ghost cells)
   SendXl[:][:][:][:][:][:] = Field[1:Nq][NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][NkyLlD:NkyLD][NxLlD  :4];
   SendXu[:][:][:][:][:][:] = Field[1:Nq][NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][NkyLlD:NkyLD][NxLuD-3:4];

   // Not required in Y
       
   // Z-Boundary (N z-we need to connect the magnetic field lines)
   // (For 2-d space simulations (x,y) boundaries are not required)
   if(Nz > 1) omp_for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int x=NxLlD; x<= NxLuD;x++) { 

        const CComplex a = ((CComplex) 0. + 1.j) *  (2.*M_PI * (2.* M_PI/Ly) * y_k);
      
        // NzLlD == NzGlD -> Connect only physcial boundaries after mode made one loop 
        SendZl[:][:][:][:][y_k][x-3] = Field[1:Nq][NsLlD:NsLD][NmLlD:NmLD][NzLlD  :2][y_k][x] * cexp( ((NzLuD == NzGuD) ? a : 0.) * geo->nu(x));
        SendZu[:][:][:][:][y_k][x-3] = Field[1:Nq][NsLlD:NsLD][NmLlD:NmLD][NzLuD-1:2][y_k][x] * cexp(-((NzLlD == NzGlD) ? a : 0.) * geo->nu(x));

   } }
   
   // Parallized version is using SendRecv call exchange ghostcells
   parallel->updateNeighbours(Fields::SendXl, Fields::SendXu,  Fields::SendYl, Fields::SendYu, Fields::SendZl,  Fields::SendZu, 
                              Fields::RecvXl, Fields::RecvXu,  Fields::RecvYl, Fields::RecvYu, Fields::RecvZl,  Fields::RecvZu); 

   // X-Boundary
   Field[1:Nq][NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][NkyLlD:NkyLD][NxLlB-2:4] = RecvXl[:][:][:][:][:][:];
   Field[1:Nq][NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][NkyLlD:NkyLD][NxLuD+1:4] = RecvXu[:][:][:][:][:][:];
 
   // Z-Boundary
   if( Nz > 1 ) {

   Field[1:Nq][NsLlD:NsLD][NmLlD:NmLD][NzLlB  :2][NkyLlD:NkyLD][NxLlD:NxLD] = RecvZl[:][:][:][:][:][:];
   Field[1:Nq][NsLlD:NsLD][NmLlD:NmLD][NzLuD+1:2][NkyLlD:NkyLD][NxLlD:NxLD] = RecvZu[:][:][:][:][:][:];
   
   }

   return; 

};





//////////////////////// Data I/O //////////////////


void Fields::initDataOutput(Setup *setup, FileIO *fileIO) {
    
   // Set sizes : Note, we use fortran ordering for field variables 
   hsize_t field_dim[]       = { grid->NzGD, grid->NkyGD, grid->NxGD  ,             1};
   hsize_t field_chunkBdim[] = { NzLB      , grid->NkyGD, NxLB+4      ,             1};
   hsize_t field_chunkdim[]  = { NzLD      , NkyLD      , NxLD        ,             1};
   hsize_t field_maxdim[]    = { grid->NzGD, grid->NkyGD, grid->NxGD  , H5S_UNLIMITED};
   hsize_t field_moffset[]   = { 2, 0, 4, 0 };
   hsize_t field_offset[]    = {NzLlB-1, NkyLlB, NxLlB-1, 0 }; 
     
   bool phiWrite = (parallel->Coord[DIR_VMS] == 0);
     
   hid_t fieldsGroup = check(H5Gcreate(fileIO->getFileID(), "/Fields",H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), DMESG("Error creating group file for Phi : H5Gcreate"));
     
   FA_phi      = new FileAttr("Phi" , fieldsGroup, 4, field_dim , field_maxdim   , field_chunkdim   , field_moffset    ,  field_chunkBdim  , field_offset, phiWrite && plasma->nfields >= 1, fileIO->complex_tid);
   FA_Ap       = new FileAttr("Ap"  , fieldsGroup, 4, field_dim , field_maxdim   , field_chunkdim   , field_moffset    ,  field_chunkBdim  , field_offset, phiWrite && plasma->nfields >= 2, fileIO->complex_tid);
   FA_Bp       = new FileAttr("Bp"  , fieldsGroup, 4, field_dim , field_maxdim   , field_chunkdim   , field_moffset    ,  field_chunkBdim  , field_offset, phiWrite && plasma->nfields >= 3, fileIO->complex_tid);
   FA_phiTime  = fileIO->newTiming(fieldsGroup);
        
   H5Gclose(fieldsGroup);
      
   dataOutputFields        = Timing(setup->get("DataOutput.Fields.Step", -1)       , setup->get("DataOutput.Fields.Time", -1.));

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


