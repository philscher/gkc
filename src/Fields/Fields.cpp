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




int GC2, GC4, Nq;


Fields::Fields(Setup *setup, Grid *_grid, Parallel *_parallel, FileIO *fileIO, Geometry *_geo)  : 

grid(_grid), parallel(_parallel), geo(_geo), solveEq(0)

{
   GC2 = 2; GC4 = 4, Nq=plasma->nfields;
 
   // Note : Fixed fields are initialized in Init.cpp
   solveEq |=  ((plasma->nfields >= 1) && (setup->get("Init.FixedPhi", ".0") == ".0")) ? Field::phi : 0;
   solveEq |=  ((plasma->nfields >= 2) && (setup->get("Init.FixedAp" , ".0") == ".0")) ? Field::Ap  : 0;
   solveEq |=  ((plasma->nfields >= 3) && (setup->get("Init.FixedBp" , ".0") == ".0")) ? Field::Bpp : 0;


   // for phi terms
   ArrayField = nct::allocate(nct::Range(1,Nq), grid->RsLD, grid->RmLD, grid->RzLB, grid->RkyLD, nct::Range(NxLlD-4, NxLD+8));
   ArrayField(&Field);

   ArrayField0 = nct::allocate(nct::Range(1,Nq), grid->RzLD, grid->RkyLD, grid->RxLD);
   ArrayField0(&Q, &Qm, &Field0);

   // Allocate boundary conditions, allocate Send/Recv buffers, note we have 4 ghost cells for X, 0 for Y
   ArrayBoundX = nct::allocate(nct::Range(0, 4 * NkyLD * NzLD * NmLD * NsLD * Nq));
   ArrayBoundX(&SendXl, &SendXu, &RecvXl, &RecvXu);
   
   ArrayBoundZ = nct::allocate(nct::Range(0, NxLD * NkyLD * 2 * NmLD * NsLD * Nq));
   ArrayBoundZ(&SendZl, &SendZu, &RecvZl, &RecvZu);
       
      
   //  brackets should be 1/2 but due to numerical errors, we should calculate it ourselves, see Dannert[2] 
   Yeb = (1./sqrt(M_PI) * __sec_reduce_add(pow2(V[NvLlD:NvLD]) * exp(-pow2(V[NvLlD:NvLD]))) * dv) * geo->eps_hat * plasma->beta; 

   initDataOutput(setup, fileIO);
} 

Fields::~Fields() 
{
   closeData() ;
}


void Fields::solve(CComplex *f0, CComplex *f, Timing timing)
{
  
  // calculate source terms  Q (Q is overwritten in the first iteration )
  for(int s = NsLlD, loop=0; s <= NsLuD; s++) { for(int m = NmLlD; m <= NmLuD; m++, loop++) {

      if(solveEq & Field::phi) calculateChargeDensity               ((A6zz) f0, (A6zz) f, (A4zz) Field0,               m, s);
      if(solveEq & Field::Ap ) calculateParallelCurrentDensity      ((A6zz) f0, (A6zz) f, (A4zz) Field0, V, m, s);
      if(solveEq & Field::Bpp) calculatePerpendicularCurrentDensity ((A6zz) f0, (A6zz) f, (A4zz) Field0, M, m, s);
  
      // thus uses AllReduce, Reduce is more effective (with false flag...)
      parallel->reduce(ArrayField0.data(Field0), Op::SUM, DIR_V, ArrayField0.getNum(), true); 
      
      // OPTIM : Normally we would decompose in m&s, thus no need for Qm                         
      // backward-transformation from gyro-center to drift-center 
      if (parallel->Coord[DIR_V] == 0) {

            gyroAverage((A4zz) Field0, (A4zz) Qm, m, s, false);           

            // Lambda function to integrate over m ( for source terms), If loop=0, we overwrite value of Q
            [=] (CComplex Q[Nq][NzLD][NkyLD][NxLD], CComplex Qm[Nq][NzLD][NkyLD][NxLD])
            {
               if  (loop == 0) Q[1:Nq][NzLlD:NzLD][NkyLlD:NkyLD][NxLlD:NxLD]   = Qm[1:Nq][NzLlD:NzLD][NkyLlD:NkyLD][NxLlD:NxLD];
               else            Q[1:Nq][NzLlD:NzLD][NkyLlD:NkyLD][NxLlD:NxLD]  += Qm[1:Nq][NzLlD:NzLD][NkyLlD:NkyLD][NxLlD:NxLD];

            } ((A4zz) Q, (A4zz) Qm);
      }
      
   }  } // for m, s

   /////////////////////////////// Solve for the corresponding fields ////////////////////////////////
   // Note :  Fields are only solved by root nodes  (X=0, V=0, S=0), Gyro-averaging is done for (V=0)
   if(parallel->Coord[DIR_V] == 0) {

      // integrate over mu-space and over species
      parallel->reduce(ArrayField0.data(Q), Op::SUM, DIR_MS, ArrayField0.getNum(), true); 

      // Solve field equation in drift coordinates
      // This routine is solved only on rood nodes, thus efficiency is crucial for scalability
      if(parallel->Coord[DIR_MS] == 0) solveFieldEquations((A4zz) Q, (A4zz) Field0);

      parallel->bcast(ArrayField0.data(Field0), DIR_MS, ArrayField0.getNum()); 

      // Gyro-averaging procedure for each species and magnetic moment ( drift-coord -> gyro-coord )
      // OPTIM : We can skip foward transform after first call
      for(int s = NsLlD; s <= NsLuD; s++) { for(int m = NmLlD; m <= NmLuD; m++) {

           // Field has complicated stride, thus cannot be easily used in FFTSolver.
           // Use temporary Qm, and copy to Fields afterwards.

           // forward-transformation from drift-center -> gyro-center 
           gyroAverage((A4zz) Field0, (A4zz) Qm, m, s, true);

           // Copy result in temporary array to field
           [=] (CComplex Qm[Nq][NzLD][NkyLD][NxLD], CComplex Field[Nq][NsLD][NmLD][NzLB][NkyLD][NxLB+4]) 
           { 
                  Field[1:Nq][s][m][NzLlD:NzLD][NkyLlD:NkyLD][NxLlD:NxLD] 
                =    Qm[1:Nq]      [NzLlD:NzLD][NkyLlD:NkyLD][NxLlD:NxLD]   ;

           } ((A4zz) Qm, (A6zz) Field);
        
      } }
   
      updateBoundary();
   
   }   
  
   parallel->bcast(ArrayField.data(Field), DIR_V, ArrayField.getNum());
        
   return;

}


void Fields::calculateChargeDensity(const CComplex f0         [NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB],
                                    const CComplex f          [NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB],
                                    CComplex       Field0             [Nq][NzLD][NkyLD][NxLD]      ,
                                    const int m, const int s) 
{
 
   // In case of a full-f simulation the Maxwellian is subtracted

   const double pqnB_dvdm = M_PI * plasma->species[s].q * plasma->species[s].n0 * plasma->B0 * dv * 
                            (plasma->species[s].doGyro ? grid->dm[m] : 1.);
    
   omp_C2_for(int z=NzLlD; z<= NzLuD;z++) {  for(int y_k=NkyLlD; y_k<= NkyLuD; y_k++) { for(int x=NxLlD; x<= NxLuD;x++) {

              Field0[Q::rho][z][y_k][x] = ( __sec_reduce_add(f [s][m][z][y_k][x][NvLlD:NvLD]) 
                        - (plasma->global ? __sec_reduce_add(f0[s][m][z][y_k][x][NvLlD:NvLD]) : 0)) * pqnB_dvdm;
     
   } } } // z, y_k, x
   
   return; 
}



void Fields::calculateParallelCurrentDensity(const CComplex f0   [NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB],
                                             const CComplex f    [NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB],
                                             CComplex       Field0       [Nq][NzLD][NkyLD][NxLD]      ,
                                             const double V[NvGB], const int m, const int s          ) 
{
  
   const double qa_dvdm = plasma->species[s].q * plasma->species[s].alpha  * plasma->B0 * M_PI * dv * grid->dm[m] ;
   
   omp_C2_for(int z=NzLlD; z<= NzLuD;z++) {  for(int y_k=NkyLlD; y_k<= NkyLuD; y_k++) { for(int x=NxLlD; x<= NxLuD;x++) {

                Field0[Q::jp][z][y_k][x] = -__sec_reduce_add(V[NvLlD:NvLD] * f[s][m][z][y_k][x][NvLlD:NvLD]) * qa_dvdm;

   } } } // z, y_k, x

   return ; 

}

  
// Note : Did not checked correctness of this function
void Fields::calculatePerpendicularCurrentDensity(const CComplex f0     [NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB],
                                                  const CComplex f      [NsLD][NmLD][NzLB][NkyLD][NxLB][NvLB],
                                                  CComplex       Field0         [Nq][NzLD][NkyLD][NxLD]      ,
                                                  const double M[NmGB], const int m, const int s          ) 
{
   
   const double qan_dvdm = - plasma->species[s].q * plasma->species[s].alpha   * plasma->B0 * M_PI * dv * grid->dm[m] ;
   
   for(int z=NzLlD; z<= NzLuD;z++) { omp_for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int x=NxLlD; x<= NxLuD;x++){

      Field0[Q::jo][z][y_k][x] =  M[m] * __sec_reduce_add(f[s][m][z][y_k][x][NvLlD:NvLD]) * qan_dvdm;
            
   } } } // z, y_k, x

   return;
}

   
// @todo this needs some cleanup work  and retesting for 3-dimensional case.
// only support parallelized version (assuming we run always parallized).
void Fields::updateBoundary()
{
  updateBoundary((A6zz ) Field, 
                  (A6zz ) SendXl   , (A6zz) SendXu, (A6zz) RecvXl, (A6zz) RecvXu,
                  (A6zz ) SendZl   , (A6zz) SendZu, (A6zz) RecvZl, (A6zz) RecvZu);
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
   if(Nz > 1) {
      
    omp_for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int x=NxLlD; x<= NxLuD;x++) { 

        const CComplex a = ((CComplex) 0. + 1.j) *  (2.*M_PI * (2.* M_PI/Ly) * y_k);
      
        // NzLlD == NzGlD -> Connect only physcial boundaries after mode made one loop 
        SendZl[:][:][:][:][y_k-NkyLlD][x-NxLlD] = Field[1:Nq][NsLlD:NsLD][NmLlD:NmLD][NzLlD  :2][y_k][x] * cexp( ((NzLuD == NzGuD) ? a : 0.) * geo->nu(x));
        SendZu[:][:][:][:][y_k-NkyLlD][x-NxLlD] = Field[1:Nq][NsLlD:NsLD][NmLlD:NmLD][NzLuD-1:2][y_k][x] * cexp(-((NzLlD == NzGlD) ? a : 0.) * geo->nu(x));

     } }
   }
   
   // Exchange ghostcells between processors [ SendXu (CPU 1) ->RecvXl (CPU 2) ]
  parallel->updateBoundaryFields(Fields::SendXl, Fields::SendXu, Fields::RecvXl, Fields::RecvXu, ArrayBoundX.getNum(),
                                 Fields::SendZl, Fields::SendZu, Fields::RecvZl, Fields::RecvZu, ArrayBoundZ.getNum()); 

   // Back copy X-Boundary cell data
   Field[1:Nq][NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][NkyLlD:NkyLD][NxLlB-2:4] = RecvXl[:][:][:][:][:][:];
   Field[1:Nq][NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][NkyLlD:NkyLD][NxLuD+1:4] = RecvXu[:][:][:][:][:][:];
 
   // Back copy Z-Boundary cell data
   if( Nz > 1 ) {

   Field[1:Nq][NsLlD:NsLD][NmLlD:NmLD][NzLlB  :2][NkyLlD:NkyLD][NxLlD:NxLD] = RecvZl[:][:][:][:][:][:];
   Field[1:Nq][NsLlD:NsLD][NmLlD:NmLD][NzLuD+1:2][NkyLlD:NkyLD][NxLlD:NxLD] = RecvZu[:][:][:][:][:][:];
   
   }

   return; 

};





///////////////////////////////////////////// Data I/O ///////////////////////////////


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

void Fields::writeData(const Timing &timing, const double dt) 
{
   if (timing.check(dataOutputFields, dt)       )   {
      //FA_phi->write(&Field0[Fields::phi][NzLlD][NkyLlD][NzLlD]);
      //FA_Ap ->write(Field0[Fields::Ap ][NzLlD][NkyLlD][NzLlD]);
      //FA_Bp ->write(Field0[Fields::Bp ][NzLlD][NkyLlD][NzLlD]);
      FA_phiTime->write(&timing);
      
      parallel->print("Wrote fields data ... "); 
   }
} 
      
void Fields::closeData() 
{

   delete FA_phi;
   delete FA_Ap;
   delete FA_Bp;
   delete FA_phiTime;

}


void Fields::printOn(std::ostream &output) const 
{

         output   << "Poisson    |  " << "Base class" << std::endl;
         output   << "Ampere     |  " << ((plasma->nfields >= 2) ? "beta :  " + Setup::num2str(plasma->beta) : " --- no electromagnetic effects ---") << std::endl;
         output   << "B_parallel |  " << ((plasma->nfields >= 3) ? "beta :  " + Setup::num2str(plasma->beta) : " --- no electromagnetic effects ---") << std::endl;
         output   << "           |  Pert. : " << ((!(solveEq & Field::Ap) && (plasma->nfields >= 2))     ? ApPerturbationStr  : " ") << std::endl;
}


