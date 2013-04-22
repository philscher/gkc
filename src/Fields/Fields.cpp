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

// global parameters defined in Global.h
int GC2, GC4;

Fields::Fields(Setup *setup, Grid *_grid, Parallel *_parallel, FileIO *fileIO, Geometry *_geo)  : 
grid(_grid), parallel(_parallel), geo(_geo), solveEq(0)

{
  GC2 = 2; GC4 = 4, Nq=plasma->nfields;
 
  // Note : Fixed fields are initialized in Init.cpp
  solveEq |=  ((Nq >= 1) && (setup->get("Init.FixedPhi", ".0") == ".0")) ? Field::Iphi : 0;
  solveEq |=  ((Nq >= 2) && (setup->get("Init.FixedAp" , ".0") == ".0")) ? Field::IAp  : 0;
  solveEq |=  ((Nq >= 3) && (setup->get("Init.FixedBp" , ".0") == ".0")) ? Field::IBp  : 0;

  // Allocate variables for calculating source terms and field terms
  ArrayField  = nct::allocate(nct::Range(0,Nq), grid->RsLD, grid->RmLD , grid->RzLB, grid->RkyLD, nct::Range(NxLlD-4, NxLD+8))(&Field);
  ArrayField0 = nct::allocate(nct::Range(0,Nq), grid->RzLD, grid->RkyLD, grid->RxLD)(&Q, &Qm, &Field0);

  // Allocate boundary conditions, allocate Send/Recv buffers, note we have 4 ghost cells for X, 0 for Y
  ArrayBoundX = nct::allocate(nct::Range(0,    4 * Nky * NzLD * NmLD * NsLD * Nq))(&SendXl, &SendXu, &RecvXl, &RecvXu);
  ArrayBoundZ = nct::allocate(nct::Range(0, NxLD * Nky *    2 * NmLD * NsLD * Nq))(&SendZl, &SendZu, &RecvZl, &RecvZu);
       
  // should equal 1/2. To avoid cancellation error due to numerical errors, calculate it numerically, see Dannert[2] 
  Yeb = (1./sqrt(M_PI) * __sec_reduce_add(pow2(V[NvLlD:NvLD]) * exp(-pow2(V[NvLlD:NvLD]))) * dv) * geo->eps_hat * plasma->beta; 

  initData(setup, fileIO);

} 

Fields::~Fields() 
{
  closeData() ;

}


void Fields::solve(const CComplex *f0, CComplex *f, Timing timing)
{
  // calculate source terms  Q (Q is overwritten in the first iteration )
  for(int s = NsLlD, loop=0; s <= NsLuD; s++) { for(int m = NmLlD; m <= NmLuD; m++, loop++) {

    // calculate drift-kinetic terms (use OpenMP schedule(static) to allow cache fusion)
    if(solveEq & Field::Iphi) calculateChargeDensity               ((A6zz) f0, (A6zz) f, (A4zz) Field0, m, s);
    if(solveEq & Field::IAp ) calculateParallelCurrentDensity      ((A6zz) f0, (A6zz) f, (A4zz) Field0, m, s);
    if(solveEq & Field::IBp ) calculatePerpendicularCurrentDensity ((A6zz) f0, (A6zz) f, (A4zz) Field0, m, s);
    #pragma omp barrier

    // Integrate over velocity space through different CPU's
    #pragma omp single
    parallel->reduce(ArrayField0.data(Field0), Op::sum, DIR_V, ArrayField0.getNum(), false); 
      
    // OPTIM : Normally we would decompose in m&s, thus no need for Qm                         
    // backward-transformation from gyro-center to drift-center 
    if (parallel->Coord[DIR_V] == 0) {

      gyroAverage((A4zz) Field0, (A4zz) Qm, m, s, false);           

      // Lambda function to integrate over m ( for source terms), If loop=0, we overwrite value of Q
      [=] (CComplex Q[Nq][NzLD][Nky][NxLD], CComplex Qm[Nq][NzLD][Nky][NxLD])
      {
        #pragma omp for collapse(2)
        for(int z = NzLlD; z <= NzLuD; z++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { 

          if  (loop == 0) Q[0:Nq][z][y_k][NxLlD:NxLD]   = Qm[0:Nq][z][y_k][NxLlD:NxLD];
          else            Q[0:Nq][z][y_k][NxLlD:NxLD]  += Qm[0:Nq][z][y_k][NxLlD:NxLD];

        } }

      } ((A4zz) Q, (A4zz) Qm);
    }
      
  } } // for m, s

  /////////////////////////////// Solve for the corresponding fields ////////////////////////////////
  // Note :  Fields are only solved by root nodes  (X=0, V=0, S=0), Gyro-averaging is done for (V=0)
  if(parallel->Coord[DIR_V] == 0) {

    // integrate over mu-space and over species
    #pragma omp single
    parallel->reduce(ArrayField0.data(Q), Op::sum, DIR_MS, ArrayField0.getNum(), false); 

    // Solve field equation in drift coordinates
    // This routine is solved only on root nodes, thus efficiency is crucial for scalability
    // Any chance to make some other useful work here ?
    if(parallel->Coord[DIR_MS] == 0) solveFieldEquations((A4zz) Q, (A4zz) Field0);

    #pragma omp single
    parallel->bcast(ArrayField0.data(Field0), DIR_MS, ArrayField0.getNum());

    // Gyro-averaging procedure for each species and magnetic moment ( drift-coord -> gyro-coord )
    // OPTIM : We can skip forward transform after first call
    for(int s = NsLlD; s <= NsLuD; s++) { for(int m = NmLlD; m <= NmLuD; m++) {

      // Fields has a complicated stride, thus cannot be easily be used in FFTSolver.
      // Use temporary Qm, and copy to Fields afterwards.

      // forward-transformation from drift-center -> gyro-center 
      gyroAverage((A4zz) Field0, (A4zz) Qm, m, s, true);

      // Copy result in temporary array to field
      [=] (CComplex Qm[Nq][NzLD][Nky][NxLD], CComplex Field[Nq][NsLD][NmLD][NzLB][Nky][NxLB+4]) 
      { 
        #pragma omp for collapse(2)
        for(int z = NzLlD; z <= NzLuD; z++) {  for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) {

          Field[:][s][m][z][y_k][NxLlD:NxLD] = Qm[:][z][y_k][NxLlD:NxLD];

        } }

      } ((A4zz) Qm, (A6zz) Field);
        
    } } // s, m
 
    #pragma omp single
    updateBoundary();

  }  // (parallel->Coord[DIR_V] == 0)

  #pragma omp single
  parallel->bcast(ArrayField.data(Field), DIR_V, ArrayField.getNum());
   
  return;

}


void Fields::calculateChargeDensity(const CComplex f0      [NsLD][NmLD][NzLB][Nky][NxLB][NvLB],
                                    const CComplex f       [NsLD][NmLD][NzLB][Nky][NxLB][NvLB],
                                          CComplex Field0          [Nq][NzLD][Nky][NxLD]      ,
                                    const int m, const int s) 
{
  
  // In case of a full-f simulation the Maxwellian is subtracted
  const double pqnB_dvdm = M_PI * species[s].q * species[s].n0 * plasma->B0 * dv * grid->dm[m] ;

  #pragma omp for collapse(2) nowait, schedule(static)
  for(int z = NzLlD; z <= NzLuD; z++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { 
  for(int x = NxLlD; x <= NxLuD; x++) {

    Field0[Field::phi][z][y_k][x] = ( __sec_reduce_add(f [s][m][z][y_k][x][NvLlD:NvLD]) 
                  - (plasma->global ? __sec_reduce_add(f0[s][m][z][y_k][x][NvLlD:NvLD]) : 0)) * pqnB_dvdm;
     
  } } } // z, y_k, x
   
  return; 

}



void Fields::calculateParallelCurrentDensity(const CComplex f0    [NsLD][NmLD][NzLB][Nky][NxLB][NvLB],
                                             const CComplex f     [NsLD][NmLD][NzLB][Nky][NxLB][NvLB],
                                                   CComplex Field0        [Nq][NzLD][Nky][NxLD]      ,
                                             const int m, const int s                                  ) 
{
  const double qa_dvdm = species[s].q * species[s].n0  * species[s].alpha * M_PI * dv * grid->dm[m] ;
   
  #pragma omp for collapse(2) nowait, schedule(static)
  for(int z = NzLlD; z <= NzLuD; z++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { 
  for(int x = NxLlD; x <= NxLuD; x++) {

    Field0[Field::Ap][z][y_k][x] = -__sec_reduce_add(V[NvLlD:NvLD] * f[s][m][z][y_k][x][NvLlD:NvLD]) * qa_dvdm;

  } } } // z, y_k, x

  return;

}

  
// Note : Did not checked correctness of this function
void Fields::calculatePerpendicularCurrentDensity(const CComplex f0     [NsLD][NmLD][NzLB][Nky][NxLB][NvLB],
                                                  const CComplex f      [NsLD][NmLD][NzLB][Nky][NxLB][NvLB],
                                                  CComplex       Field0         [Nq][NzLD][Nky][NxLD]      ,
                                                  const int m, const int s                                    ) 
{
   
  const double qan_dvdm = - species[s].q * species[s].alpha * plasma->B0 * M_PI * dv * grid->dm[m] ;
   
  #pragma omp for collapse(2) nowait, schedule(static)
  for(int z = NzLlD; z <= NzLuD; z++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { 
  for(int x = NxLlD; x <= NxLuD; x++) {

    Field0[Field::Bp][z][y_k][x] =  M[m] * __sec_reduce_add(f[s][m][z][y_k][x][NvLlD:NvLD]) * qan_dvdm;
            
  } } } // z, y_k, x

  return;
}

   
void Fields::updateBoundary()
{
 
  [&](
         CComplex Field [Nq][NsLD][NmLD][NzLB][Nky][NxLB+4], 
         CComplex SendXl[Nq][NsLD][NmLD][NzLD][Nky][GC4   ], CComplex SendXu[Nq][NsLD][NmLD][NzLD][Nky][GC4 ],
         CComplex RecvXl[Nq][NsLD][NmLD][NzLD][Nky][GC4   ], CComplex RecvXu[Nq][NsLD][NmLD][NzLD][Nky][GC4 ],
         CComplex SendZl[Nq][NsLD][NmLD][GC2 ][Nky][NxLD  ], CComplex SendZu[Nq][NsLD][NmLD][GC2 ][Nky][NxLD],
         CComplex RecvZl[Nq][NsLD][NmLD][GC2 ][Nky][NxLD  ], CComplex RecvZu[Nq][NsLD][NmLD][GC2 ][Nky][NxLD])
  {
   
    // X-Boundary (we have extended BC  - 4 ghost cells)
    SendXl[:][:][:][:][:][:] = Field[0:Nq][NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][:][NxLlD  :4];
    SendXu[:][:][:][:][:][:] = Field[0:Nq][NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][:][NxLuD-3:4];

    // Boundaries not required in Y (as we use Fourier modes)

    // (For two-dimensional simulations (x,y) boundaries are not required ) thus Nz > 1
    for(int y_k = NkyLlD; y_k <= NkyLuD && (Nz > 1); y_k++) { for(int x = NxLlD; x <= NxLuD; x++) { 

      // Z-Boundary (we need to connect the magnetic field lines)
      // NzLlD == NzGlD -> Connect only physical boundaries after mode made one loop 
      const CComplex shift_l = (NzLlD == NzGlD) ? cexp(-_imag * 2.* M_PI * geo->nu(x)) : 1.;
      const CComplex shift_u = (NzLuD == NzGuD) ? cexp(+_imag * 2.* M_PI * geo->nu(x)) : 1.;
      
      SendZl[:][:][:][:][y_k][x-NxLlD] = Field[0:Nq][NsLlD:NsLD][NmLlD:NmLD][NzLlD  :2][y_k][x] * shift_l;
      SendZu[:][:][:][:][y_k][x-NxLlD] = Field[0:Nq][NsLlD:NsLD][NmLlD:NmLD][NzLuD-1:2][y_k][x] * shift_u;

    } }
   
    // Exchange ghost cells between processors [ SendXu (CPU 1) -> RecvXl (CPU 2) ]
    parallel->updateBoundaryFields(Fields::SendXl, Fields::SendXu, Fields::RecvXl, Fields::RecvXu, ArrayBoundX.getNum(),
                                   Fields::SendZl, Fields::SendZu, Fields::RecvZl, Fields::RecvZu, ArrayBoundZ.getNum()); 

    // Back copy X-Boundary cell data
    Field[0:Nq][NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][:][NxLlB-2:4] = RecvXl[:][:][:][:][:][:];
    Field[0:Nq][NsLlD:NsLD][NmLlD:NmLD][NzLlD:NzLD][:][NxLuD+1:4] = RecvXu[:][:][:][:][:][:];
 
    // Back copy Z-Boundary cell data
    if( Nz > 1 ) {
      Field[0:Nq][NsLlD:NsLD][NmLlD:NmLD][NzLlB  :2][:][NxLlD:NxLD] = RecvZl[:][:][:][:][:][:];
      Field[0:Nq][NsLlD:NsLD][NmLlD:NmLD][NzLuD+1:2][:][NxLlD:NxLD] = RecvZu[:][:][:][:][:][:];
    }
  
    }((A6zz ) Field, 
      (A6zz ) SendXl, (A6zz) SendXu, (A6zz) RecvXl, (A6zz) RecvXu,
      (A6zz ) SendZl, (A6zz) SendZu, (A6zz) RecvZl, (A6zz) RecvZu);
  
  return; 

};





///////////////////////////////////////////// Data I/O ///////////////////////////////


void Fields::initData(Setup *setup, FileIO *fileIO) 
{
    
  // Set sizes : Note, we use Fortran ordering for field variables 
  hsize_t dim  [] = { Nq, Nz     , Nky, Nx     ,             1 };
  hsize_t mdim [] = { Nq, Nz     , Nky, Nx     , H5S_UNLIMITED };
  hsize_t cBdim[] = { Nq, NzLD   , Nky, NxLD   ,             1 };
  hsize_t cdim [] = { Nq, NzLD   , Nky, NxLD   ,             1 };
  hsize_t moff [] = { 0 , 0      , 0  , 0      ,             0 };
  hsize_t off  [] = { 0 , NzLlB-1, 0  , NxLlB-1,             0 }; 
     
  bool write = (parallel->Coord[DIR_VMS] == 0);
     
  hid_t fieldsGroup = fileIO->newGroup("Fields");
     
  FA_fields      = new FileAttr("Phi" , fieldsGroup, fileIO->file, 5, dim, mdim, cdim, moff, cBdim, off, write, fileIO->complex_tid);
  FA_fieldsTime  = fileIO->newTiming(fieldsGroup);
        
  H5Gclose(fieldsGroup);
      
  dataOutputFields = Timing(setup->get("DataOutput.Fields.Step", -1), setup->get("DataOutput.Fields.Time", -1.));

}   

void Fields::writeData(const Timing &timing, const double dt) 
{
  if (timing.check(dataOutputFields, dt)       )   {
    FA_fields->write(ArrayField0.data(Field0));
    FA_fieldsTime->write(&timing);
      
    parallel->print("Wrote fields data ... "); 
  }
} 
      
void Fields::closeData() 
{
   delete FA_fields;
   delete FA_fieldsTime;
}


void Fields::printOn(std::ostream &output) const 
{
  output   << "Fields     |  λD : "  << sqrt(plasma->debye2) 
                       << "   β : "  << Setup::num2str(plasma->beta)  
                       << "  A∥ : " << (Nq >= 2 ? "-on-" : "-off-")
                       << "  B∥ : " << (Nq >= 3 ? "-on-" : "-off-")
                       << std::endl;

}

