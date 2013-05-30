/*
 * =====================================================================================
 *
 *       Filename: Auxiliary.cpp
 *
 *    Description: Auxiliary diagnostic functions for gkc
 *
 *         Author: Paul P. Hilscher (2013), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include "Auxiliary.h"

Auxiliary::Auxiliary(Setup *setup, FileIO *fileIO, Parallel *_parallel, Fields *_fields, Vlasov *_vlasov, Grid *grid, FFTSolver *_fft)
: vlasov(_vlasov), fields(_fields), fft(_fft), parallel(_parallel)
{
  
  zonalFlow     = setup->get("Aux.ZonalFlow", 0); 
  dataOutputZF  = Timing( setup->get("DataOutput.IslandZF.Step", -1), setup->get("DataOutput.IslandZF.Time", -1.));

  ArrayZF   = nct::allocate(nct::Range(0, 3), grid->RsLD, grid->RkyLD, grid->RxLD)(&ZFProd);
  ArrayZFGy = nct::allocate(nct::Range(0,Nq), grid->RzLD, grid->RkyLD, grid->RxLD)(&ZF_Gyro_In, &ZF_Gyro_Out);

  // initialize data 
  hsize_t dim  [] = { 3, Ns     , Nky, Nx     , 1 };
  hsize_t mdim [] = { 3, Ns     , Nky, Nx     , H5S_UNLIMITED };
  hsize_t cdim [] = { 3, NsLD   , Nky, NxLD   ,             1 };
  hsize_t moff [] = { 0, 0      , 0  , 0      ,             0 };
  hsize_t off  [] = { 0, NsLlB-1, 0  , NxLlB-1,             0 }; 
    
  bool write = (parallel->Coord[DIR_VM] == 0);
   
  auxGroup = fileIO->newGroup("Aux");

  FA_ZFProd      = new FileAttr("ZonalFlow", auxGroup, fileIO->file, 5, dim, mdim, cdim, moff, cdim, off, write, fileIO->complex_tid);
  FA_ZFProdTime  = fileIO->newTiming(auxGroup);
}

Auxiliary::~Auxiliary()
{
  delete FA_ZFProd;
  delete FA_ZFProdTime;
   
  H5Gclose(auxGroup);
}

void Auxiliary::calculateZonalFlowProduction(const double dt)
{
  const double _kw_12_dx = 1./(12.*dx);

  [=](const CComplex fs [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
      const CComplex f0 [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
      const CComplex Fields[Nq][NsLD][NmLD][NzLB][Nky][NxLB+4],
      CComplex      ZFProd[3][NsLD][Nky][NxLD],
      CComplex ZF_In[Nq][NzLD][Nky][NxLD], CComplex ZF_Out[Nq][NzLD][Nky][NxLD])
  {
   
    CComplex ZF_T[NzLD][Nky][NxLD][3]; 

  for(int s = NsLlD; s <= NsLuD; s++) { 
    
    const double alpha = species[s].alpha;
    const double sigma = species[s].sigma;
    
  for(int m = NmLlD; m <= NmLuD; m++) {  
  
    ZF_T[:][:][:][:] = 0.;
    
  for(int z = NzLlD; z <= NzLuD; z++) { 
  
  //#pragma omp for
  for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) {

    const CComplex ky = _imag *  std::abs(fft->ky(y_k));

  for(int x   = NxLlD ; x   <= NxLuD ;   x++) {  
  
    // Calculate Poisson bracket by hand, no need for anti-aliasing etc
    const CComplex dphi_dx =  8.* (Fields[Field::phi][s][m][z][y_k][x+1] - Fields[Field::phi][s][m][z][y_k][x-1]) 
                                - (Fields[Field::phi][s][m][z][y_k][x+2] - Fields[Field::phi][s][m][z][y_k][x-2]) * _kw_12_dx;
  
    const CComplex dAp_dx  =  8.* (Fields[Field::Ap ][s][m][z][y_k][x+1] - Fields[Field::Ap ][s][m][z][y_k][x-1]) 
                                - (Fields[Field::Ap ][s][m][z][y_k][x+2] - Fields[Field::Ap ][s][m][z][y_k][x-2]) * _kw_12_dx;
  
    const CComplex phi     = Fields[Field::phi][s][m][z][y_k][x];
    const CComplex Ap      = Fields[Field::Ap ][s][m][z][y_k][x];
      
  for(int v = NvLlD; v <= NvLuD; v++) {
                           
    const CComplex f0_     = f0[s][m][z][y_k][x][v];
    const CComplex f1      = fs[s][m][z][y_k][x][v];
          
    const CComplex dfs_dx  = 8.*(fs[s][m][z][y_k][x+1][v] - fs[s][m][z][y_k][x-1][v])  
                              - (fs[s][m][z][y_k][x+2][v] - fs[s][m][z][y_k][x-2][v]) * _kw_12_dx;


    // [\phi, f_1]
    ZF_T[z-NzLlD][y_k][x-NxLlD][0] +=    (    dphi_dx  * conj(ky * f1) -     (ky * phi) * conj(dfs_dx)
                                       + conj(dphi_dx) *     (ky * f1) - conj(ky * phi) *      dfs_dx) * dv;
    if(Nq > 1) 
    {
    // v [A_\parallel, f_1]
    ZF_T[z-NzLlD][y_k][x-NxLlD][1] +=  - alpha * V[v]  *
                                       (      dAp_dx   * conj(ky * f1) -     (ky * Ap) * conj(dfs_dx)
                                       + conj(dAp_dx)  *     (ky * f1) - conj(ky * Ap) *      dfs_dx) * dv;
         
    // v [A_\parallel, \phi_1]
    ZF_T[z-NzLlD][y_k][x-NxLlD][2] +=  - alpha * V[v]  * sigma * f0_ * 
                                       (      dAp_dx   * conj(ky * phi) -     (ky * Ap) * conj(dphi_dx)
                                       + conj(dAp_dx)  *     (ky * phi) - conj(ky * Ap) *      dphi_dx) * dv;
    }
  } } } } // v, x, y_k, z
  
  // Correct for gyro-averaging
  for(int n = 0; n < 3; n++) {
  
    ZF_In[0][NzLlD:NzLD][:][NxLlD:NxLD] = ZF_T[:][:][:][n];
    fields->gyroAverage(ZF_In, ZF_Out, m, s, true);
    ZFProd[n][s][:][NxLlD:NxLD] += ZF_Out[0][NzLlD][:][NxLlD:NxLD];
  }
 
  } } // m, s
  }( (A6zz) vlasov->f, (A6zz) vlasov->f0, (A6zz) fields->Field, (A4zz) ZFProd,
     (A4zz) ZF_Gyro_In, (A4zz) ZF_Gyro_Out);
}

void Auxiliary::writeData(const Timing &timing, const double dt)
{

  if (timing.check(dataOutputZF, dt)       )   {

  if(zonalFlow) {
    calculateZonalFlowProduction(dt); 
    
    parallel->reduce(ArrayZF.data(ZFProd), Op::sum, DIR_VM,  ArrayZF.getNum()); 
    FA_ZFProd->write(ArrayZF.data((CComplex *) ZFProd));
    FA_ZFProdTime->write(&timing);
     
    // reset to zero
    [] (CComplex      ZFProd[3][NsLD][Nky][NxLD]) {

      ZFProd[:][NsLlD:NsLD][:][NxLlD:NxLD] = 0.;

    }( (A4zz) ZFProd);
  }
  }
}

void Auxiliary::printOn(std::ostream &output) const
{ 
  output   << "           |  Zonal flow production : " << (zonalFlow ? "Yes" : "No") << std::endl;
}

