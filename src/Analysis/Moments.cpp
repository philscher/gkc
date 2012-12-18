/*
 * =====================================================================================
 *
 *       Filename: Moments.cpp
 *
 *    Description: Calculation of the Moments of the Phase-space function
 *
 *         Author: Paul P. Hilscher (2012-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include "Analysis/Moments.h"
  
  
Moments::Moments(Setup *setup, Vlasov *_vlasov, Fields *_fields, Grid *_grid, Parallel *_parallel) 
: vlasov(_vlasov), fields(_fields), grid(_grid), parallel(_parallel) 
{

   doFieldCorrections = setup->get("Moments.includeFieldCorrections", 1);


};



void Moments::getMoments(const CComplex     f    [NsLD][NmLD][NzLB][NkyLD][NxLB  ][NvLB],
                         const CComplex Field0[Nq][NzLD][NkyLD][NxLD],
                         CComplex Mom00[NsLD][NzLD][NkyLD][NxLD], 
                         CComplex Mom20[NsLD][NzLD][NkyLD][NxLD],
                         CComplex Mom02[NsLD][NzLD][NkyLD][NxLD])
{
  
  Mom20[:][:][:][:] = 0.;
  Mom02[:][:][:][:] = 0.;
  Mom00[:][:][:][:] = 0.;

  // temporary arrays for integration over \mu
  CComplex Mom20_m[1][NzLD][NkyLD][NxLD], Mom02_m[1][NzLD][NkyLD][NxLD], 
           Mom00_m[1][NzLD][NkyLD][NxLD]; 


  for(int s = NsLlD; s <= NsLuD; s++) {
  
  // integrate over the first adiabat \mu
  for(int m = NmLlD; m <= NmLuD; m++) {
  
  const double d6Z_00 = M_PI * plasma->species[s].n0 * plasma->species[s].T0 * plasma->B0 * dv * grid->dm[m] * grid->dXYZ;
  const double d6Z_20 = M_PI * plasma->species[s].n0 * plasma->species[s].T0 * plasma->B0 * dv * grid->dm[m] * grid->dXYZ *
                               plasma->species[s].scale_v;
  const double d6Z_02 = M_PI * plasma->species[s].n0 * plasma->species[s].T0 * plasma->B0 * dv * grid->dm[m] * grid->dXYZ;

  // calculate drift-kinetic moments (in gyro-coordinates) 
  for(int z = NzLlD; z <= NzLuD; z++) {  for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { 
  for(int x = NxLlD; x <= NxLuD; x++) { 
     
    Mom00_m[0][z-NzLlD][y_k][x-NxLlD] = __sec_reduce_add(f[s][m][z][y_k][x][NvLlD:NvLD]);
    Mom20_m[0][z-NzLlD][y_k][x-NxLlD] = __sec_reduce_add(plasma->species[s].m * pow2(V[NvLlD:NvLD]) * f[s][m][z][y_k][x][NvLlD:NvLD]);
    Mom02_m[0][z-NzLlD][y_k][x-NxLlD] = __sec_reduce_add(M[m] * f[s][m][z][y_k][x][NvLlD:NvLD]);
      
  } } } // z, y_k, x

    // gyro-average Moments (push back to drift-coordinates)
    fields->gyroAverage(Mom00_m, Mom00_m, m, s, false, true);
    fields->gyroAverage(Mom20_m, Mom20_m, m, s, false, true);
    fields->gyroAverage(Mom02_m, Mom02_m, m, s, false, true);

    Mom00[s-NsLlD][:][:][:] += Mom00_m[0][:][:][:] * d6Z_00;
    Mom20[s-NsLlD][:][:][:] += Mom20_m[0][:][:][:] * d6Z_20;
    Mom02[s-NsLlD][:][:][:] += Mom02_m[0][:][:][:] * d6Z_02;

  } // m
  
  parallel->reduce(&Mom00[0][0][0][0], Op::sum, DIR_M, NsLD * NzLD * NkyLD * NxLD); 
  parallel->reduce(&Mom20[0][0][0][0], Op::sum, DIR_M, NsLD * NzLD * NkyLD * NxLD); 
  parallel->reduce(&Mom02[0][0][0][0], Op::sum, DIR_M, NsLD * NzLD * NkyLD * NxLD); 


  //////////////////////////////////////////////////////////////////////////
  // add field corrections (necessary for k_perp^2 > 1)
  
  if(doFieldCorrections) {

    CComplex AAphi_m_0[NzLD][NkyLD][NxLD],  phi0[NzLD][NkyLD][NxLD],
             AAphi_m_1[NzLD][NkyLD][NxLD];

    phi0[:][:][:] = Field0[Field::phi][NzLlD:NzLD][NkyLlD:NkyLD][NxLlD:NxLD];

    // calculate double gyro-average 
    fields->doubleGyroExp(phi0, AAphi_m_0, 0, s);
    fields->doubleGyroExp(phi0, AAphi_m_1, 1, s);

    // add field corrections from gyro-kinetic effects
    Mom00[s-NsLlD][:][:][:] -= plasma->species[s].q *       (phi0[:][:][:] - AAphi_m_0[:][:][:]);
    Mom20[s-NsLlD][:][:][:] -= plasma->species[s].q * 0.5 * (phi0[:][:][:] - AAphi_m_0[:][:][:]);
    Mom02[s-NsLlD][:][:][:] -= plasma->species[s].q *       (phi0[:][:][:] - AAphi_m_1[:][:][:]);

  } // doFieldCorrections

  }
    

  return;
}


/* 
void Moments::calculateMoment()
{
  // add field corrections (necessary for k_perp^2 > 1) 
  CComplex AAphi_m_0[NzLD][NkyLD][NxLD],  phi0[NzLD][NkyLD][NxLD],
           AAphi_m_1[NzLD][NkyLD][NxLD];

  phi0[:][:][:] = Field0[Field::phi][NzLlD:NzLD][NkyLD:NkyLD][NxLlD:NxLD];

  // calculate double gyro-average 
  fields->doubleGyroExp(phi0, AAphi_m_0, 0, s);
  fields->doubleGyroExp(phi0, AAphi_m_1, 1, s);

  // calculate seperate moments   
  calculateMoment00(M_00);


  Mom20[:][:][:] -= plasma->species[s].q * 0.5 * (phi0[:][:][:] - AAphi_m_0[:][:][:]);
  Mom02[:][:][:] -= plasma->species[s].q *       (phi0[:][:][:] - AAphi_m_1[:][:][:]);

};

void Moments::calculateMoment00(Mom00, phi0, AAphi_m_0, bool addCorrections) 
{
  CComplex Mom00_m[NzLD][NkyLD][NxLD];

  Mom_00[:][:][:] = 0.;

  // integrate over the first adiabat \mu
  for(int m = NmLlD; m <= NmLuD; m++) {
  
  const double d6Z = M_PI * plasma->species[s].n0 * plasma->species[s].T0 * plasma->B0 * dv * grid->dm[m] * grid->dXYZ;

  // calculate (drift-kinetic) moments  
  for(int z = NzLlD; z <= NzLuD; z++) {  for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { 
  for(int x = NxLlD; x <= NxLuD; x++) { 
     
    Mom00_m[0][z-NzLlD][y_k][x-NxLlD] = __sec_reduce_add(f[s][m][z][y_k][x][NvLlD:NvLD]);
      
  } } }
    
  fields->gyroAverage(Mom00_m, Mom00_m, m, s, false, true);
    
    Mom00[:][:][:] += Mom00_m[0][:][:][:] * d6Z;
  
  }

  // add field corrections
  if(addCorrections) Mom00[:][:][:] -= plasma->species[s].q * (phi0[:][:][:] - AAphi_m_0[:][:][:]);
 
}


void Moments::calculateMoment20() 
{







}


void Momnets::calculateMoment02() 
{







}
*/
