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

   doFieldCorrections = setup->get("Moments.FieldCorrections", 1);

}



void Moments::getMoments(const CComplex     f    [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                         const CComplex Field0[Nq][NzLD][Nky][NxLD],
                         CComplex Mom00[NsLD][NzLD][Nky][NxLD], 
                         CComplex Mom20[NsLD][NzLD][Nky][NxLD],
                         CComplex Mom02[NsLD][NzLD][Nky][NxLD])
{
  
  Mom20[:][:][:][:] = 0.;
  Mom02[:][:][:][:] = 0.;
  Mom00[:][:][:][:] = 0.;

  // temporary arrays for integration over \mu
  // Need Nq in order to let gyroaveraging work for Nq > 1
  CComplex Mom20_m[Nq][NzLD][Nky][NxLD], Mom02_m[Nq][NzLD][Nky][NxLD], 
           Mom00_m[Nq][NzLD][Nky][NxLD];

  for(int s = NsLlD; s <= NsLuD; s++) {
    
  // Scaling terms (1) 
  const double rho_L_ref = plasma->rho_ref / plasma->L_ref; 
  const double d_00 = rho_L_ref * plasma->n_ref * species[s].n0;
  const double d_20 = rho_L_ref * plasma->n_ref * species[s].n0 * pow2(plasma->c_ref * species[s].v_th);
  const double d_02 = rho_L_ref * plasma->n_ref * species[s].n0 * pow2(plasma->c_ref * species[s].v_th);
  
  // integrate over the first adiabat \mu
  for(int m = NmLlD; m <= NmLuD; m++) {
  
  const double d6Z_DK_00 = d_00 * M_PI * dv * grid->dm[m];
  const double d6Z_DK_20 = d_20 * M_PI * dv * grid->dm[m];
  const double d6Z_DK_02 = d_02 * M_PI * dv * grid->dm[m] * plasma->B0;
  
  // calculate drift-kinetic moments (in gyro-coordinates) 
  for(int z = NzLlD; z <= NzLuD; z++) {  for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { 
  for(int x = NxLlD; x <= NxLuD; x++) { 
     
    Mom00_m[0][z-NzLlD][y_k][x-NxLlD] = __sec_reduce_add(f[s][m][z][y_k][x][NvLlD:NvLD]);
    Mom20_m[0][z-NzLlD][y_k][x-NxLlD] = __sec_reduce_add(pow2(V[NvLlD:NvLD]) * f[s][m][z][y_k][x][NvLlD:NvLD]);
    Mom02_m[0][z-NzLlD][y_k][x-NxLlD] = __sec_reduce_add(M[m] * f[s][m][z][y_k][x][NvLlD:NvLD]);
      
  } } } // z, y_k, x
  
    // BUG : we may have decomposition in V

    // gyro-average Moments (push back to drift-coordinates)
    fields->gyroAverage(Mom00_m, Mom00_m, m, s, false, true);
    fields->gyroAverage(Mom20_m, Mom20_m, m, s, false, true);
    fields->gyroAverage(Mom02_m, Mom02_m, m, s, false, true);

    Mom00[s-NsLlD][:][:][:] += Mom00_m[0][:][:][:] * d6Z_DK_00;
    Mom20[s-NsLlD][:][:][:] += Mom20_m[0][:][:][:] * d6Z_DK_20;
    Mom02[s-NsLlD][:][:][:] += Mom02_m[0][:][:][:] * d6Z_DK_02;

  } // m

  parallel->reduce(&Mom00[s-NsLlD][0][0][0], Op::sum, DIR_M, NzLD * Nky * NxLD); 
  parallel->reduce(&Mom20[s-NsLlD][0][0][0], Op::sum, DIR_M, NzLD * Nky * NxLD); 
  parallel->reduce(&Mom02[s-NsLlD][0][0][0], Op::sum, DIR_M, NzLD * Nky * NxLD); 


  //////////////////////////////////////////////////////////////////////////
  // add field corrections (necessary for k_perp^2 > 1)
  if(doFieldCorrections) {
 

    // missing e-m corrections br_ref T_0s / B_02 j_0par / q v Y(a+1)
    const double Y_00 = 1.;
    const double Y_20 = 0.5;
    const double Y_02 = 1.;

    // pre-factors for field corrections (note : x-dependence of T neglected)
    const double d_FC_00 = d_00 * species[s].n0/species[s].T0 * Y_00;
    const double d_FC_20 = d_20 * species[s].n0/species[s].T0 * Y_20;
    const double d_FC_02 = d_02 * species[s].n0/species[s].T0 * Y_02;

    CComplex AAphi_m_0[Nq][NzLD][Nky][NxLD],  phi0[Nq][NzLD][Nky][NxLD],
             AAphi_m_1[Nq][NzLD][Nky][NxLD];


    phi0[0][:][:][:] = Field0[Field::phi][NzLlD:NzLD][:][NxLlD:NxLD];

    // calculate double gyro-average 
    fields->doubleGyroExp(phi0, AAphi_m_0, 0, s);
    fields->doubleGyroExp(phi0, AAphi_m_1, 1, s);

    // add field corrections from gyro-kinetic effects ( rhs missing b_0/T_ps^{b/2+1} * AAphi_
    Mom00[s-NsLlD][:][:][:] -= d_FC_00 * (species[s].q * (phi0[:][:][:] - AAphi_m_0[:][:][:]));
    Mom20[s-NsLlD][:][:][:] -= d_FC_20 * (species[s].q * (phi0[:][:][:] - AAphi_m_0[:][:][:]));
    Mom02[s-NsLlD][:][:][:] -= d_FC_02 * (species[s].q * (phi0[:][:][:] - AAphi_m_1[:][:][:]));

  } // doFieldCorrections


  } // s
    

  return;
}



void Moments::calculateMoment_Ap(const CComplex     f    [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                                const CComplex Field0[Nq][NzLD][Nky][NxLD],
                                CComplex Mom10[NsLD][NzLD][Nky][NxLD], 
                                CComplex Mom30[NsLD][NzLD][Nky][NxLD],
                                CComplex Mom12[NsLD][NzLD][Nky][NxLD])
{
  
  Mom10[:][:][:][:] = 0.;
  Mom30[:][:][:][:] = 0.;
  Mom12[:][:][:][:] = 0.;

  // temporary arrays for integration over \mu
  CComplex Mom10_m[Nq][NzLD][Nky][NxLD], Mom30_m[Nq][NzLD][Nky][NxLD], 
           Mom12_m[Nq][NzLD][Nky][NxLD]; 
  
  for(int s = NsLlD; s <= NsLuD; s++) {
  
  // integrate over the first adiabat \mu
  for(int m = NmLlD; m <= NmLuD; m++) {
  
  const double d6Z_00 = M_PI * species[s].n0 * species[s].T0 * plasma->B0 * dv * grid->dm[m] * grid->dXYZ;
  const double d6Z_20 = M_PI * species[s].n0 * species[s].T0 * plasma->B0 * dv * grid->dm[m] * grid->dXYZ *
                               species[s].v_th;
  const double d6Z_02 = M_PI * species[s].n0 * species[s].T0 * plasma->B0 * dv * grid->dm[m] * grid->dXYZ;

  // calculate drift-kinetic moments (in gyro-coordinates) 
  for(int z = NzLlD; z <= NzLuD; z++) {  for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { 
  for(int x = NxLlD; x <= NxLuD; x++) { 
     
    Mom10_m[0][z-NzLlD][y_k][x-NxLlD] = __sec_reduce_add(species[s].m * V[NvLlD:NvLD] * f[s][m][z][y_k][x][NvLlD:NvLD]);
    Mom30_m[0][z-NzLlD][y_k][x-NxLlD] = __sec_reduce_add(species[s].m * V[NvLlD:NvLD] * f[s][m][z][y_k][x][NvLlD:NvLD]);
    Mom12_m[0][z-NzLlD][y_k][x-NxLlD] = __sec_reduce_add(sqrt(M[m]) * f[s][m][z][y_k][x][NvLlD:NvLD]);
      
  } } } // z, y_k, x
  
    // BUG : we may have docomposition in V

    // gyro-average Moments (push back to drift-coordinates)
    fields->gyroAverage(Mom10_m, Mom10_m, m, s, false, true);
    fields->gyroAverage(Mom30_m, Mom30_m, m, s, false, true);
    fields->gyroAverage(Mom12_m, Mom12_m, m, s, false, true);

    Mom10[s-NsLlD][:][:][:] += Mom10_m[0][:][:][:] * d6Z_00;
    Mom30[s-NsLlD][:][:][:] += Mom30_m[0][:][:][:] * d6Z_20;
    Mom12[s-NsLlD][:][:][:] += Mom12_m[0][:][:][:] * d6Z_02;

  } // m
  
  parallel->reduce(&Mom10[s-NsLlD][0][0][0], Op::sum, DIR_M, NzLD * Nky * NxLD); 
  parallel->reduce(&Mom30[s-NsLlD][0][0][0], Op::sum, DIR_M, NzLD * Nky * NxLD); 
  parallel->reduce(&Mom12[s-NsLlD][0][0][0], Op::sum, DIR_M, NzLD * Nky * NxLD); 


  //////////////////////////////////////////////////////////////////////////
  // add field corrections (necessary for k_perp^2 > 1)
  if(doFieldCorrections) {

    CComplex AAphi_m_0[NzLD][Nky][NxLD],  phi0[NzLD][Nky][NxLD],
             AAphi_m_1[NzLD][Nky][NxLD];

    phi0[:][:][:] = Field0[Field::phi][NzLlD:NzLD][NkyLlD:Nky][NxLlD:NxLD];

    // calculate double gyro-average 
 //   fields->doubleGyroExp(phi0, AAphi_m_0, 0, s);
 //   fields->doubleGyroExp(phi0, AAphi_m_1, 1, s);

    // add field corrections from gyro-kinetic effects
    Mom10[s-NsLlD][:][:][:] -= species[s].q *       (phi0[:][:][:] - AAphi_m_0[:][:][:]);
    Mom30[s-NsLlD][:][:][:] -= species[s].q * 0.5 * (phi0[:][:][:] - AAphi_m_0[:][:][:]);
    Mom12[s-NsLlD][:][:][:] -= species[s].q *       (phi0[:][:][:] - AAphi_m_1[:][:][:]);

  } // doFieldCorrections


  } // s
    

  return;
}


void Moments::calculateMoment_Bp(const CComplex     f    [NsLD][NmLD][NzLB][Nky][NxLB  ][NvLB],
                                const CComplex Field0[Nq][NzLD][Nky][NxLD],
                                CComplex Mom02[NsLD][NzLD][Nky][NxLD], 
                                CComplex Mom22[NsLD][NzLD][Nky][NxLD],
                                CComplex Mom04[NsLD][NzLD][Nky][NxLD])
{
  
  Mom02[:][:][:][:] = 0.;
  Mom22[:][:][:][:] = 0.;
  Mom04[:][:][:][:] = 0.;

  // temporary arrays for integration over \mu
  CComplex Mom02_m[Nq][NzLD][Nky][NxLD], Mom22_m[Nq][NzLD][Nky][NxLD], 
           Mom04_m[Nq][NzLD][Nky][NxLD]; 
  
  for(int s = NsLlD; s <= NsLuD; s++) {
  
  // integrate over the first adiabat \mu
  for(int m = NmLlD; m <= NmLuD; m++) {
  
  const double d6Z_00 = M_PI * species[s].n0 * species[s].T0 * plasma->B0 * dv * grid->dm[m] * grid->dXYZ;
  const double d6Z_20 = M_PI * species[s].n0 * species[s].T0 * plasma->B0 * dv * grid->dm[m] * grid->dXYZ *
                               species[s].v_th;
  const double d6Z_02 = M_PI * species[s].n0 * species[s].T0 * plasma->B0 * dv * grid->dm[m] * grid->dXYZ;

  // calculate drift-kinetic moments (in gyro-coordinates) 
  for(int z = NzLlD; z <= NzLuD; z++) {  for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { 
  for(int x = NxLlD; x <= NxLuD; x++) { 
     
    Mom02_m[0][z-NzLlD][y_k][x-NxLlD] = __sec_reduce_add(species[s].m * V[NvLlD:NvLD] * f[s][m][z][y_k][x][NvLlD:NvLD]);
    Mom22_m[0][z-NzLlD][y_k][x-NxLlD] = __sec_reduce_add(species[s].m * V[NvLlD:NvLD] * f[s][m][z][y_k][x][NvLlD:NvLD]);
    Mom04_m[0][z-NzLlD][y_k][x-NxLlD] = __sec_reduce_add(sqrt(M[m]) * f[s][m][z][y_k][x][NvLlD:NvLD]);
      
  } } } // z, y_k, x
  
    // BUG : we may have docomposition in V

    // gyro-average Moments (push back to drift-coordinates)
    fields->gyroAverage(Mom02_m, Mom02_m, m, s, false, true);
    fields->gyroAverage(Mom22_m, Mom22_m, m, s, false, true);
    fields->gyroAverage(Mom04_m, Mom04_m, m, s, false, true);

    Mom02[s-NsLlD][:][:][:] += Mom02_m[0][:][:][:] * d6Z_00;
    Mom22[s-NsLlD][:][:][:] += Mom22_m[0][:][:][:] * d6Z_20;
    Mom04[s-NsLlD][:][:][:] += Mom04_m[0][:][:][:] * d6Z_02;

  } // m
  
  parallel->reduce(&Mom02[0][0][0][0], Op::sum, DIR_M, NsLD * NzLD * Nky * NxLD); 
  parallel->reduce(&Mom22[0][0][0][0], Op::sum, DIR_M, NsLD * NzLD * Nky * NxLD); 
  parallel->reduce(&Mom04[0][0][0][0], Op::sum, DIR_M, NsLD * NzLD * Nky * NxLD); 


  //////////////////////////////////////////////////////////////////////////
  // add field corrections (necessary for k_perp^2 > 1)
  if(doFieldCorrections) {

    CComplex AAphi_m_0[NzLD][Nky][NxLD],  phi0[NzLD][Nky][NxLD],
             AAphi_m_1[NzLD][Nky][NxLD];

    phi0[:][:][:] = Field0[Field::phi][NzLlD:NzLD][:][NxLlD:NxLD];

    // calculate double gyro-average 
//    fields->doubleGyroExp(phi0, AAphi_m_0, 0, s);
//    fields->doubleGyroExp(phi0, AAphi_m_1, 1, s);

    // add field corrections from gyro-kinetic effects
    Mom02[s-NsLlD][:][:][:] -= species[s].q *       (phi0[:][:][:] - AAphi_m_0[:][:][:]);
    Mom22[s-NsLlD][:][:][:] -= species[s].q * 0.5 * (phi0[:][:][:] - AAphi_m_0[:][:][:]);
    Mom04[s-NsLlD][:][:][:] -= species[s].q *       (phi0[:][:][:] - AAphi_m_1[:][:][:]);

  } // doFieldCorrections


  }
    

  return;
}


