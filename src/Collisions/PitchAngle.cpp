/*
 * =====================================================================================
 *
 *       Filename: PitchAngle.cpp
 *
 *    Description: Implementation of the Pitch angle scattering 
 *
 *         Author: Paul P. Hilscher (2012), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include "Collisions/PitchAngle.h"


Collisions_PitchAngle::Collisions_PitchAngle(Grid *grid, Parallel *parallel, Setup *setup, FileIO *fileIO, Geometry *geo) 
: Collisions(grid, parallel, setup, fileIO, geo) 
{

  nu    = setup->get("Vlasov.Collision.Nu", 0.e0);
  
  double lambda_C = setup->get("Collisions.CoulombLogarithm", 1.e-4);
     /* 
  double  e = 1.; 
  
  vc = M_PI * log(lambda_C) * pow4(e) * plasma->n_ref * plasma->L_ref / (cbrt(2.) * pow2(plasma->T_ref));

  pref_Cxy.resize(grid->RsGD); pref_Cxy = 0.;
  B.resize(grid->RsGD, Range(1,4)); B = 0.;
  
  // we need to calculate the pre-factors
  //
  //        for(int s=NsGlD; s <= NsGuD; s++)  pref_Cxy(s) += - vc * species(j).n * pow2(species(j).q)* species(j).T * pow(species[s].m, 3./2.) 

  // calculate B's for conservation term
        double v2 = pow2(V[v]) + M[m]; 
 
        B(s,1) = parallel->sum(sum(F0), DIR_ALL);
        B(s,1) = parallel->sum(sum(F0), DIR_ALL);
        :1
        B(s,1) = parallel->sum(sum(F0), DIR_ALL);
        B(s,1) = parallel->sum(sum(F0), DIR_ALL);

      * */ 
  
  // normalized velocity : double x = v/v_th;
}

void Collisions_PitchAngle::C_xy(const double df_dp, const int x, const int y, const int z, const int v, const int m, const int s) 
{
  /* 
  const double xc = 0.;
  
  return pref_Cxy(s) * 1./pow5(V[v]) * ( V[v] * F1(xc) + 3. * plasma->B0 * M[m] * df_dp / species[s].m * F2(df_dp));
  */
}


void Collisions_PitchAngle::C_vm(const double df_dv, const double df_dm, const double f, const int x, const int y, const int z, const int v, const int m, const int s) 
{
  /* 
  int j = s;
  
  const double xc = species[s].v_th / species(j).v_th * V[v];
  
  // parallel velocity
  const double C_v = ((2. * plasma->B0 * M[m])/ species[s].m * F1(xc) + pow2(V[v]) * F3(x)) * df_dv + (6. * M[m] * pow2(V[v]) * F2(xc)) * df_dm - 2. * F3(xc)/pow3(V[v]) * V[v] * f;
  const double C_m = ((6. * M[m] * V[v] * F2(xc)) * df_dv  + 4./plasma->B0 * pow2(V[v]) * M[m] * F1(xc)) * df_dv + (4. * pow2(M[m]) * F3(xc)) * df_dm - 2. * F3(xc)/pow3(V[v]) * 2. * M[m] * f;
        
  // BUG : this is wrong, we need to differentiate over v
  const double Cvm = (C_v + C_m) / dv;

  return Cvm;

  */
}

void Collisions_PitchAngle::solve(Fields *fields, const CComplex  *f, const CComplex *f0, CComplex *Coll, double dt, int rk_step)
{

  [=](const CComplex f   [NsLD][NmLD][NzLB][Nky][NxLB][NvLB],  // Phase-space function for current time step
      const CComplex f0  [NsLD][NmLD][NzLB][Nky][NxLB][NvLB],  // Background Maxwellian
            CComplex Coll[NsLD][NmLD][NzLB][Nky][NxLB][NvLB]   // Collisional term
     ) 
  {


  CComplex C_v[NmLD][NvLD];
  CComplex C_m[NmLD][NvLD];

  const double _kw_12_dv = 1./(12.*dv);

 // Take care m-space may be not equidistant discretized
  for(int s = NsLlD; s <= NsLuD; s++) {
  for(int z = NzLlD; z <= NzLuD; z++) { for(int y_k = 0     ; y_k <= Nky-2; y_k++) {
  for(int x = NxLlD; x <= NxLuD; x++) { 

  for(int m = NmLlD; m <= NmLuD; m++) { const double _kw_12_dm = 1./(12.*grid->dm[m]);
  for(int v = NvLlD; v <= NvLuD; v++) { 

    const CComplex df_dv = (8.*(f[s][m][z][y_k][x][v+1] - f[s][m][z][y_k][x][v-1]) - (f[s][m][z][y_k][x][v+2] - f[s][m][z][y_k][x][v-2]))/(12.*dv);
    const CComplex df_dm = (8.*(f[s][m+1][z][y_k][x][v] - f[s][m-1][z][y_k][x][v]) - (f[s][m+2][z][y_k][x][v] - f[s][m-2][z][y_k][x][v]))/(12.*dm);

    C_v[m-NmLlD][v-NvLlD] =   (2. * plasma->B0 * M[m]) * df_dv -   ( 4. * M[m] * V[v]                 ) * df_dm; 
    C_m[m-NmLlD][v-NvLlD] = - (4. * M[m] * V[v]      ) * df_dv +   ( 8./plasma->B0 * pow2(V[v]) * M[m]) * df_dm;
  } }

  for(int m = NmLlD; m <= NmLuD; m++) { for(int v = NvLlD; v <= NvLuD; v++) { 
        
      // use 4-th order central differences
      Coll[s][m][z][y_k][x][v] =  ( 8. * (C_v[v+1][m] +  C_v[v-1][m]) - (C_v[v+2][m] - C_v[v-2][m])) / (12.*dv) 
                                + ( 8. * (C_m[v][m+1] +  C_v[v][m-1]) - (C_v[v][m+1] - C_v[v][m-1])) / (12.*dm);
  } }

  } } } }
  
  } ((A6zz) f  , (A6zz) f0, (A6zz) Coll); 
}

void Collisions_PitchAngle::initData(hid_t fileID) 
{

}   
     

void Collisions_PitchAngle::printOn(std::ostream &output) const 
{
  output   << "Collisions |  Pitch-Angle " << ((nu > 0.) ? Setup::num2str(nu) :  "off") << std::endl;
}
