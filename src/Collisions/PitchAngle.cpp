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

#ifndef COLLISIONS__PITCHANGLE_H
#define COLLISIONS__PITCHANGLE_H

#include "Collisions/PitchAngle.h"

void Collision_PitchAngle::solve(Fields *fields, const CComplex  *f, const CComplex *f0, CComplex *Coll, double dt, int rk_step)
{

  Complex df_dm[NmLD];
  Complex df_dv[NmLD];

  for(int s = NsLlD; s <= NsLuD; s++) { for(int m  = NmLlD ;   m <= NmLuD ;   m++) { 
  for(int z = NzLlD; z <= NzLuD; z++) { for(int y_k= NkyLlD; y_k <= NkyLuD; y_k++) {
  for(int x = NxLlD; x <= NxLuD; x++) { for(int v  = NvLlD ;   v <= NvLuD ;   v++) { 


    const Complex df_dv = (8.  *(f[s][m][z][y_k][x][v+1] - f[s][m][z][y_k][x][v-1]) - (f[s][m][z][y_k][x][v+2] - f[s][m][z][y_k][x][v-2]))/(12.*dv);
    df_dm[v][m] = (8.  *(f[s][m+1][z][y_k][x][v] - f[s][m-1][z][y_k][x][v]) - (f[s][m+2][z][y_k][x][v] - f[s][m-2][z][y_k][x][v]))/(12.*dm);


    
    C_v[v][m] =   (2. * plasma->B0 * M[m]) * df_dv -   ( 4. * M[m] * V[v]                 ) * df_dm 
    C_m[v][m] = - (4. * M[m] * V[v]      ) * df_dv +   ( 8./plasma->B0 * pow2(V[v]) * M[m]) * df_dm
  } }

  for(int v = NvLlD; v <= NvLuD; v++) { for(int m  = NmLlD ; m <= NmLuD ; m++) { 
        
      // use 4-th order central differences
      df_dt += ( 8. * (C_v[v+1][m] +  C_v[v-1][m]) - (C_v[v+2][m] - C_v[v-2][m])) / (12.*dv);
      df_dt += ( 8. * (C_m[v][m+1] +  C_v[v][m-1]) - (C_v[v][m+1] - C_v[v][m-1])) / (12.*dm);
      
  }} 
  } } } }
  
}


Collision_PitchAngle::Collision_PitchAngle(Setup *setup,  Grid *grid) 
{

  nu    = setup->get("Vlasov.Collision.Nu", 0.e0);
  
  double lambda_C = setup->get("Collisions.CoulombLogarithm", 1.e-4);
  double  e = 1.; 
  
  vc = M_PI * log(lambda_C) * pow4(e) * plasma->n_ref * plasma->L_ref / (cbrt(2.) * pow2(plasma->T_ref));

  pref_Cxy.resize(grid->RsGD); pref_Cxy = 0.;
  B.resize(grid->RsGD, Range(1,4)); B = 0.;
  
  // we need to calculate the pre-factors
  //
  //        for(int s=NsGlD; s <= NsGuD; s++)  pref_Cxy(s) += - vc * species(j).n * pow2(species(j).q)* species(j).T * pow(species[s].m, 3./2.) 

  // calculate B's for conservation term
     /* 
        double v2 = pow2(V[v]) + M[m]; 
 
        B(s,1) = parallel->sum(sum(F0), DIR_ALL);
        B(s,1) = parallel->sum(sum(F0), DIR_ALL);
        :1
        B(s,1) = parallel->sum(sum(F0), DIR_ALL);
        B(s,1) = parallel->sum(sum(F0), DIR_ALL);

      * */ 
  
  // normalized velocity : double x = v/v_th;
}

double Collision_PitchAngle::C_xy(const double df_dp, const int x, const int y, const int z, const int v, const int m, const int s) 
{
  const double xc = 0.;
  
  return pref_Cxy(s) * 1./pow5(V[v]) * ( V[v] * F1(xc) + 3. * plasma->B0 * M[m] * df_dp / species[s].m * F2(df_dp));
}


inline double C_vm(const double df_dv, const double df_dm, const double f, const int x, const int y, const int z, const int v, const int m, const int s) 
{
  int j = s;
  
  const double xc = species[s].v_th / species(j).v_th * V[v];
  
  // parallel velocity
  const double C_v = ((2. * plasma->B0 * M[m])/ species[s].m * F1(xc) + pow2(V[v]) * F3(x)) * df_dv + (6. * M[m] * pow2(V[v]) * F2(xc)) * df_dm - 2. * F3(xc)/pow3(V[v]) * V[v] * f;
  const double C_m = ((6. * M[m] * V[v] * F2(xc)) * df_dv  + 4./plasma->B0 * pow2(V[v]) * M[m] * F1(xc)) * df_dv + (4. * pow2(M[m]) * F3(xc)) * df_dm - 2. * F3(xc)/pow3(V[v]) * 2. * M[m] * f;
        
  // BUG : this is wrong, we need to differentiate over v
  const double Cvm = (C_v + C_m) / dv;

  return Cvm;
}

void Collision_PitchAngle::initData(hid_t fileID) 
{

}   
     
virtual void Collision_PitchAngle::printOn(std::ostream &output) const 
{
  output   << "Collisions |  Pitch-Angle " << ((nu > 0.) ? Num2String(nu) :  "off") << std::endl;
}

