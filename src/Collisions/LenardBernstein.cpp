/*
 * =====================================================================================
 *
 *       Filename: Collision_LenardBernstein.cpp
 *
 *    Description: Lenard-Bernstein collisional operator including
 *                 energy & momentum conservational terms
 *
 *         Author: Paul P. Hilscher (2012-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include "Collisions/LenardBernstein.h"
     
#include "Vlasov/Vlasov.h"

Collisions_LenardBernstein::Collisions_LenardBernstein(Grid *grid, Parallel *parallel, Setup *setup, FileIO *fileIO, Geometry *geo) 
: Collisions(grid, parallel, setup, fileIO, geo) 

{
  
  beta        = setup->get("Collisions.Beta", 0.e0);
  consvMoment = setup->get("Collisions.ConserveMoments", 0);
   
  // allocate arrays
  ArrayPreFactors     = nct::allocate(grid->RsLD, grid->RmLD , grid->RvLD)(&a , &b , &c , &nu);
  ArrayCorrectionTerm = nct::allocate(grid->RzLD, grid->RkyLD, grid->RxLD)(&dn, &dP, &dE);
      
  // as some terms include complicated functions we pre-calculate them 
  calculatePreTerms((A3rr) a, (A3rr) b, (A3rr) c, (A3rr) nu);

}

 
Collisions_LenardBernstein::~Collisions_LenardBernstein()
{


}
 

void Collisions_LenardBernstein::calculatePreTerms(double a[NsLD][NmLD][NvLD], double  b[NsLD][NmLD][NvLD], 
                                                   double c[NsLD][NmLD][NvLD], double nu[NsLD][NmLD][NvLD])
{

  for(int s = NsLlD; s <= NsLuD; s++) { 
  for(int m = NmLlD; m <= NmLuD; m++) { simd_for(int v = NvLlD; v <= NvLuD; v++) {

    // v = v_\parallel^2 + 2 \mu  / v_{\sigma, th}^2
    const double v_ = ( pow2(V[v]) + 2. * M[m] ) / pow2(species[s].v_th);

    nu[s][m][v] = v_ * dv * grid->dm[m]; // directly normalize 
    a [s][m][v] = 1. - 3. * sqrt(M_PI/2.) * ( erf(v_) - Derf(v_)) * pow(v_, -0.5);
    b [s][m][v] = V[v] * pow(v_, -3./2.) * erf(v_)    ;
    c [s][m][v] = pow(v_,-1./2.) *  (erf(v_) - Derf(v_)) ;
      
  } } } // s, m, v 

  return;
}



void Collisions_LenardBernstein::solve(Fields *fields, const CComplex  *f, const CComplex *f0, CComplex *Coll, double dt, int rk_step) 
{

  // Don't calculate collisions if collisionality is set to zero
  if (beta == 0.) return;

  [=](const CComplex f   [NsLD][NmLD][NzLB][Nky][NxLB][NvLB],  // Phase-space function for current timestep
      const CComplex f0  [NsLD][NmLD][NzLB][Nky][NxLB][NvLB],  // Background Maxwellian
            CComplex Coll[NsLD][NmLD][NzLB][Nky][NxLB][NvLB],  // Collisional term
            CComplex dn              [NzLD][Nky][NxLD]      ,
            CComplex dP              [NzLD][Nky][NxLD]      ,
            CComplex dE              [NzLD][Nky][NxLD]      ,
            const double a [NsLD][NmLD][NvLD],
            const double b [NsLD][NmLD][NvLD],
            const double c [NsLD][NmLD][NvLD],
            const double nu[NsLD][NmLD][NvLD]
     ) 
  {

    for(int s = NsLlD; s <= NsLuD; s++) {  if(consvMoment) {
    for(int m = NmLlD; m <= NmLuD; m++) { 

      // (1) Calculate Moments (note : can we recycle Fields  n, P ?)
      const double pre_dvdm = species[s].n0 * plasma->B0 * dv * grid->dm[m] ;

      #pragma omp for collapse(2)
      for(int z = NzLlD; z <= NzLuD; z++) { for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { 
      for(int x = NxLlD; x <= NxLuD; x++) {

        dn[z][y_k][x] = __sec_reduce_add(f[s][m][z][y_k][x][NvLlD:NvLD]                       ) * pre_dvdm;
        dP[z][y_k][x] = __sec_reduce_add(f[s][m][z][y_k][x][NvLlD:NvLD] * V       [NvLlD:NvLD]) * pre_dvdm;
        dE[z][y_k][x] = __sec_reduce_add(f[s][m][z][y_k][x][NvLlD:NvLD] * nu[s][m][NvLlD:NvLD]);

        
      } } } }
          
      // need to communicate with other CPU's ? 
    }
          
    
    // b-cast dn, dP, dE 
    const double _kw_12_dv_dv = 1./(12. * pow2(dv)); // use from e.g. Grid ?
    const double _kw_12_dv    = 1./(12. * dv )     ;

    
    #pragma omp for collapse(3)
    for(int   m = NmLlD ; m   <= NmLuD ; m++  ) { for(int z = NzLlD; z <= NzLuD; z++) {  
    for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { for(int x = NxLlD; x <= NxLuD; x++) {
           
    simd_for(int v = NvLlD; v <= NvLuD; v++) {

      // Velocity derivaties for Lennard-Bernstein Collisional Model (shouln't I use G ?)
      const CComplex f_      = f[s][m][z][y_k][x][v];
      const CComplex df_dv   = (8. *(f[s][m][z][y_k][x][v+1] - f[s][m][z][y_k][x][v-1]) 
                                  - (f[s][m][z][y_k][x][v+2] - f[s][m][z][y_k][x][v-2])) * _kw_12_dv;
      const CComplex ddf_dvv = (16.*(f[s][m][z][y_k][x][v+1] + f[s][m][z][y_k][x][v-1]) 
                                  - (f[s][m][z][y_k][x][v+2] + f[s][m][z][y_k][x][v-2]) 
                              -  30.*f[s][m][z][y_k][x][v]) * _kw_12_dv_dv;
      const double v2_rms    = 1.;//pow2(alpha)
        
      Coll[s][m][z][y_k][x][v] =
                 
        beta  * (f_  + V[v] * df_dv + v2_rms * ddf_dvv)           ///< Lennard-Bernstein Collision term
        // add conservation terms 
        + (consvMoment ?
          (a[s][m][v] * dn[z][y_k][x] +                           ///< Density  correction   
           b[s][m][v] * dP[z][y_k][x] +                           ///< Momentum correction
           c[s][m][v] * dE[z][y_k][x]) * f0[s][m][z][y_k][x][v]   ///< Energy   correction
                       : 0.);

    } // v

    } } } } // x, y_k, z, m 
         
   
  } } ((A6zz) f  , (A6zz) f0, (A6zz) Coll, 
       (A3zz) dn , (A3zz) dP, (A3zz) dE, 
       (A3rr) a  , (A3rr) b , (A3rr)  c, (A3rr) nu);
   
  return;

};


 void Collisions_LenardBernstein::printOn(std::ostream &output) const 
{

  output   << "Collisions |  Drift-Kinetic Lenard-Bernstein  β = " << beta 
           << " Corrections : " << (consvMoment  ? "Yes" : "No") <<  std::endl;
  
}

void Collisions_LenardBernstein::initData(hid_t fileID) 
{
  hid_t collisionGroup = check(H5Gcreate(fileID, "/Collisions",H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), DMESG("Error creating group file for Collision : H5Gcreate"));
     
  check(H5LTset_attribute_string(collisionGroup, ".", "Model", "Lenard-Bernstein"), DMESG("H5LTset_attribute"));
  check(H5LTset_attribute_double(collisionGroup, ".", "Beta" , &beta, 1), DMESG("H5LTset_attribute"));
            
  H5Gclose(collisionGroup);

}

