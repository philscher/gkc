/*
 * =====================================================================================
 *
 *       Filename: Collision_HyperDiffusion.cpp
 *
 *    Description: Hyper-Diffusion collisional operator to test collapse
 *                 of CvK modes
 *
 *         Author: Paul P. Hilscher (2013-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include "Collisions/HyperDiffusion.h"
#include "Vlasov/Vlasov.h"

Collisions_HyperDiffusion::Collisions_HyperDiffusion(Grid *grid, Parallel *parallel, Setup *setup, FileIO *fileIO, Geometry *geo) 
: Collisions(grid, parallel, setup, fileIO, geo) 

{
  // Get beta_C for each species
  for(int s = 0; s <= NsGuD; s++) beta[s] = setup->get("Collisions.Species" + Setup::num2str(s) + ".Beta", 0.e0);
  
  initData(setup, fileIO);
}

 
Collisions_HyperDiffusion::~Collisions_HyperDiffusion()
{


}
 

void Collisions_HyperDiffusion::solve(Fields *fields, const CComplex  *f, const CComplex *f0, CComplex *Coll, double dt, int rk_step) 
{
 

  // Don't calculate collisions if collisionality is set to zero
  if (__sec_reduce_add(std::abs(beta[NsGlD:Ns])) == 0.) return;
  
  const double _kw_dv4 = 1./pow4(dv);

  [=](const CComplex f   [NsLD][NmLD][NzLB][Nky][NxLB][NvLB],  // Phase-space function for current timestep
      const CComplex f0  [NsLD][NmLD][NzLB][Nky][NxLB][NvLB],  // Background Maxwellian
            CComplex Coll[NsLD][NmLD][NzLB][Nky][NxLB][NvLB]   // Collisional term
     ) 
  {
    for(int s = NsLlD; s <= NsLuD; s++) {
    
    #pragma omp for collapse(3) nowait
    for(int   m = NmLlD ; m   <= NmLuD ; m++  ) { for(int z = NzLlD; z <= NzLuD; z++) {  
    for(int y_k = NkyLlD; y_k <= NkyLuD; y_k++) { for(int x = NxLlD; x <= NxLuD; x++) {


    simd_for(int v = NvLlD; v <= NvLuD; v++) {

      //const int    sc_idx[]  = { v-2, v-1,    v, v+1, v+2};  // stencil index
      //const double sc_val[]  = { -1., 16., -30., 16., -1.};  // stencil values
      
      const int    sc_idx[]  = { v-2, v-1,    v, v+1, v+2};  // stencil index
      const double sc_val[]  = { 1., -4., 6., -4., 1.};  // stencil values

      Coll[s][m][z][y_k][x][v] = -beta[s] * _kw_dv4 * __sec_reduce_add(sc_val[:] * f[s][m][z][y_k][x][sc_idx[:]]);
                 
    } // v

    } } } } // x, y_k, z, m 
    
    } // s
   
  } ((A6zz) f, (A6zz) f0, (A6zz) Coll); 
}


void Collisions_HyperDiffusion::printOn(std::ostream &output) const 
{

  auto arr2str = [=](const double *val, const int len) -> std::string {

    int prec = 2;
    std::ostringstream ss;
    ss << std::setprecision(prec) << std::scientific;

    // only add " / "  between two numbers
    for(int n = 0; n < len; n++) ss << val[n] << ( n == len-1 ? "" : " / ");

    return ss.str();
  };

  output   << "Collisions |  Hyper-Diffusion (Order : " << hyp_order <<")  β = " << arr2str(&beta[1], Ns) << std::endl;
}

void Collisions_HyperDiffusion::initData(Setup *setup, FileIO *fileIO)
{
  hid_t collisionGroup = fileIO->newGroup("Collisions");
     
  check(H5LTset_attribute_string(collisionGroup, ".", "Model", "Hyper-Diffusion"), DMESG("H5LTset_attribute"));
  check(H5LTset_attribute_double(collisionGroup, ".", "Beta" ,  beta, Ns+1), DMESG("H5LTset_attribute"));
            
  H5Gclose(collisionGroup);
}

