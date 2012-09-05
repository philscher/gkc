/*
 * =====================================================================================
 *
 *       Filename:  Transport.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  03/27/2012 01:07:17 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */


#include "Transport.h"


Transport::Transport(Setup *setup, Grid *grid, Parallel *parallel, FileIO *fileIO, Analysis *analysis)
{





}


void Transport::update(Vlasov *vlasov) {
/*
 *
    // calculate new temperature distribution and update maxwellian
    // moments, mea, 
    //
    // mean temperature, skeenes, and max

   for(int s = NsLlD; s <= NsLuD; s++) {
   for(int z=NzLlD; z<= NzLuD;z++) { 
   #pragma omp parallel for
   for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { 
     for(int x=NxLlD; x<= NxLuD;x++) { 
                   
    // add old Maxwellian
    vlasov-f(x,y_z,v,m,s) -= vlasof->f0(x,y,z,v,m,s);

      for(int m = NmLlD; m <= NmLuD; m++) { for(int v = NvLlD; v <= NvLuD; v++) {
            
        T = (vlasov->f0(x,y_k,z,v,m,s) + vlasov->f(x,y_k,z,v,m,s)) * (pow2(V(v)) + M(m));

    }

      for(int m = NmLlD; m <= NmLuD; m++) { for(int v = NvLlD; v <= NvLuD; v++) {
        vlasov->f0(x, y_k, z, v, m, s)  =  (n / pow( M_PI*T, 1.5) * exp(-pow2(V(v))) * 
                    ((plasma->species(s).doGyro == true) ?   exp(- M(m)    * plasma->B0/T) :  T/(plasma->B0)));
    } }

    // get new Maxwellian
    vlasov->f(x,y,z,v,m,s) -= vlasov->f0(x,y_k,z,v,m,s);

    }}}
*/
};

