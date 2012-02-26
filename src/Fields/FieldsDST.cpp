/*
 * =====================================================================================
 *
 *       Filename:  FieldsDST.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  01/18/2012 05:59:18 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */



#include "FieldsDST.h"


FieldsDST::FieldsDST(Setup *setup, Grid *grid, Parallel *parallel, FileIO * fileIO, Geometry<HELIOS_GEOMETRY> *geo, FFTSolver *fftsolver) : 
    FieldsFFT(setup, grid, parallel, fileIO, geo, fftsolver)
{

        // check if we support DST



};


FieldsDST::~FieldsDST()
{

};
  

Array3z FieldsDST::solvePoissonEquation(Array3z rho, Timing timing) {

        FieldsFFT::solvePoissonEquation(rho, timing);
        
        for(int y_k=NkyLlD; y_k<= NkyLuD;y_k++) { for(int z=NzLlD; z<=NzLuD;z++) { for(int x_k=fft->K1xLlD; x_k<= fft->K1xLuD;x_k++) {

          fft->kXIn(x_k,y_k,z, Field::phi) = cmplxd(0., imag(fft->kXIn(x_k,y_k,z, Field::phi)));

        }}}

        return rho;


}
