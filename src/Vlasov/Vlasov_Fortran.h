/*
 * =====================================================================================
 *
 *       Filename: Vlasov_Fortran.h
 *
 *    Description: Implementation of GK Vlasov's equation using 
 *                 Fortran.
 *
 *         Author: Paul P. Hilscher (2010), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef __VLASOV_FORTRAN_H
#define __VLASOV_FORTRAN_H

#include "Vlasov.h"

class VlasovFortran : public Vlasov {
        
        void Vlasov_EM_Linear   (Array6z fs, Array6z fss, Array6z ft, Fields *fields, double dt, int rk_step);
        void Vlasov_EM_NonLinear(Array6z fs, Array6z fss, Array6z ft, Fields *fields, double dt, int rk_step);
        void Vlasov_2D(Array6z fs, Array6z fss, Array6z ft, Fields *fields, double dt, int rk_step) {
            check(-1, DMESG("Not implemented"));

        };
  public:
        VlasovFortran(Grid *_grid, Parallel *_parallel, Setup *_setup, Geometry<HELIOS_GEOMETRY> *_geo, FFTSolver *fft); 
        int solve(std::string equation_tyoe, Fields *fields, Array6z fs, Array6z fss, double dt, int rk_step, int user_boundary_type=BOUNDARY_CLEAN) {
          check(-1, DMESG("FIX ME"));
        return HELIOS_SUCCESS; 
        }
  
  protected :

        void printOn(ostream &output) const;

};


#endif //  __VLASOV_FORTRAN_H

