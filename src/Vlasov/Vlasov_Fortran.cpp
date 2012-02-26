/*
 * =====================================================================================
 *
 *       Filename: Vlasov_Fortran.cpp
 *
 *    Description: Implementation of GK Vlasov's equation using 
 *                 Fortran.
 *
 *         Author: Paul P. Hilscher (2010), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */


#include "Vlasov.h"
#include "Vlasov/Vlasov_Fortran.h"
/* 
extern "C" {
void Fortran_Vlasov_EM_Linear(double *f, double *fs, double *fss,  double *ft, double *f0, double *G, double *Xi, double *phi, double *Ap, double &dt, int &rk_step, Species *species, double *Vel, double *Mu, double *B0, double *beta, double *eps, double *max_Xi, double *hyper_visc_z); 
void Fortran_Vlasov_EM_NonLinear(double *f, double *fs, double *fss,  double *ft, double *f0, double *G, double *Xi, double *phi, double *Ap, double &dt, int &rk_step, Species *species, double *Vel, double *Mu, double *B0, double *beta, double *eps, double *max_Xi, double *hyper_visc_z); 
void Morinishi_SA (double *v, double *z, double *f, double *fs, double *fss,  double *ft, double *f0, double *Ex, double *Ey, double *Ez, double &dt, int &rk_step, Species *species, double *Mu, double *B0, double *shear, double *R0, double *q0, double *epsilon_trap, int *calcTerms); 

};

 * */


void VlasovFortran::Vlasov_EM_Linear(Array6z _fs, Array6z _fss, Array6z _ft, Fields *fields, double dt, int rk_step) {
  double hyper_visc = 0;
/* 
     Fortran_Vlasov_EM_Linear(f0.data(), f.data(), _fs.data(), _fss.data(), ft.data(), G.data(), Xi.data(),  
                      fields->phi.data(), fields->Ap.data(), dt, rk_step, plasma->species.data(), V.data(), M.data(),
                      &plasma->B0, &plasma->beta, &geo->eps_hat, Xi_max.data(), &hyper_visc);
 * */
}

void VlasovFortran::Vlasov_EM_NonLinear(Array6z _fs, Array6z _fss, Array6z _ft, Fields *fields, double dt, int rk_step) {
/* 
  double hyper_visc = 0;
     Fortran_Vlasov_EM_NonLinear(f0.data(), f.data(), _fs.data(), _fss.data(), ft.data(), G.data(), Xi.data(),  
                      fields->phi.data(), fields->Ap.data(), dt, rk_step, plasma->species.data(), V.data(), M.data(),
                      &plasma->B0, &plasma->beta, &geo->eps_hat, Xi_max.data(), &hyper_visc);
 * */
}

VlasovFortran::VlasovFortran(Grid *_grid, Parallel *_parallel, Setup *_setup, Geometry<HELIOS_GEOMETRY> *_geo, FFTSolver *fft)    : Vlasov(_grid, _parallel, _setup, _geo, fft) 
{
}

void VlasovFortran::printOn(ostream &output) const
        {

            output << "Vlasov     |   Fortran wrapper with Mornishi scheme for Non-Lin"  << std::endl;
            Vlasov::printOn(output);

        };



