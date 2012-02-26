/*
 * =====================================================================================
 *
 *       Filename:  PoissonFEM.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/25/2010 11:42:31 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */


#ifndef __FIELDS_FEM_H
#define __FIELDS_FEM_H

#include "Global.h"

#ifdef HELIOS_FEM


#include "Fields.h"
#include "FEMSolver.h"
#include "Plasma.h" 

#include "quadrature_gauss.h"
#include "quadrature_grid.h"

//class PoissonSystem : public LinearImplicitSystem {
class PoissonSystem : public LinearImplicitSystem {
  Array3d n;
QBase *qrule;
AutoPtr<FEBase> fe;
    int z;
public:
    PoissonSystem(EquationSystems &es, const std::string &name, const unsigned int number);
    void solved(Array3d _n, Array3d phi0) {

//        phi.reference(_phi);
        n.reference(_n);
        
        
        for(z = NzLlD; z <= NzLuD; z++) {
          PoissonSystem::solve();
     
          
        const MeshBase &mesh = get_mesh();
        
        for(unsigned int i = 0; i < mesh.n_nodes(); i++) {
            const unsigned int dof_nr = mesh.node(i).dof_number(0,0,0);
            const Real x = mesh.node(i)(0);
            const Real y = mesh.node(i)(1);

            int x_idx = (int) (x/Lx * Nx) + 3;
            int y_idx = (int) (y/Ly * Ny) + 3;


            phi0(x_idx, y_idx, z) = current_solution(dof_nr);
        } 
        }
//        std::cout << phi0 << std::endl;
    };
    void assemble();
};





class PoissonFEM : public Poisson {
Mesh mesh;
EquationSystems *equation_systems;
PoissonSystem *fields_system;
  LibMeshInit *init;
virtual Array1z calcFluxSurfAvrg(Array3z rho_k);
virtual Array3d solvePoissonEquation(Timing timing);

FEMSolver *fem;


public:
 PoissonFEM(Setup *setup, Grid *grid, Parallel *parallel, Geometry *geo, FEMSolver *femsolver);
//* Destructor
~PoissonFEM();
 Array3d gyroAverage(Array3d A3, int m, int s) { return A3; };

protected:
        virtual void printOn(ostream &output) const {
         output   << "Poisson   |  FEM " << std::endl;
        }

};


#endif // __FIELDS_FEM_H
#endif // __FIELDS_FEM_H
