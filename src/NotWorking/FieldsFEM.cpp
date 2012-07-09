/*
 * =====================================================================================
 *
 *       Filename:  PoissonFEM.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/25/2010 11:39:54 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include "PoissonFEM.h"



 // use it to instrument the matrix assembly code and look
 // for bottlenecks where we should focus optimization efforts.
 //
 // This example also shows how to extend example 3 to run in
 // parallel.  Notice how litte has changed!  The significant
 // differences are marked with "PARALLEL CHANGE".

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>
#include <iostream>
#include <fstream>
// Basic include file needed for the mesh functionality.
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "exodusII_io.h"
#include "linear_implicit_system.h"
#include "equation_systems.h"

// Define the Finite Element object.
#include "fe.h"

// Define Gauss quadrature rules.
#include "periodic_boundaries.h"
// Define the DofMap, which handles degree of freedom
// indexing.
#include "dof_map.h"

// Define useful datatypes for finite element
// matrix and vector components.
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "quadrature_monomial.h"


// The definition of a geometric element
#include "elem.h"

#include "getpot.h"

#include "blitz/array.h"

using namespace libMesh;
int signI(double T) { 
  if (T > 0.) return 1;
  if (T < 0.) return -1;
  return 0;
}
  

PoissonSystem::PoissonSystem(EquationSystems &es, const std::string &name, const unsigned int number) : LinearImplicitSystem(es, name, number) {

    add_variable("phi", FIRST, LAGRANGE);
//    std::cout << get_linear_solver().solver_type();
  //  std::cout << get_linear_solver().get_info(); 

//  Utility::string_to_enum<Order> get_linear_solver().preconditioner_type();
     
}




void PoissonSystem::assemble()
{
  const MeshBase& mesh = get_equation_systems().get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

 // AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));
    DofMap& dof_map = get_dof_map();
    FEType fe_type = dof_map.variable_type(0);
    //AutoPtr<FEBase> fe
    //(FEBase::build(dim, fe_type));
    fe  = (FEBase::build(2, fe_type));
    qrule = new QGrid(2, FIRST);
//    qrule = new QMonomial(2, FIRST);
    fe->attach_quadrature_rule (qrule);

  const std::vector<Real>& JxW = fe->get_JxW();
  const std::vector<Point>& q_point = fe->get_xyz();
  const std::vector<std::vector<Real> >& phi = fe->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
    //DofMap& dof_map = get_dof_map();
  

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  std::vector<unsigned int> dof_indices;

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
            


  for ( ; el != end_el; ++el)
    {

      const Elem* elem = *el;

      dof_map.dof_indices (elem, dof_indices);
      fe->reinit (elem);
      Ke.resize (dof_indices.size(),  dof_indices.size());
      Fe.resize (dof_indices.size());
      
      for (unsigned int qp=0; qp<qrule->n_points(); qp++) {
            
           const Real x = q_point[qp](0);
           const Real y = q_point[qp](1);
            //std::cout << "XYZ : " << xyz[qp] << std::endl;
            // add GhostCells + Offset
            int x_idx = (int) (x/Lx * Nx) + 3;
            int y_idx = (int) (y/Ly * Ny) + 3;

            double f1 = (x - X(x_idx))/dx;
            double f2 = (y - Y(y_idx))/dy;
            int x_off = 0, y_off = 0;
            
            x_off = signI(f1);
            y_off = signI(f2);
            double fac = max(x_off , y_off)/2.;
//            std::cout << fac;
            //std::cout << f1 << "  " << f2 << std::endl; 
  //          double value = ((1.-fac) * n(x_idx, y_idx, z) + fac * n(x_idx+x_off, y_idx+y_off, z));
      double        value = -n(x_idx, y_idx, z);
            //std::cout << x << "/" << X(x_idx) << " " <<  y << "/" << Y(y_idx) << std::endl;
//        const double norm   = plasma->species(1).n0 * pow2(plasma->species(1).q)/plasma->species(1).T0;
        const double norm   = plasma->species(1).n0 * pow2(plasma->species(1).q)/plasma->species(1).T0;
        const double rho_t2 = plasma->species(1).T0 * plasma->species(1).m / (pow2(plasma->species(1).q) * plasma->B0);
        const double adiab  = plasma->species(0).n0 * pow2(plasma->species(0).q)/plasma->species(0).T0;
        

        // Matrix assembly
        for (unsigned int i=0; i<phi.size(); i++) {
          
             Fe(i) += JxW[qp]*value*phi[i][qp];
         
          // Set Stiffness matrix 
            for (unsigned int j=0; j<phi.size(); j++) Ke(i,j) += JxW[qp]*((plasma->debye2 + rho_t2 * norm) * dphi[i][qp]*dphi[j][qp] + adiab * phi[j][qp]*phi[i][qp]);
          
        } 
      
      
      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
      
      matrix->add_matrix (Ke, dof_indices);
      rhs->add_vector    (Fe, dof_indices);
      }
    }

};


PoissonFEM::PoissonFEM(Setup *setup, Grid *grid, Parallel *parallel, Geometry *geo, FEMSolver *femsolver) : Poisson(setup, grid, parallel,geo), fem(femsolver)
{
  init = new LibMeshInit(setup->argc, setup->argv, parallel->Comm[DIR_ALL]);
  // LibMeshInit blib(setup->argc, setup->argv, parallel->Comm[DIR_ALL]);
  //GetPot command_line (setup->argc, setup->argv);
  
  MeshTools::Generation::build_square(mesh, Nx-1,Ny-1, X(NxGlD), X(NxGuD), Y(NyGlD), Y(NyGuD),QUAD4);
  equation_systems =  new EquationSystems(mesh);
  fields_system   = &equation_systems->add_system<PoissonSystem>("Poisson");

  // set here the periodic boundary condition
  if(setup->get("Helios.Poisson.Boundary.X", "Periodic") == "Periodic") {
    PeriodicBoundary pbc1;
    pbc1.myboundary = 0;
    pbc1.pairedboundary = 2;
    //RealVectorValue boundary_translation1(0., Ly, 0.);
    RealVectorValue boundary_translation1(0., Y(NyGuD)-Y(NyGlD), 0.);
    pbc1.translation_vector = boundary_translation1;
    fields_system->get_dof_map().add_periodic_boundary(pbc1);
  } else check(-1, DMESG("Othher Boundaries not implemented yet"));

  if(setup->get("Helios.Poisson.Boundary.Y", "Periodic") == "Periodic") {
     // set here the periodic boundary condition
    PeriodicBoundary pbc2;
    pbc2.myboundary = 1;
    pbc2.pairedboundary = 3;
    //RealVectorValue boundary_translation2(-Lx, 0., 0.);
    RealVectorValue boundary_translation2(-(X(NxGuD)-X(NxGlD)), 0., 0.);
    pbc2.translation_vector = boundary_translation2;
    fields_system->get_dof_map().add_periodic_boundary(pbc2);
  } else check(-1, DMESG("Othher Boundaries not implemented yet"));


   
  equation_systems->init();



}

PoissonFEM::~PoissonFEM()
{

delete fields_system;
delete equation_systems;
delete init;

}



Array1z PoissonFEM::calcFluxSurfAvrg(Array3z rho_k)
{




}


Array3d PoissonFEM::solvePoissonEquation(Timing timing)
{
  //* Electric potential and Electric fields
   // for(int z=NzLlD; z<= NzLuD;z++){ for(int y=NyLlD; y<= NyLuD;y++){ for(int x=NxLlD; x<= NxLuD;x++) { n(x,y,z) = 1. + 10.e-10 * sin(X(x)*2*M_PI/20.) *sin(Y(y)*2*M_PI/30.)  ; }}}
//  n = 1.e10 * n;
  fields_system->solved(n, phi0);
//  phi0 = 1.e-10 * phi0;
  return phi0;

}






