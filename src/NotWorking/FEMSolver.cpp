/*
 * =====================================================================================
 *
 *       Filename:  FEMSolver.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/17/2010 08:15:59 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include "FEMSolver.h"

#include "string_to_enum.h"

FEMSolver::FEMSolver(Setup *_setup, Parallel *_parallel) : setup(_setup), parallel(_parallel)
{
return;
#ifdef PARALLEL_MPI
  LibMeshInit init (setup->argc, setup->argv, parallel->Comm[DIR_ALL]);
#else
  LibMeshInit init (setup->argc, setup->argv);//, parallel->Comm[DIR_ALL]);
#endif
  std::string order  = setup->get("LibMesh.Order" , "FIRST"); 
  std::string family = setup->get("LibMesh.Family", "LAGRANGE");

  check((LIBMESH_DIM < 2) ? -1 : 1, DMESG("LibMesh was conly compiiled with 1-d support, 2-D required as minimum"));
  
  
  
  Mesh mesh_xy;
 
  //MeshTools::Generation::build_square (mesh_xy, Nx, Ny, X(NxGlD), X(NxGuD), Y(NyGlD), Y(NyGuD), QUAD4);
  MeshTools::Generation::build_square (mesh_xy, Nx, Ny, 0., Lx, 0., Ly, QUAD4);
  // Print information about the mesh to the screen.
  mesh_xy.print_info();
  
  // Create an equation systems object.
  equation_systems  = new EquationSystems (mesh_xy);
 

// 1-D FEM

  equation_systems->add_system<LinearImplicitSystem> ("Poisson");
  equation_systems->get_system("Poisson").add_variable("u", Utility::string_to_enum<Order> (order), Utility::string_to_enum<FEFamily> (family));

  equation_systems->get_system("Poisson").attach_assemble_function(assemble_fields);
  equation_systems->init();
  equation_systems->print_info();
  mesh_xy.print_info();

   equation_systems->get_system("Poisson").solve();
}

  
FEMSolver::~FEMSolver() 
{




};

  


//
//
//
// We now define the matrix assembly function for the
// Poisson system.  We need to first compute element
// matrices and right-hand sides, and then take into
// account the boundary conditions, which will be handled
// via a penalty method.void FEMSolver::assemble_fields(EquationSystems& es, const std::string& system_name, N)
void FEMSolver::assemble_fields(EquationSystems& es, const std::string& system_name)
{
  libmesh_assert (system_name == "Poisson");

  //PerfLog perf_log ("Matrix Assembly");
  
  const MeshBase& mesh = es.get_mesh();

  const unsigned int dim = mesh.mesh_dimension();

  LinearImplicitSystem& system = es.get_system<LinearImplicitSystem>("Poisson");
  const DofMap& dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(0);

  // Build a Finite Element object of the specified type.  Since the
  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, CONSTANT);
  fe->attach_quadrature_rule (&qrule);
  AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));
  QGauss qface(dim-1, CONSTANT);
  fe_face->attach_quadrature_rule (&qface);
  const std::vector<Real>& JxW = fe->get_JxW();
  const std::vector<Point>& q_point = fe->get_xyz();
  const std::vector<std::vector<Real> >& phi = fe->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
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
      Ke.resize (dof_indices.size(), dof_indices.size());
      Fe.resize (dof_indices.size());

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        for (unsigned int i=0; i<phi.size(); i++)
          for (unsigned int j=0; j<phi.size(); j++)
            Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
      

      // Important part
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          // fxy is the forcing function for the Poisson equation.
          // In this case we set fxy to be a finite difference
          // Laplacian approximation to the (known) exact solution.
          //
          // We will use the second-order accurate FD Laplacian
          // approximation, which in 2D on a structured grid is
          //
          //
          // Since the value of the forcing function depends only
          // on the location of the quadrature point (q_point[qp])
          // we will compute it here, outside of the i-loop          
          const Real x = q_point[qp](0);
          const Real y = q_point[qp](1);
          const Real z = q_point[qp](2);
          const Real eps = 1.e-3;

          const Real uxx = 1.;
          const Real uyy = 1.;
          

          //Real fxy = n(x,y,z) - 1.;
          Real fxy = 1.;
         

          


          // Add the RHS contribution
          for (unsigned int i=0; i<phi.size(); i++)
            Fe(i) += JxW[qp]*fxy*phi[i][qp];          
        }
      
      // Stop logging the right-hand-side computation
      //perf_log.pop ("Fe");

      {
        
        for (unsigned int side=0; side<elem->n_sides(); side++)
          if (elem->neighbor(side) == NULL)
            {
              const Real penalty = 1.e10;

              const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
              const std::vector<Real>& JxW_face = fe_face->get_JxW();
              const std::vector<Point >& qface_point = fe_face->get_xyz();
              fe_face->reinit(elem, side);
              for (unsigned int qp=0; qp<qface.n_points(); qp++)
              {
                // The location on the boundary of the current
                // face quadrature point.
                const Real xf = qface_point[qp](0);
                const Real yf = qface_point[qp](1);


                // The boundary value.
                const Real value = 2.;

                // Matrix contribution of the L2 projection. 
                for (unsigned int i=0; i<phi_face.size(); i++)
                  for (unsigned int j=0; j<phi_face.size(); j++)
                    Ke(i,j) += JxW_face[qp]*penalty*phi_face[i][qp]*phi_face[j][qp];

                // Right-hand-side contribution of the L2
                // projection.
                for (unsigned int i=0; i<phi_face.size(); i++)
                  Fe(i) += JxW_face[qp]*penalty*value*phi_face[i][qp];
              } 
            }
      } 
      
      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      //perf_log.push ("matrix insertion");
      
      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);

      //perf_log.pop ("matrix insertion");
    }
}


