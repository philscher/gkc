
#ifndef __FEMSolver_H
#define __FEMSolver_H


#include "config.h"

#ifdef HELIOS_FEM

#include "Global.h"
#include "Parallel.h"
#include "Setup.h"

// Function prototype.  This is the function that will assemble
// the linear system for our Poisson problem.  Note that the
// function will take the  EquationSystems object and the
// name of the system we are assembling as input.  From the

// C++ include files that we need
#include <algorithm>
#include <math.h>

// Basic include file needed for the mesh functionality.
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "linear_implicit_system.h"
#include "equation_systems.h"

#include "fe.h"
#include "quadrature_gauss.h"




// Define the DofMap, which handles degree of freedom
// indexing.

// Define useful datatypes for finite element
// matrix and vector components.
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "dense_matrix.h"
#include "dense_vector.h"

#include "elem.h"
#include "dof_map.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;


class FEMSolver {
 
  EquationSystems *equation_systems;
protected:
  Parallel *parallel;
  Setup   *setup;
  
  static void assemble_fields(EquationSystems& es,  const std::string& system_name);

public:
  
  FEMSolver(Setup *_setup, Parallel *_parallel);
  ~FEMSolver(); 


};

#else // HELIOS_FEM

class FEMSolver {
 
public:
  
  FEMSolver(Setup *_setup, Parallel *_parallel) {};
  ~FEMSolver() {}; 

};

#endif // HELIOS_FEM


#endif // __FEMSOLVER_H
