cimport openmp

import numpy as np
cimport numpy as np
cimport cython

#from cython.parallel import prange
#@cython.boundscheck(False)
#@cython.wraparound(False)

cdef extern from "iCode-func.h":
    void setupMatrix(double complex w, double ky, double* X, double* K,  double Ls,  double Ln, int Nx, double complex* A,  double h, double q, double mass, double T, double eta, double rho) 


def setupMatrixPy(species, double complex w, double ky, np.ndarray[double, ndim=1] X, np.ndarray[double, ndim=1] K,  
                  double Ls,  double Ln, int Nx, np.ndarray[double complex, ndim=2] A, double dx, double lambda_D2):
   

    cdef int n
    #cdef double q,T,mass,rho,eta

    for n in range(Nx):
      
      # Add adiabatic response term
      # Add finite Debye length effect      
      A[n,n] = \
      - species[0].q**2 * species[0].n / species[0].T \
      + lambda_D2 * (K[n]**2 + ky**2) /2. 

    # Add kinetic species effects
    for sp in species[1:]:
       
        # Call C++ function to add up matrix
        setupMatrix(w, ky, <double*> X.data, <double*> K.data, Ls,  Ln, Nx, \
                   <double complex *> A.data, dx, sp.q, sp.m, sp.T, sp.eta, sp.n)
    
    return A


