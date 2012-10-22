from numpy import *
    ########################### Eigenvalue calculations ##################
import scipy.linalg

def getMinAbsEigenvaluLA(K):
   eigvals =  scipy.linalg.eigvals(K)
   #return min(abs(eigvals))
   idx = argmin(abs(eigvals))
   return eigvals[idx]

def getMinAbsEigenvalue(K):
    
    import sys, slepc4py
    slepc4py.init(sys.argv)

    from petsc4py import PETSc
    from slepc4py import SLEPc

    N = len(K[0,:])

    opts = PETSc.Options()

    A = PETSc.Mat(); A.createDense((N,N))

    # copy matrices to petsc, is there are more direct way ?
    for i in range(N):
        for j in range(N):
            A[i,j] = K[i,j]


    A.assemble()

    E = SLEPc.EPS(); E.create()

    E.setOperators(A)
    # Non-hermite problem
    E.setProblemType(SLEPc.EPS.ProblemType.NHEP)
    E.setWhichEigenpairs(SLEPc.EPS.Which.SMALLEST_MAGNITUDE)
    E.setTolerances(tol=1.e-11, max_it=1000)
    E.setDimensions(1)
    
    E.solve()

    #nPrint = PETSc.Sys.Print

    eps_type = E.getType()
    #Print( "Solution method: %s" % eps_type )

    nconv = E.getConverged()
    Print( "Number of converged eigenpairs %d" % nconv )

    val = []
    #    # Create the results vectors
    vr, wr = A.getVecs()
    vi, wi = A.getVecs()
    #    #
    for i in range(nconv): val.append(E.getEigenpair(i, vr, vi))
    #        val.append(k)
    #        error = E.computeRelativeError(i)
    #    #E.getEigenvector(0, vr, vi)
    #   #return vr[:]+ vi[:]
    # print "Eigenvalues: ", val
    vals = array(val)
    idx = argmin(abs(vals))
    return vals[idx]

   #return array(val)[-1]
      
def getDeterminant(K):
    import sys, slepc4py
    slepc4py.init(sys.argv)

    from petsc4py import PETSc
    from slepc4py import SLEPc

    N = len(K[0,:])

    opts = PETSc.Options()

    A = PETSc.Mat(); A.createDense((N,N))
    #A.setFromOptions()

    # copy matrices to petsc, is there are more direct way ?
    for i in range(N):
        for j in range(N):
            A[i,j] = K[i,j]


    A.assemble()
    
    
    A.LU()

    print A
    return mul(A.getDiagonal())

      
def getDeterminant(K):
     f, logDet = linalg.slogdet(K)
     return f * exp(logDet)


