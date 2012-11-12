

def getEigenvalues(func, args, N=64, EigvType="All"):

        import sys, slepc4py
        slepc4py.init(sys.argv)

        from petsc4py import PETSc
        from slepc4py import SLEPc
        import numpy as np

    
        # Equation
        class A(object):
                def __init__(self): pass

                def mult(self, A, x, y):
                  y[...] = func(x[...], args)
            
        
        context = A()
        A = PETSc.Mat().createPython((N,N), context)
        A.setUp()
        """
            def setupBMatrix():
                B = PETSc.Mat(dtype=complex); B.create()
                B.setSizes([N,N])
                B.setFromOptions()
                for n in range(N):
                    for m in range(N):
                        if m == n : B[n,m] =  - eps *  2. * omega[n]
                        else      : B[n,m] = 0.
                B.assemble()
                return B

            B = setupBMatrix() 
        """

        
        S = SLEPc.EPS().create()
        
        S.setOperators(A)
        #S.setOperators(A,B)
        S.setFromOptions()
        A.setFromOptions()
        S.setProblemType(SLEPc.EPS.ProblemType.NHEP)
        S.setBalance()
        S.setDimensions(N)
      
        #F1 = PETSc.Vec().createSeq(Nv)
        #F1.setValues(range(Nv), getInitialF1(ky))
        #S.setInitialSpace(F1)

        S.setTolerances(tol=1.e-15, max_it=500000)
    
        S.solve()
        
        its = S.getIterationNumber()
        sol_type = S.getType()
        nev, ncv, mpd = S.getDimensions()
        tol, maxit = S.getTolerances()
        nconv = S.getConverged()
        
        print("Number of converged eigenpairs: %d" % nconv)
        
        eigval = []
        eigvec = []
        
        if nconv > 0:
                ### (1) Get All eigenvalues ###
                xr, tmp = A.getVecs()
                xi, tmp = A.getVecs()
                for i in range(nconv):
                    k = S.getEigenpair(i, xr, xi)
                    error = S.computeRelativeError(i)
                    #print(" %9f%+9f j  %12g" % (k.real, k.imag, error))
                    eigval.append(k.real + k.imag * 1.j)
                    eigvec.append(xr[...] + 1.j * xi[...])
      
        return np.array(eigval), np.array(eigvec)



