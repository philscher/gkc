import numpy as np

def cartesian(arrays, out=None):
    """
    Generate a cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

    """

    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out

def shootingZero2D(F, args, RangeX, RangeY):
    """
        Find minumum
    """
    w_err_min= 1.e99
    res_vec_min = 0.0 + 0.j

    grid = cartesian((RangeX, RangeY))

    for vec_w in grid:
      res_vec = F(vec_w, args[0], args[1])
      res     = res_vec[0] + 1.j * res_vec[1]

      if(abs(res) < w_err_min): 
        w_err_min = abs(res)
        res_vec_min = res_vec
    


    return (res, w_err_min)


"""
File          zermuller.py
Author      Ernesto P. Adorio,Ph.D.
email        ernesto.adorio@gmail.com
               UPDEPP at Clark Field
               Pampanga, the Philippines
license     educational use
 
Description Zero finding using muller's method.
Return      array of roots found
Reference   Conte and DeBoor, pp.120-122.
Revisions   04.03.2009 first Python release based on old 1998 C++ code by
                 the translator
"""
 
 
from cmath import *
def zermuller(Func, ky, S, xinit, ztol, ftol=8.0e-2, dx=0.1, maxiter=500, wantreal=False, nroots=1, verbose=False):
    nmaxiter = 0
    retflag  = 0
    roots = []


    for j in range(nroots):
        if(verbose) : print "j=",  j
        x1  = xinit
        x0  = x1 - dx
        x2  = x1 + dx
 
        f0, undeflate  = deflate(Func, x0, ky, S,  j, roots)
        f1, undeflate  = deflate(Func, x1, ky, S,  j, roots)
        f2, undeflate  = deflate(Func, x2, ky, S,  j, roots)

        h21 = x2 - x1
        h10 = x1 - x0
        f21 = (f2 - f1) / h21
        f10 = (f1 - f0) / h10
 
        for i in range(maxiter):
            if(x0 != x0) : return [float('nan')]
            if(verbose) : print "iter", i
            f210 = (f21 - f10) / (h21+h10) 
            b    = f21 + h21 * f210
            t    = b*b- 4.0 * f2 * f210
 
            if (wantreal) :     	# force real roots ? #
               if (real(t) < 0.0):
                   t = 0.0
               else :
                   t =  real(t)
            Q = sqrt(t)
            D = b + Q
            E = b - Q
 
            if (abs(D) < abs(E)) :
                D = E
 
 
            if (abs(D) <= ztol) :      # D is nearly zero ? #
                xm = 2 * x2 - x1
                hm = xm - x2
            else :
                hm = -2.0 * f2 / D
                xm = x2 + hm
 
 
            # compute deflated value of function at xm.  #
            fm, undeflate = deflate(Func, xm, ky, S, j, roots)
 
 
            # Divergence control #
            absfm = abs(fm)
            absf2 = 100. * abs(f2)
            # Note: Originally this was a while() block but it
            #       causes eternal cycling for some polynomials.
            #       Hence, adjustment is only done once in our version.
            if (absf2 > ztol and absfm >= absf2) :
                hm    = hm * 0.5
                xm    = x2 + hm
                fm    = Func(xm, ky, S)
                absfm = abs(fm)
 
 
            # function or root tolerance using original function
            if (abs(undeflate) <= ftol or abs(hm) <= ztol) :
                if (i > nmaxiter) :
                    nmaxiter = i
                    retflag = 0
                    break
 
            # Update the variables #
            x0  = x1
            x1  = x2
            x2  = xm
 
            f0  = f1
            f1  = f2
            f2  = fm
 
            h10 = h21
            h21 = hm
            f10 = f21
            f21 = (f2 - f1) / h21
 
 
        if (i > maxiter) :
                nmaxiter = i
                retflag  = 2
                break
 
        xinit = xm
        if(verbose) : print "a root is ", xinit
        roots.append(xinit)
 
        # initial estimate should be far enough from latest root.
        xinit    = xinit + 0.85
 
    maxiter = nmaxiter
    return roots
 
 
 
 
def deflate(f,z, ky, S, kroots, roots):
    """
    Arguments 
      f                 Input: complex<double> function whose root is desired
      z                 Input: test root
      kroots            Input: number of roots found so far
      roots             Input/Output: saved array of roots
 
    Return value
      Deflated value of f at z.
 
    Description
      This routine is local to zermuller.
      Basically, it divides the complex<double> function f by the product
           		(z - root[0]) ... (z - root[kroots - 1]).
      where root[0] is the first root, root[1] is the second root, ... etc.
    """
    undeflate = t = f(z,ky, S)
    nroots = len(roots)
    for i in range(nroots):
        denom = z - roots[i]
        while (abs(denom) < 1e-8):# avoid division by a small number #
            denom += 1.0e-8
        t = t / denom
    return t, undeflate
 






