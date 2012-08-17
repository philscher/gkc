import gkcData
import gkcStyle

import pylab
import numpy as np
import scipy.optimize

import numpy.linalg as la




def getPOD(data, calc_type='eig'):
    if calc_type == 'svd':
        data = array(orig_data)
        data = (data - data.mean(axis=0)) / data.std(axis=0)
        u, s, v = linalg.svd(data)

        #normalize
        data /= (len(data[:,0])-1)
        print s #should be s**2 instead!
        print v
        eig_val = s
        eig_vec = v

    elif calc_type == 'eigen':
        # usually a good idea to mean center your data first:
        #data -= np.mean(data, axis=0)
        # calculate the covariance matrix 
        data = abs(data)*sign(arctan2(imag(data),real(data)))
        C = np.corrcoef(data)#, rowvar=0)
        # returns an m x m matrix, or here a 5 x 5 matrix)
        # now get the eigenvalues/eigenvectors of C:
        eig_val, eig_vec = la.eig(C)

    else : raise TypeError("No such calc_type ('svd' or 'eigen' allowed)")
       
    return eig_val, eig_vec

def plotCFL(fileh5):
  clf()
  semilogy(fileh.root.cfl.cols.time, fileh.root.cfl.cols.Total[:], "r")
  semilogy(fileh.root.cfl.cols.time, fileh.root.cfl.cols.Fx[:]   , "c")
  semilogy(fileh.root.cfl.cols.time, fileh.root.cfl.cols.Fy[:]   , "m")
  semilogy(fileh.root.cfl.cols.time, fileh.root.cfl.cols.Fz[:]   , "g")
  semilogy(fileh.root.cfl.cols.time, fileh.root.cfl.cols.Fv[:]   , "b")

  xlabel("Time")
  ylabel("Possible timestep $\\Delta t$")
  leg = legend(("Total", "$E_y \\frac{\\partial f}{\\partial x}$", "$E_x \\frac{\\partial f}{\\partial y}$", \
          "$E_z \\frac{\\partial f}{\\partial v_\parallel}$", "$v_\\parallel \\frac{\\partial f}{\\partial z}$"), ncol=2)
  leg.draw_frame(0)



