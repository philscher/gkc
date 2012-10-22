from pylab import *
import scipy 
#import math
import mpmath as mp
import traceback
import random
import numpy
import Dispersion_ConstTheta 
#import fpectl
#fpectl.turnon_sigfpe()
import scipy.linalg as la
import scipy.sparse.linalg as sla
import SlepcDet
import gkcStyle

import iCode


class Species:
  def __init__(self, m=0., q=1., T=1., n=0., eta=0., name="Unamed"):
    self.m    = m
    self.q    = q
    self.T    = T
    self.n    = n
    self.eta  = eta
    self.name = name  


############################## Settings for Integral Mode ######################################
Nx = 65

# Gao case Ln, Ls, Lx, Ly, theta, lambda_D2 = 1., 40., 12., 32., 0.1, 1.




# My setup
species = [  Species(m=0.,q=-1.,T=1.,n=1., eta=0.,name= "Adiab"),   Species(1.,1.,1.,1., 5., "Ion") ]#,   Species(1./1836.,-1.,1.,1., 4., "Electron")]
#species = [  Species(name= "Adiab"),   Species(m=1.,q=1.,T=1.,n=1.,eta=5., name="Ion"),  Species(m=0.0025,q=-1.,T=1.,n=1., eta=0., name="Electron") ]
Ln, Ls, Lx, Ly, lambda_D2, ky_list = 1., 1./0.2, 64., 64., 0., [0.5]


## Gao Setup
species = [  Species(name= "Adiab"),   Species(m=1836.,q=1.,T=1.,n=1.,eta=0., name="Ion"),  Species(m=1.,q=-1.,T=1.,n=1., eta=3., name="Electron") ]
Ln, Ls, Lx, Ly, lambda_D2, ky_list = 1., 0.025, 60., 64., 0., 2.*pi/64. * arange(1, 8)
#Ln, Ls, Lx, Ly, lambda_D2, ky_list = 1., 0.025, 60., 64., 0., [0.3]

######################## Setup Grid ######################

kx_list  = 2*pi/Lx * linspace(-Nx/2., Nx/2., Nx)
X        = linspace(-Lx/2, Lx/2, Nx)

dx, dk  = Lx/Nx, 2.*pi/Lx

dx, dk = dx * dk *dk  , 1.

fig = figure(figsize=(30,10))



global w_min, D_min
w_min = 0.+0.j
D_min = 1e99 + 1.j*1.e99 

sol = []


def solveDispersion(ky):
    A = zeros((Nx,Nx), dtype=complex)
    
    def setupA(w):
        A[:,:] = 0.

        iCode.setupMatrixPy(species, w, ky, X, kx_list, Ls, Ln, Nx, A, dk*dx, lambda_D2)
        return A


    def solveEquation(w):
        global D_min, w_min

        A = setupA(complex(w))
        #print A 
        #val = SlepcDet.getMinAbsEigenvalue(A)
        val = SlepcDet.getMinAbsEigenvaluLA(A)

        #val = det(A)   
        #(sign, logdet) = np.linalg.slogdet(A)
        #val = sign * logdet

        if abs(val) < abs(D_min) : w_min = complex(w) 

        print ky, " w   :  %.3f+%.3f j" % (real(complex(w)), imag(complex(w))) ,  "  Determinant : %.2e " %  abs(val)

        if val != val: return 0. + 0.j
        return val	
    
    try :
            w0=  -0.01 + 0.02j
            w_damp = complex(mp.findroot(solveEquation, (w0, w0-0.005j, w0+0.005), solver='muller', tol=1.e-8, ftol=1.e-8, maxsteps=5000))
            #w_damp = PyPPL.getZero(solveEquation, init=(w0, w0+0.01j, w0+0.02), solver='muller', tol=1.e-9, ftol=1.e-6, maxsteps=5000)
    except:
          traceback.print_exc(file=sys.stdout)
    try:
        #for n in range(Nx):
        n = 0
        global w_min
        print "-----------------> ", w_min
        
        # solution found for w0, get solution vector
        werr = solveEquation(w_min)
        A = setupA(w_min)

        #print A
        
    
        #S = solve(A, append(1.+0.j,zeros(Nx-1, dtype=complex)))
        #S = solve(A, append(1.+0.j, append(zeros(Nx-2, dtype=complex), 1.+0.j)))
        #b = append(0., append(1.+0., zeros(Nx-2, dtype=complex)))
        #b = zeros(Nx, dtype=complex)
        #b = ones(Nx, dtype=complex)
        #b[:] = 0. ; 
        #b[0] = 1.
        #S, err = solve(A, b), 0.
        #S, err = sla.lgmres(A,b, tol=1.e-9)

        # We found our eigenvalue w_min, now we use the
        # inverse iteration to find the closest eigenvector
  
        I = (1.+0.j) * eye(Nx)
        b = (1.+1.j) * ones(Nx, dtype=complex)

        for n in range(4096):
            b_prev = b
            b = solve(A - w_min * I, b)
            # RESCALE
            b  = b / sum(abs(b))
            if (abs(sum( sqrt(sum(b**2)/sum(b_prev**2))   * b_prev - b    )) < 1.e-10) : break
            #print("Eigv Error : %.2e Abs : %.2e " % (abs(sum( sqrt(sum(b**2)/sum(b_prev**2))   * b_prev - b    )), sum(abs(b))) )
         
        #print "Sol : " , b

        clf()
        
        gkcStyle.newFigure(ratio='2.33:1', basesize=12)
        
        subplot(131)
        fig.suptitle("$\omega_0$ = %.4f %.4fi  $\pm$ %.2e %.2e i" % (real(w_min), imag(w_min), real(werr), imag(werr)))
        
        ###################### Plot Fourier Modes ##########3
        b = -real(b) + 1.j * imag(b)
        plot(kx_list, real(b), 'r.-', label="real")
        plot(kx_list, imag(b),  '.-', label="imag", color=gkcStyle.color_indigo)
        xlim((min(kx_list), max(kx_list))) 

        xlabel("$k_x$")
        ylabel("$\phi(k_x)$")
        legend(ncol=2).draw_frame(0)
        
        ################### Plot real modes ########3
        subplot(132)
        
        # Remove High frequency modes
        #b[:3] = 0.;
        #b[-4:] = 0.;
        # We have to transform to FFTW format 
        F = append(append(b[Nx/2], b[Nx/2+1:]), b[:Nx/2])
        print "F--------------->", F 
        plot(X,real(np.fft.ifft(F)), 'r.-', label='real')
        plot(X,imag(np.fft.ifft(F)),  '.-', label='imag', color=gkcStyle.color_indigo)
        
        xlim((min(X), max(X)))
        xlabel("$x$")
        ylabel("$\phi(x)$")
        legend(ncol=2).draw_frame(0)
        
        ################ Plot Contour
        subplot(133)
        y = linspace(0., Ly, 128)
        KxKy = zeros((Nx, 65), dtype=complex)
        nky = ky * Ly / (2.*pi)
        KxKy[:,nky] = np.fft.ifft(F)
        XY = np.fft.irfft(KxKy, axis=1)
        
        xlabel("$x$")
        ylabel("$y$")
        contourf(X, y, XY.T, 20, vmin=-abs(XY).max(), vmax=abs(XY).max())
        colorbar()
        
        savefig("Plot2_" + str(ky) + ".png", bbox_inches='tight')

        # append and normalize
        sol.append(np.fft.ifft(b/abs(b).max()))
    except:
#species = [  Species(name= "Adiab"),   Species(m=1.,q=1.,T=1.,n=1.,eta=5., name="Ion"),  Species(m=0.0025,q=-1.,T=1.,n=1., eta=0., name="Electron") ]
          traceback.print_exc(file=sys.stdout)
    return w_min, abs(solveEquation(w_min))


w_list1 = []
def plotMode():
  for ky in ky_list: 
    wn, err = solveDispersion(ky)
    w_list1.append (wn)

def plotContours():

    ky = 0.5
    R = linspace(-1.5, 0.5, 16)
    I = linspace(0., 10., 16)
    V = zeros((len(R),len(I)), dtype=complex)
    for r in range(len(R)):
     for i in range(len(I)):

      A = zeros((Nx,Nx), dtype=complex)
      iCode.setupMatrixPy(species, R[r]+1.j*I[i], ky, X, kx_list, Ls, Ln, Nx, A, dk*dx, lambda_D2)
      val = det(A)
      #(sign, logdet) = np.linalg.slogdet(A)
      #val = sign * logdet
      V[r,i] = val
      print "x, y", R[r], I[i] , " r : ", val
   
    subplot(131)
    contourf(R,I,real(V), 100)
    colorbar()
    subplot(132)
    contourf(R,I,imag(V), 100)
    colorbar()
    subplot(133)
    contourf(R,I,abs(V), 100)
    colorbar()
    #pcolor(R,I,imag(V))
    savefig("Contour.png")
    
    #print "(Integral)  Solution  is w : ",w0
    #print "(local)     Solution  is w : ",w_Local

#plotContours()
plotMode()

################################## Plot Figures ############################


### Plot 
clf()
ky_list = array(ky_list)

plot(ky_list, real(w_list1), 'o-', label='real')
plot(ky_list, imag(w_list1), 'o-', label='imag')
legend(ncol=2, loc='best').draw_frame(0)
xlim((min(ky_list), max(ky_list)))

savefig("Results.png")

"""

# Make 3D Plot kx, ky, z
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt

clf()
Z = array(sol)
ax = fig.add_subplot(121, projection='3d')

_X,_ky = np.meshgrid(X,ky_list)
ax.plot_surface(_X, _ky, real(Z), rstride=1, cstride=1, cmap=cm.jet)
ax = fig.add_subplot(122, projection='3d')
ax.plot_surface(_X, _ky, imag(Z), rstride=1, cstride=1, cmap=cm.jet)
#ax.set_zlim3d(0, 1)
ax.set_xlabel(r'$\phi_\mathrm{real}$')
ax.set_ylabel(r'$\phi_\mathrm{im}$')
ax.w_yaxis.set_scale("log")
savefig("Results_3D.png")
"""



