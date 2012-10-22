from pylab import *
import scipy 
import mpmath as mp
import traceback
import numpy
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
Nx = 129

# My setup
#species = [ Species(m=0.,q=-1.,T=1.,n=1., eta=0.,name= "Adiab"),   Species(1.,1.,1.,1., 5., "Ion") ]#,   Species(1./1836.,-1.,1.,1., 4., "Electron")]
species = [  Species(name= "Adiab"),   Species(1.,1.,1.,1.,5., "Ion"),  Species(1./1837.,-1.,1.,1., 2., "Electron") ]
#Ln, Ls, Lx, Ly, lambda_D2, ky_list, w0 = 1., 1./0.2, 32., 64., 0., [0.4908] ,   -0.25 + 0.05j
# adiab Ln, Ls, Lx, Ly, lambda_D2, ky_list, w0 = 1., 1./0.2, 12., 64., 0., [2.55] ,   -0.25 + 0.05j
Ln, Ls, Lx, Ly, lambda_D2, ky_list, w0 = 1., 1./0.2, 12., 64., 0., [2.55] ,   -0.3 + 0.06j
#Ln, Ls, Lx, Ly, lambda_D2, ky_list, w0 = 1., 1./0.2, 32., 64., 0., logspace(-1, log10(2*pi), 32 ),   -0.25 + 0.05j

## Gao Setup
#species = [  Species(name= "Adiab"),   Species(m=1836.,q=1.,T=1.,n=1.,eta=0., name="Ion"),  Species(m=1.,q=-1.,T=1.,n=1., eta=3., name="Electron") ]
#Ln, Ls, Lx, Ly, lambda_D2, ky_list, w0 = 1., 0.025, 30., 256., 1., [0.3], -0.01 + 0.025j

######################## Setup Grid ######################

dx, kx_list, dk  = Lx/Nx, 2.*pi/Lx * linspace(-(Nx-1)/2., (Nx-1)/2., Nx), 2.*pi/Lx
X                = linspace(-Lx/2, Lx/2, Nx)

sol = []
        
def getMinEigenvalue(M):
            eigvals =  scipy.linalg.eigvals(M)
            idx = argmin(abs(eigvals))
            return eigvals[idx]
            """
            # get minimal non-negative eigenvalue
            eigvals_lz = eigvals[where(real(eigvals) > 0.)]
            idx = argmin(real(eigvals_lz))
            return eigvals_lz[idx]
            """

def solveDispersion(ky):

    L = zeros((Nx,Nx), dtype=complex)
    
    def setup_L(w):
        L[:,:] = 0.
        iCode.setupMatrixPy(species, w, ky, X, kx_list, Ls, Ln, Nx, L, dk*dx, lambda_D2)
        return L


    def solveEquation(w):

        w = complex(w)

        L = setup_L(w)
        
        # Get Minium absolute eigenvalue. 
        # 
        # The non-linear eigenvalue problem obeys 
        # det(L) = 0. This requirement is equal
        # to an eigenvalue 0 for the linear eigenvalue
        # problem det(L'- lambda I) = 0, which is identical
        # to the non-linear requirenement once lambda = 0.
        # 
        # [ This is less sensitive (and thus numerically less demanding) 
        # than caluculating directly the determinant
        
        
        val = getMinEigenvalue(L)

        #val = det(A)   
        #(sign, logdet) = np.linalg.slogdet(L)
        #val = sign * logdet

        print ky, " w   :  %.8f+%.8f j" % (real(complex(w)), imag(complex(w))) ,  "  Determinant : %.2e " %  abs(val)

        return val	
    
    try :
            #omega = complex(mp.findroot(solveEquation, (w0, w0+0.05, w0-0.005j), solver='muller', tol=1.e-15, ftol=1.e-15, maxsteps=5000))
            omega = complex(mp.findroot(solveEquation, (w0, w0-0.01, w0+0.001j), solver='muller', tol=1.e-15, ftol=1.e-15, maxsteps=5000))
    except:
            return float('nan') + 1.j * float('nan')
            #traceback.print_exc(file=sys.stdout)
    try:
        # n = 0
        
        # solution found for w0, get solution vector
        werr = solveEquation(omega)
        
        L = setup_L(omega)

        # We found our eigenvalue omega, now we use the
        # inverse iteration to find the closest eigenvector
        # to the eigenvalue
        L__lambda_I = L - omega * eye(Nx)

        # Start with random inital eigenvector
        b = 1. * rand(Nx) + 1.j * rand(Nx)

        # Convergence is fast thus large iteration number not required
        # However, how to best check the error ? 
        # e.g. A dot phi
        # Tested also, Rayleigh-Quotient Iteration but less successfull.
        for n in range(128):

            # rescale b
            b = b/real(sqrt(vdot(b,b)))
            # inverse iteration 
            b = solve(L__lambda_I, b)

            # calculate error
            r = dot(L,b)
            residual = real(sqrt(vdot(r,r)))
            print("I-I Residual : %.2e " % residual )
            if (residual < 1.e-9) : break
       
        clf()
        
        fig = gkcStyle.newFigure(ratio='1.41:1', basesize=9)
        
        
        fig.suptitle("$\omega_0$ = %.4f %.4fi  $\pm$ %.2e %.2e i" % (real(omega), imag(omega), real(werr), imag(werr)))
        
        ###################### Plot Fourier Modes ##########3
        plot(kx_list, real(b), 'r.-', label="real")
        plot(kx_list, imag(b),  '.-', label="imag", color=gkcStyle.color_indigo)
        
        xlim((min(kx_list), max(kx_list))) 

        xlabel("$k_x$")
        ylabel("$\phi(k_x)$")
        legend(ncol=2).draw_frame(0)
        
        savefig(str(ky) + "_Plot1.pdf", bbox_inches='tight')
        ################### Plot real modes ########3
        clf()

        # We have to transform to FFTW format, which has
        # form of [ k=0, k=1, ..., k = N/2, k = -(N/2-1), ..., k=-1 ] 
        F = append(append(b[Nx/2], b[Nx/2+1:]), b[:Nx/2])
        K = np.fft.ifft(F)
        K = append(append(K[Nx/2], K[Nx/2+1:]), K[:Nx/2])


        # Correct for phase
        K = K * exp(-1.j*arctan2(imag(sum(K)),real(sum(K))))
        #print " Phase : " , arctan2(imag(sum(K)),real(sum(K)))

        #F_list = append(append(kx_list[Nx/2], kx_list[Nx/2+1:]), kx_list[:Nx/2])
        #print "fft kx_list -----------> " , F_list
        
        plot(X, real(K), 'r.-', label='real')
        plot(X, imag(K),  '.-', label='imag', color=gkcStyle.color_indigo)
        plot(X,  abs(K),  'g-', label='abs', linewidth=5., alpha=0.5)
        
        xlim((min(X), max(X)))

        xlabel("$x$")
        ylabel("$\phi(x)$")
        legend(ncol=2).draw_frame(0)
        
        savefig(str(ky) + "_Plot2.pdf", bbox_inches='tight')
        ################ Plot Contour
        clf()
        y = linspace(0., Ly, 512)
        KxKy = zeros((Nx, 65), dtype=complex)
        nky = ky * Ly / (2.*pi)
        KxKy[:,nky] = K
        #np.fft.ifft(F)
        XY = np.fft.irfft(KxKy, axis=1, n=512)
        
        xlabel("$x$")
        ylabel("$y$")
        contourf(X, y, XY.T, 20, vmin=-abs(XY).max(), vmax=abs(XY).max())
        #contourf(X, y, XY.T, 20, vmin=-abs(XY).max(), vmax=abs(XY).max())
        colorbar()
        
        savefig(str(ky) + "_Plot3.pdf", bbox_inches='tight')

        # append and normalize
        sol.append(np.fft.ifft(b/abs(b).max()))
    except:
          traceback.print_exc(file=sys.stdout)
    return omega


w_list1 = []
def plotMode():
  for ky in ky_list: 
    wn = solveDispersion(ky)
    w_list1.append (wn)

def plotContours():

    ky = 0.5
    R = linspace(-0.5, 0.5, 8)
    I = linspace(-.3, 0.3, 8)
    V = zeros((len(R),len(I)), dtype=complex)

    fig = figure(figsize=(30,10))
   
    n = 0
    for r in range(len(R)):
     for i in range(len(I)):

      A = zeros((Nx,Nx), dtype=complex)
      iCode.setupMatrixPy(species, R[r]+1.j*I[i], ky, X, kx_list, Ls, Ln, Nx, A, dk*dx, lambda_D2)
      #val =  getMinEigenvalue(A)
      (sign, logdet) = np.linalg.slogdet(A)
      val = sign * logdet
      V[r,i] = val
      #print "x, y", R[r], I[i] , " r : ", val
      print n, "/", len(R) * len(I)
      n = n+1
  
    """
    #subplot(131)
    #norm = mpl.colors.Normalize(vmin = -1., vmax = 1.)
    #contourf(R,I,real(V), 100, vmin=-1., vmax=1., norm = norm)
    xlabel("Real")
    ylabel("Imag")
    cb = colorbar()
    cb.set_clim(vmin=-1, vmax=1)

    #subplot(132)
    #contourf(R,I,imag(V), 100, vmin=-1., vmax=1.)
    #norm = mpl.colors.Normalize(vmin = -1., vmax = 1.)
    contourf(R,I,imag(V), 100, vmin=-1., vmax=1., norm = norm)
    xlabel("Real")
    ylabel("Imag")
    cb = colorbar()
    cb.set_clim(vmin=-1, vmax=1)
    
    subplot(133)
    """
    pcolor(R,I,log10(abs(V)))
   
    xlabel("Real")
    ylabel("Imag")
    cb = colorbar()
    #cb.set_clim(vmin=0., vmax=1)
    
    
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

fig = figure(figsize=(30,10))
subplot(211)
semilogx(ky_list, real(w_list1), 'o-', label='real')
subplot(212)
semilogx(ky_list, imag(w_list1), 'o-', label='imag')
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



