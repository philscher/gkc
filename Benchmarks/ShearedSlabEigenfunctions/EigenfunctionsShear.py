from pylab import *
from scipy import *
from scipy.special.orthogonal import *
from PyPPL import *
from scipy.optimize import *
from numpy import *
from Execute import *
from FitFunction import *

import mpmath as mp



S = { 'L_n' : 200., 'L_T' :  25.00,  'L_s' : 2000., 'tau' : 1. }

def ShearedSlabDispersion(W, S):
      # Diamagnetic drift frequency
      eta = S['L_n']/S['L_T']
      a = (- S['L_n']**2. / (S['L_s']**2. * W**2.))**0.25
      eps = a**(-2)*(- S['ky']**2. + (1 - W)/(W+(1.+eta)/S['tau']))
      return eps-2*S['l']-1

def ShearedSlabEigenfunction(x,ky, l):

    S['ky'] = ky
    S['l' ] = l

    W = getComplexZero(ShearedSlabDispersion, S)
    
    #res = []
    #for xv in x:
    #mp.mp.dps =50 
    #    a = (- S['L_n']**2. / (S['L_s']**2. * w**2.))**0.25
    #    zeta = mp.mpc(a*xv)
    #    #res.append(complex(  mp.sqrt(a/(mp.sqrt(pi)*(2**l)*mp.factorial(l)))*mp.hermite(l,zeta)*mp.exp(-(zeta/2.)**2.)))
    #    res.append(complex(  mp.sqrt(a/(mp.sqrt(pi)*(2**l)*mp.factorial(l)))*hermitenorm(l)(zeta)*mp.exp(-(zeta/2.)**2.)))
    #return array(res)
    
    a = (- S['L_n']**2. / (S['L_s']**2. * W**2.))**0.25
    zeta = a*x
    return sqrt(a/(sqrt(pi)*2**l*math.factorial(l)))*hermitenorm(l)(zeta)*exp(-(zeta/2.)**2.)


"""

## Plot Eigenfunction
x = linspace(-100./2.,100./2.,101)
for l in range(6):
  phi_v = ShearedSlabEigenfunction(x,0.2,l)
  clf()
  fig = plt.figure(1, [5,4])
  ax = fig.add_subplot(111)
  pl1 = plot(x,real(phi_v), "b-", label="Real Part")
  pl2 = plot(x,imag(phi_v), "r-", label="Imaginary Part")
  xlabel("x")
  legend([pl1, pl2], ("Frequency", "Growthrate"))
  text(0.02,0.02, "Radial Mode Number " + str(l), {'color' : 'r', 'fontsize' : 20}, transform = ax.transAxes)
  savefig("Eigenfunctions_"+str(l) + ".png")



clf()
# Plot Growthrateiis
ky = logspace(-2,1.,101)
for l in range(12):
  w = []
  for kky in ky:
    S['l'] = l
    S['ky'] = kky
    w_e = - S['ky']/S['L_n']
    w.append(w_e*getComplexZero(ShearedSlabDispersion, S, dx=1.e-2, ztol=1.e-8, ftol=1.e-8))
  w = array(w)
  clf()
  pl1 = semilogx(ky, real(w), "b-", label="Frequency")
  ylabel("Frequency")
  xlabel("$k_y$")
  twinx()
  pl2 = semilogx(ky,imag(w), "r-", label="Growthrates")
  ylabel("Growthrates")
  legend((pl1, pl2), ("Frequency", "Growthrate")) 
  savefig("Dispersion_"+str(l)+".png")
"""
