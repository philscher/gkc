
#from HeliosPlot import *
from FitGrowthrates import *
from pylab import *
from matplotlib import patches
from EigenfunctionsShear import *
from scipy.integrate import *
import numpy as np
shearList = []

import mpmath as mp

#from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
#from mpl_toolkits.axes_grid1.inset_locator import mark_inset

X = fileh[0].root.Grid._v_attrs.X
D0 = getDomain(fileh[0])

print  D0['ky'][:8]

"""
for nky in range(len(D0['ky'][:])):
    ky = D0['ky'][nky]
    phi_x = rfft(sum(fileh[0].root.Potential.Electric[:,:,:,-1], axis=0), axis=0)[nky,:]
    plot(X,phi_x,label="Mode $k_y$ : " + str(ky))
    # now built up the eigenfucntions
    X_th   = linspace(min(X),max(X),501)
    phi_th = linspace(0.,0.,501)
    normv = []
    for l in range(len(D0['kx'][:30])):
      #phi_tht, w = phi(X,ky, l)
      phi_tht = phi(X,ky, l)
      fac = simps(real(phi_tht)*phi_x)
      normv.append(fac)
      phi_th += fac * phi(X_th,ky, l)#[0]
    print "Fitten ky ", nky, "   ky : ", ky
    fig = plt.figure(1, [5,4])
    ax = fig.add_subplot(111)
    clf()
    pl1 = plot(X,phi_x, 'g.-', label="Numerical")
    text(0.02,0.02, "Mode " + str(nky) + " $k_y$ : " + str(ky), {'color' : 'r', 'fontsize' : 20}, transform = ax.transAxes)
    #twinx()
    pl2 = plot(X_th,phi_th, 'k-', label="Theory")
    twinx()
    xlim((min(X), max(X)))
    xlabel("X")
    ylabel("phi(x)")
    legend(loc="upper left")
    # this is an inset axes over the main axes
    a = axes([.62, .6, .25, .25], axisbg='#EFEA99')
    bar(array(range(len(normv))),real(normv)/sum(normv))
    #setp(a)#, xticks=[], yticks=[])
    savefig("PlotEigenfunctions_TH3_"+str(nky)+".png")

"""

"""
for nky in range(len(D0['ky'][:])):
  values = []
  for nT  in range(len(fileh[0].root.Potential.phi[1,1,1,:])):

        ky = D0['ky'][nky]
        phi_x = rfft(sum(fileh[0].root.Potential.phi[:,:,:,nT], axis=0), axis=0)[nky,:]
        normv = []
        for l in range(len(D0['kx'][:11])):
            phi_tht = ShearedSlabEigenfunction(X,ky, l)
            fac = simps(real(phi_tht)*real(phi_x))
            normv.append(fac)
        values.append(array(normv)/sum(normv))
        print "Fitten nT : ", nT, " nky : ", nky, " ky : ", ky
  #setp(a)#, xticks=[], yticks=[])"
  values = np.array(values)
  clf()
  for nx in range(len(values[0,:])): plot(range(len(values[:,1])), values[:,nx], label="$H_%i(x)$" % nx)
  legend()
  xlabel("Time")
  ylabel("Weight")
  savefig("Eigenfunctions_ShearedSlab_"+str(nky)+".png")

"""


# make mode growth plot


for nky in range(len(D0['ky'][:])):
  n_modes = 13
  eigenmode_evolution = []
  Time = getTime(fileh[0].root.Potential.Time[:])[:,0]
  print "Fitten  nky : ", nky, " ky : ", D0['ky'][nky]
  for nT  in range(len(Time)):

        ky = D0['ky'][nky]

        phi_x = rfft(sum(fileh[0].root.Potential.phi[:,:,:,nT], axis=0), axis=0)[nky,:]
        
        normv = []
        X = fileh[0].root.Grid._v_attrs.X
        for l in range(n_modes):
            
            phi_tht = ShearedSlabEigenfunction(X,ky,l)
            # map to new basis
            a = simps(real(phi_tht)*real(phi_x))
            normv.append(a)
        eigenmode_evolution.append(normv)
        
  #setp(a)#, xticks=[], yticks=[])"
  eigenmode_evolution = np.array(eigenmode_evolution)

  #calc eigenvalues
  eigenvalue = []
  for n in range(len(eigenmode_evolution.T)):
            clf()
            Y = eigenmode_evolution.T[n]
            S['ky'] = ky
            S['l'] = n
            w_e = - S['ky']* 1./S['L_n']
            num_egv = w_e * imag(getComplexZero(ShearedSlabDispersion, S))
            try :
                semilogy(Time, Y)
                savefig("BalBlu_" + str(n) + ".png") 
            except : pass
            eigenvalue.append((num_egv, fitExpGrowthOptimize(Time,Y,len(Time)*3/4)))
  eigenvalue = array(eigenvalue) 
  clf()
  pl1 = plot(range(n_modes), eigenvalue[:,0], "r.-", label="Theory")
  pl2 = plot(range(n_modes), eigenvalue[:,1], "b^-", label="Measured")
  legend([pl1,pl2], ("Theory", "Measured"))
  xlabel("Time")
  ylabel("$\\Gamma$")
  savefig("Eigenvalues_ky_"+str(nky)+".png")





# make mode growth plot














