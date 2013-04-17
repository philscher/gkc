
from gkc import *
import pylab
import numpy as np


def plotCrossCorrelateValues(fileh5, linLimit=0, numBins=33, mom=0, field=0):

  import scipy.signal
  import cmath
   
  
  #mom   = ['M_{00}', 'M_{20}', 'M_{02}' ] 
  #field = ['\\phi', 'A_{1\\parallel}', 'B_{1\parallel}' ] 


  D = gkcData.getDomain(fileh5) 

  ky_list = D['ky'][1:-1]
  
  Nx = D['Nx']
  
  CrossPhases = zeros( (numBins-1, len(ky_list)) )

  alpha_bin = linspace(-pi, pi, numBins)
  alpha_pl  = linspace(-pi, pi, numBins-1) # Wroong ?!

  # Get non-linear cross-phases
  for m, ky in enumerate(ky_list) :
  
    # Take mean value over x
    #mean_crossphase   = mean(fileh[0].root.Analysis.Flux.CrossPhases[field, mom, 0, m+1,:,-500:-1], axis=0)
    mean_crossphase   = mean(fileh[0].root.Analysis.Flux.CrossPhases[field, mom, 0, m+1,Nx/2-5:Nx/2+5,-1200:-1], axis=0)
    bin_data, alpha   = np.histogram(mean_crossphase, alpha_bin)
    CrossPhases[:,m]  = bin_data

  ####################### Plot data ################################
 
  print shape(alpha), shape(ky_list), shape(CrossPhases)

  pylab.contourf(alpha_pl, ky_list, CrossPhases.T, 100, cmap=pylab.cm.hot)
  #pylab.pcolor(alpha_pl, ky_list, CrossPhases.T, cmap=pylab.cm.hot)
    
  pylab.gca().set_yscale("log") 
    
  pylab.gca().set_xticks([-np.pi, -np.pi/2., 0, np.pi/2., np.pi])
  pylab.gca().set_xticklabels(['$\\pi$', '$-\\frac{\pi}{2}$', '$0$', '$\\pi/2$', '$\\pi$'])
    
  # plot zero phase line
  pylab.plot(np.linspace(0., 0., 101), np.linspace(ky_list[0], ky_list[-1], 101), 'r-', linewidth=2.)
    
  pylab.xlim((-np.pi, np.pi))
  pylab.ylim((ky_list[0], ky_list[-1]))
   
  pylab.xlabel("$\\alpha$")
  pylab.ylabel("$k_y$")


for n, fileh5 in enumerate(fileh):

  fig = gkcStyle.newFigure("Thesis")#, "1:1") 
  
  clf()
  plotCrossCorrelateValues(fileh5, numBins=129, mom=0, field=0)
  savefig("CrossPhases_Phi_Mom_00.png", bbox_inches='tight')
  
  clf()
  plotCrossCorrelateValues(fileh5, numBins=129, mom=1, field=0)
  savefig("CrossPhases_Phi_Mom_20.png", bbox_inches='tight')

  clf()
  plotCrossCorrelateValues(fileh5, numBins=129, mom=2, field=0)
  savefig("CrossPhases_Phi_Mom_02.png", bbox_inches='tight')

  clf()
  plotCrossCorrelateValues(fileh5, numBins=129, mom=0, field=0)
  savefig("CrossPhases_Ap_Mom_00.png", bbox_inches='tight')
  
  clf()
  plotCrossCorrelateValues(fileh5, numBins=129, mom=1, field=0)
  savefig("CrossPhases_Ap_Mom_20.png", bbox_inches='tight')

  clf()
  plotCrossCorrelateValues(fileh5, numBins=129, mom=2, field=0)
  savefig("CrossPhases_Ap_Mom_02.png", bbox_inches='tight')

