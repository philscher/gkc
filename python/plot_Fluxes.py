""" 
Description : Write me please ! 

License     : GPLv3+
Authors     : The gkc++ Authors

"""

from gkc import *

for fileh5 in fileh:
  
  hflux = array(fileh5.root.Analysis.scalarValues.cols.HeatFlux)
  pflux = array(fileh5.root.Analysis.scalarValues.cols.ParticleFlux)
  time  =       fileh5.root.Analysis.scalarValues.cols.Time[:]

  Nq    = len(hflux[0,0,:])
  Ns    = len(hflux[0,:,0])
      
  q_leg = ['\\phi', 'A_{1\\parallel}', 'B_{1\\parallel}']

  ############ Plot particle fluxes

  for q in range(Nq):

    clf()
    

    for s in range(Ns): 
      # Use only first letter due to bug
      name = fileh5.root.Species.cols.Name[1+s][0]

      # Weird scaling issue
      mass = fileh5.root.Species.cols.Mass[1+s]

      # Plot Particle Number
      plot(time, mass * pflux[:,s,q], label=name)


    xlabel("$t$")
    ylabel("$\\Gamma_{%s}$" % q_leg[q])
    
    legend(loc='best').draw_frame(0)

    xlim((min(time), max(time) * 1.02))

    
    plot(linspace(min(time), max(time) * 1.02, 128), zeros(128), color='#333333', linewidth=1.0)
 
    savename = getFileName(fileh5) + "_ParticleFlux_%i.pdf" % q
    savefig(savename, bbox_inches='tight')

  ############ Plot heat fluxes

  for q in range(Nq):

    clf()

    for s in range(Ns): 
      # Use only first letter due to bug
      name = fileh5.root.Species.cols.Name[1+s][0]
      
      # Weird scaling issue
      mass = fileh5.root.Species.cols.Mass[1+s]

      # Plot Particle Number
      plot(time, mass * hflux[:,s,q], label=name)

    xlabel("$t$")
    ylabel("$\\chi_{%s}$" % q_leg[q])
    
    legend(loc='best').draw_frame(0)


    #plot Zero line
    plot(linspace(min(time), max(time) * 1.02, 128), zeros(128), color='#333333', linewidth=1.0)

    xlim((min(time), max(time) * 1.02))

    
    savename = getFileName(fileh5) + "_HeatFlux_%i.pdf" % q
    savefig(savename, bbox_inches='tight')


