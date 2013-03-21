""" 
Description : Plots and saves the time evolution of the total kinetic energy
              and the total particle number. Note that in order to see energy
              conservation the non-linear Landau damping needs to be included.

License     : GPLv3+
Authors     : The gkc++ Authors

"""

from gkc import *


for fileh5 in fileh:

  Ns = len(fileh5.root.Species.cols.Name[1:])
  
  # Plot Particle Number
  for s in range(Ns): 
    name = fileh5.root.Species.cols.Name[1+s][0]
    semilogy(abs(array(fileh[0].root.Analysis.scalarValues.cols.ParticleNumber)[:,s]), label=name)

  xlabel("$t$")
  ylabel("$\\delta n$")
    
  legend(loc='best').draw_frame(0)
  
  savename = getFileName(fileh5) + "_ParticleConservation.pdf"
  savefig(savename, bbox_inches='tight')

  # Plot Energy Number
  clf()

  for s in range(Ns): 
    name = fileh5.root.Species.cols.Name[1+s][0]
    semilogy(abs(array(fileh5.root.Analysis.scalarValues.cols.KineticEnergy)[:,s]), label=name)
  
  xlabel("$t$")
  ylabel("$\\delta E_{\\rm{kin}}$")
  
  legend(loc='best').draw_frame(0)

  savename = getFileName(fileh5) + "_KineticEnergy.pdf"
  savefig(savename, bbox_inches='tight')

