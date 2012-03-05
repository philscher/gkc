
def getGrowth(fileh5=fileh[0], pos=(-1,-10), dir='Y'):
    
  D = getDomain(fileh5)

  dPhi = log(fileh5.root.Analysis.PowerSpectrum.Y[0,1:,pos[0]])-log(fileh5.root.Analysis.PowerSpectrum.Y[0,1:,pos[1]])
  dT   = fileh5.root.Analysis.PowerSpectrum.Timing[pos[0]][1]-fileh5.root.Analysis.PowerSpectrum.Timing[pos[1]][1]

  return (D['ky'], dPhi/dT)



def getGrowthEnergy(fileh5=fileh[0], pos=(-1,-10)):
    
  D = getDomain(fileh5)

  dPhi = log(fileh5.root.Analysis.scalarValues.cols.phiEnergy[pos[0]]) - log(fileh5.root.Analysis.scalarValues.cols.phiEnergy[pos[1]])
  dT   = fileh5.root.Analysis.scalarValues.cols.Time[pos[0]] - fileh5.root.Analysis.scalarValues.cols.Time[pos[1]]
     
  # factor 1/2 because we are looking for phi=exp[gamma*t] but calculating from phi^2
  return 0.5*dPhi/dT



