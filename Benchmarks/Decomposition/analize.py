from numpy import *
from tables import *
import sys

# Note the first file is the reference file

fileh0 = openFile(sys.argv[1])



for file_name in sys.argv[2:]:
  try :
    fileh = openFile(file_name)
    
    #print colored(file_name, 'red'), "   : " "Phi2 : %.2f" % ((sum(fileh.root.Potential.phi[:,:,:,-1]**2)-sum(fileh0.root.Potential.phi[:,:,:,-1]**2)) / sum(fileh0.root.Potential.phi[:,:,:,-1]**2)), \
    print file_name, "   : " "Phi2 : %.4e" % ((sum(fileh.root.Potential.phi[:,:,:,-1]**2)-sum(fileh0.root.Potential.phi[:,:,:,-1]**2)) / sum(fileh0.root.Potential.phi[:,:,:,-1]**2)), \
                     "Ap2 : %.4e" % ((sum(fileh.root.Potential.Ap[:,:,:,-1]**2)-sum(fileh0.root.Potential.Ap[:,:,:,-1]**2)) / sum(fileh0.root.Potential.Ap[:,:,:,-1]**2)), \
                     " , PowX (phi): %.4e" % ((sum(fileh.root.Analysis.PowerSpectrum.X[0,:,-1]**2)-sum(fileh0.root.Analysis.PowerSpectrum.X[0,:,-1]**2))/sum(fileh0.root.Analysis.PowerSpectrum.X[0,:,-1]**2)), \
                     " , PowY (phi): %.4e" % ((sum(fileh.root.Analysis.PowerSpectrum.Y[0,:,-1]**2)-sum(fileh0.root.Analysis.PowerSpectrum.Y[0,:,-1]**2))/sum(fileh0.root.Analysis.PowerSpectrum.Y[0,:,-1]**2)), \
                     " Energy (kin.) %.4e" % ((fileh.root.Analysis.scalarValues.cols.KineticEnergy[-1] - fileh0.root.Analysis.scalarValues.cols.KineticEnergy[-1])/fileh0.root.Analysis.scalarValues.cols.KineticEnergy[-1])[0], \
                     " Energy (phi.) %.4e" % ((fileh.root.Analysis.scalarValues.cols.phiEnergy[-1] - fileh0.root.Analysis.scalarValues.cols.phiEnergy[-1])/fileh0.root.Analysis.scalarValues.cols.phiEnergy[-1])
    fileh.close() 
#    try :
#       print                " , PowX (phi): " , sum(fileh.root.Analysis.PowerSpectrum.X[1,:,-1]**2)-sum(fileh0.root.Analysis.PowerSpectrum.X[1,:,-1]**2), \
#                       " , PowY (phi): " , sum(fileh.root.Analysis.PowerSpectrum.Y[1,:,-1]**2)-sum(fileh0.root.Analysis.PowerSpectrum.Y[1,:,-1]**2)
#    except : pass
  except: print file_name, "  Failed !!!!"
fileh0.close()
