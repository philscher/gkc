import pylab
import sys
import gkcLinear 


################## gs2 Data Analysis #############
import Scientific.IO.NetCDF as NetCDF

# Open File
gs2_file =  NetCDF.NetCDFFile(sys.argv[1],'r')

# Get time evolution of 
phi2 = gs2_file.variables['phi2_by_ky'][:]
T    = gs2_file.variables['t'][:]

ky    = gs2_file.variables['ky'][:]

gamma_gs2 = gkcLinear.getGrowthrate(T, phi2.T, 10, -1)



################ gkc++ Data Analysis ###############

import tables

gkc_file = tables.openFile(sys.argv[2])

omega = getFrequencyGrowthrates(gkc_file, 10, stop=-1)

frequency, gamma = real(omega), imag(omega)


################### Plot stuff ####################

pylab.semilogx(ky, gamma_gs2, label='gs2')
pylab.semilogx(ky, gamma_gkc, label='gkc')

pylab.savefig("LinearGrowthrates.png", bbox_inches='tight')
