#from HeliosPlot import *
from FitGrowthrates import *
from pylab import *
from matplotlib import patches
import numpy as np
shearList = []

for i in range(len(fileh)):
    a0 = getFullModeSpectrum(fileh[i], pos=(-30, -1))
    n = argmax(a0[1])
    growthrate = a0[1][1:][n]*sqrt(1837.)*sqrt(2.)
    maxky      = a0[0][1:][n]
    shear      = fileh[i].root.Geometry._v_attrs['Shear'][0]
    shearList.append( (shear,  maxky,  growthrate) )

dtypes = [('Shear', float), ('MaxKy', float), ('Growth', float)]
shearList = np.array(shearList)#, dtype=dtypes)
#shearList.sort(order='Shear')                        # doctest: +SKIP
#print shearList

print shearList
shearList.view('f8,f8,f8').sort(order=['f0'], axis=0) 
print shearList


plot(shearList[:,0], shearList[:,1], "y.-")
xlabel("Shear $\\hat{s}$")
ylabel("Max unstable Mode $k_y$")
twinx()
plot(shearList[:,0], shearList[:,2], "rs-", markersize=10.)
ylabel("Max. Growthrates $k_y$")
patches.Arrow(10., 10., 10., 10., width=1.0, color="y")
patches.Arrow(40., 60., 100., 100., width=1.0, color="r")

savefig("PlotShear_1231.png")


