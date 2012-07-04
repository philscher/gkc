from matplotlib import colors, ticker
from scipy.signal import *
from pylab import *
from tables import *
from mpl_toolkits.axes_grid.parasite_axes import HostAxes, ParasiteAxes
import pylab
from Filtering import *

import os
import string

import  FitFunction 

from scipy.ndimage.filters import *
from scipy.ndimage  import *

#from FitGrowthrates import *
import streamlines


import HeliosStyle

from gkcData import *
from gkcAnalysis import *



fileh = []
D     = []

for fileh5 in sys.argv[1:]:
     try:
        fileh.append(openFile(fileh5))
     except:
       print "Cannot open : " , fileh5, " (skipping)"
#     D.append(getDomain(fileh5 = fileh[-1]))


# reload data states, (the simple way)
def reload():
   files= []  
   for file in sys.argv[1:]:
         files.append(openFile(file))


########################################3
#Physics function
##
#   phi(z,y,x)
#     f(v,z,y,x)
#


def plotCFL(fileh5):
  clf()
  semilogy(fileh.root.cfl.cols.time, fileh.root.cfl.cols.Total[:], "r")
  semilogy(fileh.root.cfl.cols.time, fileh.root.cfl.cols.Fx[:]   , "c")
  semilogy(fileh.root.cfl.cols.time, fileh.root.cfl.cols.Fy[:]   , "m")
  semilogy(fileh.root.cfl.cols.time, fileh.root.cfl.cols.Fz[:]   , "g")
  semilogy(fileh.root.cfl.cols.time, fileh.root.cfl.cols.Fv[:]   , "b")

  xlabel("Time")
  ylabel("Possible timestep $\\Delta t$")
  leg = legend(("Total", "$E_y \\frac{\\partial f}{\\partial x}$", "$E_x \\frac{\\partial f}{\\partial y}$", \
          "$E_z \\frac{\\partial f}{\\partial v_\parallel}$", "$v_\\parallel \\frac{\\partial f}{\\partial z}$"), ncol=2)
  leg.draw_frame(0)



#########################################################
# Define some plot functions
 

def contourXTime(fileh5, frame=0, data="temp", species=0, plot3D=False):
  D = getDomain(fileh5)
  
  if(data == "temp"): data = fileh5.root.Analysis.XProperties.Temperature[species,:,1:]

  rect = [0.08, 0.1, 0.75, 0.8 ] 
  conax = axes(rect, axisbg='w')

  X = getTime(fileh5.root.Analysis.XProperties.Timing)[1:,1]
  Y = linspace(0.0, D['Lx'], D['Nx'])
 
  
 # if(filter == "high") : data =  data - gaussian_filter(data, 0.5)
 # if(filter == "low")  :  data =  gaussian_filter(data, 0.5)

  if(plot3D == False):
    contourf(X,Y,data,100, rect=rect)

    xlim((X[0], X[-1]))
    ylim((Y[0], Y[-1]))
  
    xlabel("Time")
    ylabel("X [$\\rho_x$]")
    # set size
    rect = [0.84, 0.1, 0.05, 0.8 ] 
    colax = axes(rect, axisbg='w')
    cb = colorbar(cax=colax)
    cb.set_label("$\phi$")

  if(plot3D == True):
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    import matplotlib.pyplot as plt
    import numpy as np
    clf()
    fig = gcf()
    ax = Axes3D(fig)
    Z =  nan_to_num(data)
    X, Y = meshgrid(X,Y)
    norm = normalize(vmin = Z.min(), vmax=Z.max())
    #pl = ax.plot_surface(X, Y, Z, 300, cmap=cm.jet, norm=norm)
    pl = ax.plot_surface(X, Y, Z, 100,rstride=1, cstride=1, norm=norm, cmap=cm.jet, antialiased=True)
    cset = ax.contourf(X, Y, Z, 100, zdir='z', offset=Z.min())#v_minmax, cmap=grow, norm=norm)
    xlabel("Time")
    ylabel("Position")
    ax.set_zlabel("Temperature")
    fig.colorbar(pl, shrink=0.6)



def plotXProperty(fileh5, time=-1, var='n'):
  D = getDomain(fileh5)
  x = linspace(0.0, D['Lx'], D['Nx'])
  species_name = []
  pl = []
  for s in range(D['Ns']):
    species_name.append(fileh5.root.Species.cols.Name[s])
    if   var == 'n':
        pl.append(plot(x, fileh5.root.Analysis.XProperties.Density[:,s,time]))
    elif var == 'T':
        pl.append(plot(x, fileh5.root.Analysis.XProperties.Temperature[:,s,time]))
    else:
        print "No such varue"
  xlabel("X [$\\rho_x$]")
  if   var == 'n':
        ylabel("Density")
  elif var == 'T':
        ylabel("Temperature")
  
  leg = legend(pl, species_name, loc='lower right')
  leg.draw_frame(0)




def plotHeatKy(fileh5, time=-1, log=0):
  D = getDomain(fileh5)

  for s in range(D['Ns']):
    semilogx(D['ky'],fileh5.root.Analysis.XProperties.HeatKy[1:,s,time], '.-', label="Heat Transport (" + fileh5.root.Species.cols.Name[s]+ ")" )
  leg = legend(loc='lower right', ncol=2)
  leg.draw_frame(0)
  xlabel("ky")
  ylabel("$\\Xi(k_y)$")



