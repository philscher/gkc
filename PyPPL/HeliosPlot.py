#from scipy.io import read_array
from matplotlib import colors, ticker
from scipy.signal import *
from pylab import *
from tables import *
from mpl_toolkits.axes_grid.parasite_axes import HostAxes, ParasiteAxes
import pylab
from Filtering import *

import os
import string
#import GrowthrateFit

import  FitFunction 

from scipy.ndimage.filters import *
from scipy.ndimage  import *

from FitGrowthrates import *
import streamlines

from HeliosStyle import *


fileh = []
D     = []

for fileh5 in sys.argv[1:]:
     try:
        fileh.append(openFile(fileh5))
     except:
       print "Cannot open : " , fileh5, " (skipping)"
#     D.append(getDomain(fileh5 = fileh[-1]))










def Lapl3D(A) :
  B = gradient(A)
  return diff(B[0]) + diff(B[1]) + diff(B[2])

def getTime(timing):
  timing_new = []
  for t in timing:
    timing_new.append((t[0], t[1]))
  return array(timing_new)

def getDomain(fileh5=fileh[0]):
  Ny = fileh5.root.Grid._v_attrs.Nky 
  Ly = fileh5.root.Grid._v_attrs.Ly
  Nx = fileh5.root.Grid._v_attrs.Nx 
  Lx = fileh5.root.Grid._v_attrs.Lx
  Nz = fileh5.root.Grid._v_attrs.Nz 
  Lz = fileh5.root.Grid._v_attrs.Lz

  species = []
  for s in range(fileh5.root.Grid._v_attrs.Ns):
     species.append((fileh5.root.Species.cols.Name[s], fileh5.root.Species.cols.Charge[s], fileh5.root.Species.cols.Mass[s]))
     # for modes don't included DC component
  D = {    'Nx' : fileh5.root.Grid._v_attrs.Nx, \
           'Nky' : fileh5.root.Grid._v_attrs.Nky, \
           'Nz' : fileh5.root.Grid._v_attrs.Nz, \
           'Nv' : fileh5.root.Grid._v_attrs.Nv, \
           'Nm' : fileh5.root.Grid._v_attrs.Nm, \
           'Ns' : fileh5.root.Grid._v_attrs.Ns, \
           'Nkx': fileh5.root.Grid._v_attrs.Nx/2+1, \
           'Ny': fileh5.root.Grid._v_attrs.Nky*2-2, \
           'Nkp': fileh5.root.Grid._v_attrs.Nz/2+1, \
           'Lx' : fileh5.root.Grid._v_attrs.Lx, \
           'Ly' : fileh5.root.Grid._v_attrs.Ly, \
           'Lz' : fileh5.root.Grid._v_attrs.Lz, \
           'Lv' : fileh5.root.Grid._v_attrs.Lv,  \
           'Lm' : fileh5.root.Grid._v_attrs.Lm,  \
           'Tscale' : sqrt(1837.) / (2. *sqrt(2.) * 1.e-3), \
           'X'     : fileh5.root.Grid._v_attrs.X, \
           'Z'     : fileh5.root.Grid._v_attrs.Z, \
           'V'     : fileh5.root.Grid._v_attrs.V, \
           'M'     : fileh5.root.Grid._v_attrs.M, \
           'kx'   : arange(1, Nx/2+1) * 2 * pi / Lx,\
           'ky'   : arange(1, Ny) * 2 * pi / Ly,\
           'kp'   : arange(1, Nz/2+1) * 2 * pi / Lz,\
           'species': species,\
           'Debye2' : fileh5.root.Plasma._v_attrs.Debye2[0],\
           }
  D['Y'] = linspace(0., D['Ly'], D['Ny'])
  D['TScale'] =  D['Lv']/D['Lz']

  return D

            


def gaussian(height, center_x, center_y, width_x, width_y):
      """Returns a gaussian function with the given parameters"""
      width_x = float(width_x)
      width_y = float(width_y)
      return lambda x,y: height*exp(-(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)



#files = []
#for file in sys.argv[1:]:
#    files.append(openFile(file))

def plotFormat():
  rect = [0.08, 0.1, 0.75, 0.8 ] 
  conax = axes(rect, axisbg='w')
  gcf().set_size_inches( (16,12) )
  rcParams['legend.loc'] = 'best'

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

def Temperature(time, species=0, fileh5=fileh[0]):
  f = fileh5.root.Phasespace.Data[species,0,:,:,:,:,time]
  D = getDomain(fileh5)
  v = D['V']
  
  temperature = zeros( shape(f[0,:,:,:]))

  for z in range(len(f[0,:,0,0])):
    for y in range(len(f[0,0,:,0])):
      for x in range(len(f[0,0,0,:])):
        temperature[z,y,x] = (f[:,z,y,x]* v**2).sum() / f[:,z,y,x].sum()
       # .sum(axis=0)[z,y,x]


  return temperature


# return density


def plotCFL(fileh5=fileh[0]):
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



def Density(time, species=0,fileh5=fileh[0]):

  f = fileh5.root.Phasespace.Data[species,0,:,:,:,:,time]

  density = zeros(shape(f[0,:,:,:]))
  for z in range(len(f[0,:,0,0])):
    for y in range(len(f[0,0,:,0])):
      for x in range(len(f[0,0,0,:])):
        density[z,y,x] =  f[:,z,y,x].sum()
        #.sum(axis=0)[z,y,x]
        #density[z,y,x] =  sum(file.root.PhaseSpace[i,j,k]
  return density


def KolmogorovXY(filter='none', fileh5=fileh[0], frame=-1):
  X = fileh.root.Analysis.PowerSpectrum.X[:,frame]
  Y = fileh.root.Analysis.PowerSpectrum.Y[:,frame]
  lenMax = min(len(X), len(Y))

  KR = linspace(1.e-6,1.e-6,sqrt(2.) * lenMax+1)
  
  for x in range(lenMax):
    for y in range(lenMax):
      r = int(sqrt(x**2 + y**2))
      KR[r] = X[x] + Y[y] 

  print KR
  loglog(KR,'.-')



def Kolmogorov2D(time, filter='none', fileh5=fileh[0]):
  phi = fileh5.root.Potential.Potential.phi[:,:,:,time]
  E = gradient(phi)
  energy = (E[0]*E[0] + E[1]*E[1] + E[2]*E[2])
  energy = sum(energy, axis=0)
  if(filter == "high")  : energy =  energy - gaussian_filter(energy, 2.0, cvar=mean(energy))
  if(filter == "low")   : energy =  gaussian_filter(energy, 2.0, cvar=mean(energy))
  energy_k = fftn(energy[:,:])/(energy.size)
  kgorov = abs(energy_k)**2
  lenx = len(energy_k[:,0])
  leny = len(energy_k[0,:])
  len_max = min(lenx, leny)
  kmax = sqrt( 2 * len_max**2) 
  k = linspace(0.0,0.0,kmax+1)

  for x in range(len_max):
    for y in range(len_max):
          k_x = int(sqrt(x**2 + y**2))
          k[k_x] += kgorov[y,x]/(4*pi*0.5**2)
  k[0] = sum(energy_k)
  print k
  loglog(k,'.-')




def Kolmogorov3D(time, filter='none', fileh5=fileh[0]):
  #temp0= Temperature(0, fileh5)
  #temp = Temperature(time, fileh5)
  #dens0= Density(0, fileh5)
  #dens = Density(time, fileh5)

  phi = fileh5.root.Potential.phi[:,:,:,time]
  if(filter == "high")  : phi =  phi - gaussian_filter(phi, 0.6)
  if(filter == "low")   : phi =  gaussian_filter(phi, 0.6)
  E = gradient(phi)
  energy = (E[0]*E[0] + E[1]*E[1] + E[2]*E[2])

  energy_k = fftn(energy)/(energy.size)
  
  kgorov = abs(energy_k)**2

  lenx = len(kgorov[:,0,0])
  leny = len(kgorov[0,:,0])
  lenz = len(kgorov[0,0,:])

  len_max = min(lenx, leny, lenz)
  print len_max


  kmax = sqrt( 3 * len_max**2) 
  k = linspace(0.0,0.0,kmax+1)

  for x in range(len_max):
    for y in range(len_max):
       for z in range(len_max):
          k_x = int(sqrt(x**2 + y**2 + z**2))
          k[k_x] += kgorov[z,y,x]/(4*pi*0.5**2)
  k[0] = sum(kgorov)

  loglog(k,'.-')
  ylabel("$log E(k)$")
  xlabel("$log k$")
  title("Kolmogorov spectrum")


#########################################################
# Define some plot functions
 
def contourXY(frame=0, data="phi", filter='None',Z = 0, vminmax = 0, fileh5=fileh[0], clear=True, stream=True, species=0):
  fig = subplot(111)
  print "DATA : " , data
  if(data == "phi"): data = fileh5.root.Potential.phi[Z,:,:,frame]
  if(data == "Ap"):  data = fileh5.root.Potential.Ap[Z,:,:,frame]
  if(data == "n"):   data = fileh5.root.Moments.TemperatureParallel[species,Z,:,:,frame]
  if(data == "T"):   data = fileh5.root.Moments.NumberDensity[species,Z,:,:,frame]
  if(data == "charge") :
    data = Lapl3D(fileh5.root.Potential.phi[:,:,:,frame])[Z,:,:]
  print shape(data)
  D = getDomain(fileh5)

  rect = [0.08, 0.1, 0.75, 0.8 ] 
  conax = axes(rect, axisbg='w')

  x = D['X'] 
  y = D['Y'] 
  
  if(filter == "high") : data =  data - gaussian_filter(data, 0.5)
  elif(filter == "low")  :  data =  gaussian_filter(data, 0.5)
  else : print "NoSuchFilter"

  if clear == True: clf()
  contourf(x,y,data,100, rect=rect)
  cb = colorbar()#cax=colax)
  """  
  if stream == True:
    data = fileh5.root.Potential.Ap[Z,:,:,frame]
    dA_dy, dA_dx = gradient(data)
    X,Y = meshgrid(x,y)
    streamlines.Streamlines(x,y, dA_dy, -dA_dx, maxLen=500.).plotArrows(ax=conax)
  """
  xlim((x[0], x[-1]))
  ylim((y[0], y[-1]))
  
  title("phi potential at Time Step = " + str(fileh5.root.Potential.Time[frame][0]) + " Time = "  + str(fileh5.root.Potential.Time[frame][1]))
#      ) +\                         " (Filter : " + filter + ")") # +  " (at Z = %i) " % Z) 
  xlabel("X")
  ylabel("Y")
  # set size
  cb.set_label("$\phi$")


def contourXYFFT(frame=0, data="phi", filter='None',Z = 0, vminmax = 0, fileh5=fileh[0], clear=True, stream=True, species=0):
  fig = subplot(111)
  print "DATA : " , data
  if(data == "phi"): data = fileh5.root.Potential.phi[Z,:,:,frame]
  if(data == "Ap"):  data = fileh5.root.Potential.Ap[Z,:,:,frame]
  if(data == "n"):   data = fileh5.root.Moments.TemperatureParallel[species,Z,:,:,frame]
  if(data == "T"):   data = fileh5.root.Moments.NumberDensity[species,Z,:,:,frame]
  if(data == "charge") :
    data = Lapl3D(fileh5.root.Potential.phi[:,:,:,frame])[Z,:,:]
  
  D = getDomain(fileh5)

  rect = [0.08, 0.1, 0.75, 0.8 ] 
  conax = axes(rect, axisbg='w')


  # back transform to real space
  # Resize to have higher resolution in y
  print "Shape : " , shape(data)
  #Psi_1k = resize(data,(D['Nx'], max(128, D['Ny'])))
  #Psi_1k[nky:,:] = 0.

  # Fourier back transform c2r, renormaliza
  data = 128.*ifft(data, axis=0)


 
  x = D['X'] 
  y = D['Y'] 
  print "Shape3 : " , shape(data), " X : ", len(x), "  " , len (y)
  
  if(filter == "high") : data =  data - gaussian_filter(data, 0.5)
  elif(filter == "low")  :  data =  gaussian_filter(data, 0.5)
  else : print "NoSuchFilter"

  if clear == True: clf()
  contourf(x,y,data,100, rect=rect)
  cb = colorbar()#cax=colax)
  """  
  if stream == True:
    data = fileh5.root.Potential.Ap[Z,:,:,frame]
    dA_dy, dA_dx = gradient(data)
    X,Y = meshgrid(x,y)
    streamlines.Streamlines(x,y, dA_dy, -dA_dx, maxLen=500.).plotArrows(ax=conax)
  """
  xlim((x[0], x[-1]))
  ylim((y[0], y[-1]))
  
  title("phi potential at Time Step = " + str(fileh5.root.Potential.Time[frame][0]) + " Time = "  + str(fileh5.root.Potential.Time[frame][1]))
#      ) +\                         " (Filter : " + filter + ")") # +  " (at Z = %i) " % Z) 
  xlabel("X")
  ylabel("Y")
  # set size
  cb.set_label("$\phi$")



























def contourXTime(frame=0, data="temp", fileh5=fileh[0], species=0, plot3D=False):
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



def plotXProperty(fileh5=fileh[0], time=-1, var='n'):
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



def plotHeatFlux(fileh5=fileh[0], log=0, opt="", legend=1):
  D = getDomain(fileh5)

  if log==0:
    for s in range(D['Ns']):
        plot(fileh5.root.Analysis.scalarValues.cols.Time[3:], fileh5.root.Analysis.scalarValues.cols.HeatFlux[3:][:,s], label="Heat Flux (" + fileh5.root.Species.cols.Name[s]+ ")") 
  else:
    for s in range(D['Ns']):
        semilogx(fileh5.root.Analysis.scalarValues.cols.Time[3:], fileh5.root.Analysis.scalarValues.cols.HeatFlux[3:][:,s], label="Heat Flux (" + fileh5.root.Species.cols.Name[s]+ ")") 
  if legend==1:
    leg = legend(loc='lower right', ncol=2)
    leg.draw_frame(0)
  xlabel("Time")
  ylabel("Heat Flux")

def plotParticleFlux(fileh5=fileh[0], log=0, opt=""):
  D = getDomain(fileh5)
  for s in range(D['Ns']):
    plot(fileh5.root.Analysis.scalarValues.cols.Time[3:], fileh5.root.Analysis.scalarValues.cols.ParticleFlux[3:][:,s], label="Heat Flux (" + fileh5.root.Species.cols.Name[s]+ ")") 
  leg = legend(loc='lower right', ncol=2)
  leg.draw_frame(0)
  xlabel("Time")
  ylabel("Particle Flux")



def plotEnergyEvolution(fileh5=fileh[0], log=0, opt="", fak=1., delta=False):
  D = getDomain(fileh5)


  if(log == 1) : plotf = semilogy
  else         : plotf = plot

  fak = 1.
  if(log == 1) : fak = -1
  if(log == 1) : myf = abs
  else         : myf = lambda x : x

  plotf(fileh5.root.Analysis.scalarValues.cols.Time[3:], myf(fileh5.root.Analysis.scalarValues.cols.phiEnergy[3:]), opt, label="Field Energy")
  for s in range(D['Ns']):
    plotf(fileh5.root.Analysis.scalarValues.cols.Time[3:], myf(fileh5.root.Analysis.scalarValues.cols.KineticEnergy[3:][:,s]) , label="Kinetic Energy (" + fileh5.root.Species.cols.Name[s]+ ")" )
  plotf(fileh5.root.Analysis.scalarValues.cols.Time[3:], myf(myf(fileh5.root.Analysis.scalarValues.cols.phiEnergy[3:])+fak*myf(sum(fileh5.root.Analysis.scalarValues.cols.KineticEnergy[3:], axis=1))) , label="Total Energy")
  leg = legend(loc='lower right', ncol=2)
  leg.draw_frame(0)
  xlabel("Time")
  ylabel("Energy")

  if delta == True:
    twinx()
    def minmax(x):
      if x > 1.    : return 1.
      if x < 1.e-4 : return 1.e-4
      return x

    dE = abs(myf(fileh5.root.Analysis.scalarValues.cols.FieldEnergy[3:])+fak*myf(sum(fileh5.root.Analysis.scalarValues.cols.KineticEnergy[3:], axis=1)))/abs(myf(sum(fileh5.root.Analysis.scalarValues.cols.KineticEnergy[3:], axis=1)))
    semilogy(fileh5.root.Analysis.scalarValues.cols.Time[3:], map(minmax, dE), 'k', linewidth=3.)
    #plotf(fileh5.root.Analysis.scalarValues.cols.Time[3:], dE)
    ylabel("$\\frac{E_{kin}+E_{\\phi}}{E_{kin}} \%$")
    ylim((0.8e-4, 1.4))


def plotParticleEvolution(fileh5=fileh[0], log=0, fak=1.):
  D = getDomain(fileh5)
  if(log == 1) : plotf = semilogy
  else         : plotf = plot
    
  for s in range(D['Ns']):
    plotf(fileh5.root.Analysis.scalarValues.cols.Time[3:], abs(fileh5.root.Analysis.scalarValues.cols.ParticleNumber[3:][:,s]) , label="Number (" + fileh5.root.Species.cols.Name[s]+ ")" )
  if D['Ns'] > 1:
    charge = 0.
    for s in range(D['Ns']):
       charge = charge + fileh5.root.Species.cols.Charge[s]*fileh5.root.Analysis.scalarValues.cols.ParticleNumber[3:][:,s]
    plotf(fileh5.root.Analysis.scalarValues.cols.Time[3:], abs(charge) , label="Charge")
  
  leg = legend(loc='lower right', ncol=2)
  leg.draw_frame(0)
  xlabel("Time")
  ylabel("Number")


def plotHeatKy(fileh5=fileh[0], time=-1, log=0):
  D = getDomain(fileh5)

  for s in range(D['Ns']):
    semilogx(D['ky'],fileh5.root.Analysis.XProperties.HeatKy[1:,s,time], '.-', label="Heat Transport (" + fileh5.root.Species.cols.Name[s]+ ")" )
  leg = legend(loc='lower right', ncol=2)
  leg.draw_frame(0)
  xlabel("ky")
  ylabel("$\\Xi(k_y)$")



def plotCorrelate(A,B,timeStep=-1, fileh5=fileh[0]):
    D = getDomain(fileh5)
#    A = fileh5.root.Potential.phi[:,:,:,timeStep]
#    B = fileh5.root.Potential.phi[:,:,:,timeStep]

    # first correlate then average
    A = sum(A, axis=0)
    B = sum(B, axis=0)
    Corr = scipy.signal.correlate(A, B)
  
    #Corr = sum(Corr, axis=0)
    x_n  = list(-(D['kx']))
    x_n.reverse()
    x = x_n + list((D['kx'][:]))
    y = [ 0 ] + list((D['ky']))
    print "Shape Ap : ", shape(Corr), " x : " , len(x), "   y : " , len(y)
    contourf(Corr[:,:].T,100)#, locator=ticker.LogLocator())
    colorbar()
    print x 
    print y

    #xlim((x[0], x[-1]))
    #ylim((y[0], y[-1]))
    print Corr

def plotSpectrum(timeStep=-1, fileh5=fileh[0], doLog=False):
    A = fileh5.root.Potential.phi[:,:,:,timeStep]
    D = getDomain(fileh5)
    #Average over z axis
    A2 = sum(A, axis=0)
    #Perform FFT
    Ak = rfft(A2.T)
    Ap = (abs(Ak))**2
    print "Shape Ap : ", shape(Ap)
    #plot
    #subplot(121)
    #contourf(D['ky'], D['kx'], pow2(real(Ak[0:lenx,0:leny])),100, locator=ticker.LogLocator())
    #subplot(122)
    
    x_n  = list(-(D['kx']))
    x_n.reverse()
    x = x_n + list((D['kx'][:]))
    y = [ 0 ] + list((D['ky']))
    if(doLog): x = log(x)
    if(doLog): y = log(y)
    print "Shape Ap : ", shape(Ap), " x : " , len(x), "   y : " , len(y), Ap
    contourf(x, y, Ap[:,:].T,100, locator=ticker.LogLocator())
    print x 
    print y
    xlim((x[0], x[-1]))
    ylim((y[0], y[-1]))
    colorbar()
    


def plotPowerSpectrumsGraphYZ(fileh5=fileh[0]):
    #grow_Y = GrowthrateFit.getAnalysis.scalarValues.fileh5, dir="Y") 
    #grow_Z = GrowthrateFit.getAnalysis.scalarValues.fileh5, dir="Z")
    plotFormat()

    clf()
    #p1 = plot(grow_Y,'c')
    ylabel("Growthrate")
    xlabel("Mode")
    twinx()
    #p2 = plot(grow_Z,'m')
    ylabel("Growthrate")
    legend((p1,p2),("Ky", "kp"))
    title("Analysis.scalarValues.of Mode")



def HeatFluxMode(timeStep=-1, fileh5=fileh[0], doLog=False):
    D = getDomain(fileh5)
    Q = []
    for s in range(D['Ns']):
#        A = fileh5.root.Moments.HeatFluxParallel[s,:,:,:,timeStep]
        A = fileh5.root.Moments.HeatFluxOrthogonal[s,:,:,:,timeStep]
        #Average over z axis
        A2 = sum(A, axis=0)
        A2 = sum(A2, axis=1)
        #Perform FFT
        Ak = rfft(A2.T)
        #print "Shape Ap : ", shape(Ap)
        Q.append(Ak[1:])
        print shape(D['ky'])
        print shape(Ak)
        #semilogx(D['ky'], Ak[1:], label=fileh5.root.Species.cols.Name[s])
        #semilogx(D['ky'], Ak)

    Q_sum = sum(abs(Q))
    for s in range(D['Ns']):
        semilogx(D['ky'], Q[s]/Q_sum, label=fileh5.root.Species.cols.Name[s])

    #plot
    #subplot(121)
    leg = legend(loc='lower right', ncol=2)
    leg.draw_frame(0)
    xlabel("$k_y")
    ylabel("Q(k_y)")


def plotTurbulenceTime(dir='Y', pos=(-2,-1), fileh5=fileh[0], doFit='False', posT=(1,-1)):
  
  D = getDomain(fileh5)
  data = fileh5.root.Analysis.PowerSpectrum.Y[1:,:]
  T = getTime(fileh5.root.Analysis.PowerSpectrum.Timing)[:,1]
  timeEvolution = []
  for step in range(len(data[0,:])):
    timeEvolution.append(data[:,step]/sum(data[:,step]))

  contourf(log10(D['ky']), T,timeEvolution, 250, locator=ticker.LogLocator(), cmap=cm.jet)
  colorbar()


def Turbulence(A, time):
  A = A
  Ak = fft(A)
  # sum axis
  Ak = Ak.sum(axis=0)
  Ak = Ak.sum(axis=1)
  semilogx(abs(Ak)**2)



def plotTurbulenceSpectra(dir='Y', pos=(-2,-1), fileh5=fileh[0], doFit='False', posT=(1,-1), field=0):
    D = getDomain(fileh5)
    if(dir == 'X'):
      data = sum(fileh5.root.Analysis.PowerSpectrum.X[field,0:,pos[0]:pos[1]], axis=1)/abs(pos[1]-pos[0])
      loglog(D['kx'], data/sum(data), '.-')
      xlabel("$k_x$")
    elif(dir == 'Y'):
      data = sum(fileh5.root.Analysis.PowerSpectrum.Y[field, 0:,pos[0]:pos[1]], axis=1)/abs(pos[1]-pos[0])
      #data = data * (D['ky'])**2
      print sum(data), shape(data[1:]), shape(D['ky']), D['ky']
      loglog(D['ky'], data[1:]/sum(data), 'b.-')
      # Plot Zonal flow seperately
      loglog(D['ky'][0], data[0], 'ko', markersize=8.)
      xlim((min(D['ky']/1.2), 1.2*max(D['ky'])))
      xlabel("$k_y$")
    elif(dir == 'Z'):
      data = sum(fileh5.root.Analysis.PowerSpectrum.Z[field, 0:,pos[0]:pos[1]], axis=1)/abs(pos[0]-pos[1])
      loglog(D['kp'], data[1:], '.-')
      loglog(D['kp'][0], data[0], 'o')
      xlabel("$k_z$")
      
    ylabel("$|\\phi_k(k_y)|^2$")
        
    if doFit == 'True' :
        pos_a = posT[0]
        pos_b = posT[1]
        fitfunc = lambda p, x: p[0]*x + p[1] # Target function
        errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
        p0 = [1.0, 1.0, 1.0] # Initial guess for the parameters
        
        p1, success = optimize.leastsq(errfunc, p0[:], args=(log(D['ky'][pos_a:pos_b]), log(data[pos_a:pos_b])))


        # C'mon baby wanna see you again
        print p1
        loglog(D['ky'][pos_a:pos_b], (D['ky'][pos_a:pos_b]-D['ky'][0])**p1[0],'r', linewidth=9.)
        #text(D['ky'][5], data[5], "$\\propto k_y^{%2.f}$" % p1[0], ha="center", family=font, size=14)
        text(log(D['ky'][5]), log(data[5]), "$\\propto k_y^{%2.1f}$" % p1[0], ha="center", size=14)

def plotInstantGrowth(numModes=10, dir='Y', fileh5=fileh[0], order=11, filterType='hanning', field=0):
    D = getDomain(fileh5) 
    T = getTime(fileh5.root.Analysis.PowerSpectrum.Time)[:,1]

    #need custom gradient as scipy does not support variable step size
    # First ist forward difference
    # Last central difference
    def Cgradient(F, X):
        df = zeros(len(F))
        N  = len(F)
        df[0] = (F[1] -F[0])/(X[1]-X[0])
        for n in arange(1,N-1): df[n] = (F[n+1] -F[n-1])/(X[n+1]-X[n-1])
        df[N-1] = (F[N-1] -F[N-2])/(X[N-1]-X[N-2])

        return 2.*df

    if(dir == 'X'):

      grad = getGradient(log(fileh5.root.Analysis.PowerSpectrum.X[field,:numModes,:].T))
      
      plot(T, grad)
      legend_list = []
      for i in range(len(fileh5.root.Analysis.PowerSpectrum.X[field,:numModes,0])):
        legend_list.append("kx = %i" % i)
      leg = legend(legend_list, loc='best', ncol=4)
      leg.draw_frame(0)
    
    if(dir == 'Y'):
      for m in numModes:
        Var =   fileh5.root.Analysis.PowerSpectrum.Y[field,m,:]
        gamma = array(Cgradient(log10(Var), T))
        gamma = smoothFilter(gamma, filterType, window_len=order)[0:-order+1]
        print shape(T) , " g : " , shape(gamma)
        plot(T, gamma, label='ky = %.2f' % (D['ky'][m]), color='#'+colors_216[(16*m+1) % 215])
      leg = legend(loc='best', ncol=4)
      #plot(T[10+(order-1)/2:-10-(order-1)/2], grad.T)

      

      #grad = getGradient(fileh5.root.Analysis.PowerSpectrum.Y[field,:numModes,:])
      #plot(T[10+(order-1)/2:-10-(order-1)/2], grad.T)
    
      #legend_list = []
      #for i in range(len(fileh5.root.Analysis.PowerSpectrum.Y[field,:numModes,0])):
      
      leg.draw_frame(0)
    
    
    if(dir == 'Z'):
        scale = fileh5.root.Grid._v_attrs.Lz/(2. * pi)
        plot(gradient(T,scale*fileh5.root.Analysis.PowerSpectrum.Z[:numModes,:].T))
    
        legend_list = []
        for i in range(len(fileh5.root.Analysis.PowerSpectrum.Z[:numModes,0])):
            legend_list.append("n = %i" % i)
        leg = legend(legend_list, loc='lower right', ncol=2)
        leg.draw_frame(0)
    
    
    xlabel("Time")
    ylabel("Mode Growth $\gamma$")
        

def plotTimePower(fileh5=fileh[0], numModes=[], dir='Y', field=0, stride=1, leglabel='k'):
    #fig = figure()#frameon=False)
    #ax = fig.add_axes((0.35,0.2,0.55,0.7))
    #DefaultSize = gcf().get_size_inches()
    #gcf().set_size_inches( (DefaultSize[0]*2.2, DefaultSize[1]) )
                       

    
    #xlabel("Time")
    #xlim((fileh5.root.Analysis.scalarValues.cols.Time[0], fileh5.root.Analysis.scalarValues.cols.Time[-1]))
    #twiny()
    T = getTime(fileh5.root.Analysis.PowerSpectrum.Time)[2:,1]
    #T1 = getTime(fileh5.root.Analysis.PowerSpectrum.Timing)[:,1]
    #ylim((T1[0], T1[-1]))
    #twiny()
    clf()
    if(dir == 'X'):
      pl = semilogy(T, fileh5.root.Analysis.PowerSpectrum.X[field,:numModes,2:].T)
      legend_list = []
      for i in range(len(fileh5.root.Analysis.PowerSpectrum.X[field, :numModes,0])):
        legend_list.append("kx = %i" % i)
      leg = legend(legend_list, loc='lower right', ncol=2)
      leg.draw_frame(0)
    
    if(dir == 'Y'):
      scale = fileh5.root.Grid._v_attrs.Ly/(2. * pi)
      
      legend_list = []
      if numModes == []: 
        pl = semilogy(T, (scale*fileh5.root.Analysis.PowerSpectrum.Y[field, :,2:][::stride,:]).T)
        for i in range(len(fileh5.root.Analysis.PowerSpectrum.Y[field, :,0][::stride])):
            if(leglabel == 'm') : legend_list.append("m = %i" % (i*stride))
            else                : legend_list.append("ky = %.1f" % (i*stride / scale)) 
      else:
        for m in numModes:
            pl = semilogy(T, (scale*fileh5.root.Analysis.PowerSpectrum.Y[field, m,2:]).T)
            if(leglabel == 'm') : legend_list.append("m = %i" % m)
            else                : legend_list.append("ky = %.1f" % (m / scale)) 

    
      leg = legend(legend_list, loc='upper left', ncol=4, mode="expand")
      leg.draw_frame(0)
    
     
    
    if(dir == 'Z'):
        scale = fileh5.root.Grid._v_attrs.Lz/(2. * pi)
        pl = semilogy(T,scale*fileh5.root.Analysis.PowerSpectrum.Z[field, :numModes,2:].T)
    
        legend_list = []
        for i in range(len(fileh5.root.Analysis.PowerSpectrum.Z[field, :numModes,0])):
            legend_list.append("n = %i" % i)
        leg = legend(legend_list, loc='lower right', ncol=2)
        leg.draw_frame(0)
        
        
    xlabel("Time Step")
    xlabel("Time")
    if       field == 0 : ylabel("Mode Power $|\\phi|^2$")
    elif  field == 1    : ylabel("Mode Power $|A_\\parallel|^2$")
    elif  field == 2    : ylabel("Mode Power $|B_\\parallel|^2$")
    else                : print "No such field"
    return pl, leg
    # framless




def plotXPropTempDensity(time=0, Z=4, Y=6, fileh5=fileh[0]):

  markers = ['v-c', '^-y', '<-r', '>-m', '*-b', 'd-g', 'p-b', '1-r', '2-m']
  pl = []
  pl_name = []
  D = getDomain(fileh5)
  # we iterate of all species
  for s in range(len(fileh.root.Phasespace.Data[:,0,0,0,0,0,0])):
    species_name = fileh.root.Species.cols.Name[s]
    species_name = "species"
    x = D['X']
    pl.append(plot(x, Density(time, s)[Z,Y,:], markers[2*s]))
    pl_name.append("Density (" + species_name + ")")
    xlabel("Position [x]")
    ylabel("Density")
  
    # set scaling for density
    d = Density(time)[Z,Y,:]
    if((max(d) - min(d)) < 0.2): ylim(min(d)-0.1, max(d)+0.1)
  
    twinx()
    pl.append(plot(x, Temperature(time, s)[Z,Y,:],markers[2*s+1]))
    pl_name.append("Temperature (" + species_name + ")")
    ylabel("Temperature")
    title("Density/Temperature profile at Time Step = " + str(time))

    

  leg = legend(pl, pl_name, loc='lower center')
  leg.draw_frame(0)

def plotInfo(fileh5=fileh[0]):
  
  ##description = fileh5.root.Info._v_attrs.Description
  #info        = fileh5.root.Info._v_attrs.Info
  D = getDomain(fileh5)

  string  = '$L_x = %3.1f,\, L_y = %3.1f,\, L_z = %3.1f,\, L_v = %3.1f,\, L_m = %3.1f$\n' % (D['Lx'], D['Ly'], D['Lz'], D['Lv'], D['Lm']) 
  string += '$N_x = %3.1f,\, N_y = %3.1f,\, N_z = %3.1f,\, N_v = %3.1f,\, N_m = %3.1f$\n' % (D['Nx'], D['Ny'], D['Nz'], D['Nv'], D['Nm']) 

  #species 
  if(D['Ns'] > 1):
    string += 'Species\n'
    for s in range(D['Ns']):
      string +=  D['species'][s][0] + ' Charge %1.3f' % D['species'][s][1] + ' Mass %1.3f\n' % D['species'][s][2] 

  pylab.setp(pylab.gca(), frame_on=False, xticks=(), yticks=())

  pylab.text(0, 0,  string)



def Info():
    D = getDomain(fileh)
    print "\n\n"
    print "HDF5 Output Filename : ", fileh.root.Info._v_attrs.Output
    print "Collisions           : ", fileh.root.Info._v_attrs.Collisions
    print "Geometry             : ", fileh.root.Info._v_attrs.Geometry

    print "Model           : ", fileh.root.Info._v_attrs.Model
    print "Physics         : ", fileh.root.Info._v_attrs.Physics
    print "Solver          : ", fileh.root.Info._v_attrs.Solver
    print "Type            : ", fileh.root.Info._v_attrs.Type
    print "Phi Saves            : ", len(fileh.root.Potential.phi[0,0,0,:])
    print 'Lx = %3.1f, Ly = %3.1f, Lz = %3.1f, Lv = %3.1f, Lm = %3.1f' % (D['Lx'], D['Ly'], D['Lz'], D['Lv'], D['Lm']) 
    print 'Nx = %3.1f, Ny = %3.1f, Nz = %3.1f, Nv = %3.1f, Nm = %3.1f' % (D['Nx'], D['Ny'], D['Nz'], D['Nv'], D['Nm']) 
    
    #print "PSF Saves  : ", len(fileh.root.Potential.phi[0,0,0,:])
    #$print "X Prop     : ", len(fileh.root.Analysis.Xpropertie.E[0,0,0,:])



#Info()
def plotMagneticFieldShearedSlab2(fileh5=fileh[0], Z=4, frame=-1):
   D = getDomain(fileh5)

   X,Y = meshgrid(D['X'],D['Y'])
   data = fileh5.root.Visualization.phi[0,:,:,frame]
   T = getTime(fileh5.root.Visualization.Time)[frame,:]
   if frame == 0: dt = 0.
   else         : dt = (getTime(fileh5.root.Visualization.Time)[frame,1]- getTime(fileh5.root.Visualization.Time)[frame-1,1]) / (getTime(fileh5.root.Visualization.Time)[frame,0]- getTime(fileh5.root.Visualization.Time)[frame-1,0])
   contourf(X,Y,data, 100)
   colorbar()



def plotMagneticFieldShearedSlab(fileh5=fileh[0], Z=4, T=-1):

  #def tildePsi(pol_k,x) : return x**2/2.+cos(pol_k*x)

    
  D = getDomain(fileh5)
  
  shat = fileh5.root.Constants._v_attrs.shat[0]
  def addShear(X,Y) : return shat*(X/D['Lx'])**2/2.
    
  X,Y = meshgrid(D['X'],D['Y'])

  Ap = addShear(X,Y) - fileh5.root.Potential.Ap[Z,:,:,T]

  dA_dy, dA_dx = gradient(Ap)
    #B_x, B_y = Ap(X,Y)
    
  streamlines.Streamlines(D['X'],D['Y'], dA_dy, -dA_dx, maxLen=1000., spacing = 2).plotArrows()
    
  xlim((min(D['X']),max(D['X'])))
  ylim((min(D['Y']),max(D['Y'])))
  show()


"""

def plotStreamLines(ntime=-1, nz=4, fileh5=fileh[0]):

    Ap = fileh5.root.Potential.Ap[nz,:,:,ntime]
    dA_dy, dA_dx = gradient(Ap)
  
    D = getDomain(fileh5)
    rect = [0.08, 0.1, 0.75, 0.8 ] 
    x = linspace(0.0, D['Lx'], D['Nx'])
    y = linspace(0.0, D['Ly'], D['Ny'])
    conax = axes(rect, axisbg='w')

    streamplot.streamplot(x,y, dA_dy, -dA_dx)

    xlim((x[0], y[-1]))
    ylim((x[0], y[-1]))
"""




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



def plotVisualize(fileh5=fileh[0], Z=0, frame=-1, text="", interpolate=False):
   D = getDomain(fileh5)

   X,Y = meshgrid(D['X'],D['Y'])
   data = fileh5.root.Visualization.phi[0,:,:,frame].T
   T = getTime(fileh5.root.Visualization.Time)[frame,:]
   clf()
   contourf(X,Y,data, 100)
   colorbar()
   xlim((X.min(), X.max()))
   ylim((Y.min(), Y.max()))
   if printTitle == True : title("T : " + str(T[1]) + " : " + text)




def getRealFromXky(fileh5, data, modes=[], interpolate=False):

   D = getDomain(fileh5)
   # Fourier back transform c2r, renormaliza, not it stored in fftw shape
   Nky = len(data[:,0])
   Nx  = len(D['X'])
   if(interpolate == True) :
       
       new_Ny = max(512, 2*Nky)

       ex_data = zeros((new_Ny, len(D['X'])))

       # note we have no nyquest frequency
       #bring to fftw shape 
       #ex_data[0:Ny/2+1,:]   = data[0:Ny/2+1,:]
       # copy negative frequencyies
       #ex_data[-Ny/2+1:,:] = data[-Ny/2+1:,:]

       #data = ex_data
       #Ny   = new_Ny
  
   
   Y = linspace(D['Y'].min(), D['Y'].max(), D['Ny'])
   X,Y = meshgrid(D['X'], Y)

   if modes != [] :
     for m in range(Ny/2+1):
        # remove unwanted frequencies
        if m not in modes :
          data[m,:] = 0.
          # and negative freq.
          # index :
          #print "Index : ", m , "   Ny - m : ", Ny - m
          #if (m != 0) and (m != Ny/2): data[Ny-m,:] = 0.
   print "Shape data : ", shape(data)
   data = irfft(data, axis=0)

   print "New shape : " , shape(data.T), " Nx : " , len(X) , " Ny : " , len(Y)

   return X, Y, data





def plotVisualizeFFT(fileh5=fileh[0], data="Phi", Z=0, species=0, frame=-1, text="", norm=True, interpolate=False, clearPlot=True, interpolation='bilinear', modes=[], printTitle=True):
   D = getDomain(fileh5)
  
   if   data == "Phi" : 
                        data = fileh5.root.Visualization.Phi[Z,:,:,frame]
                        T = getTime(fileh5.root.Visualization.Time)[frame,:]
                        cmap = cm.jet
   elif data == "Tp"  : 
                        data = fileh5.root.Moments.Temperature_v[Z,:,:,species,frame]
                        T = getTime(fileh5.root.Moments.Time)[frame,:]
                        cmap = cm.gist_heat
   elif data == "HeatFlux"  : 
                        data = fileh5.root.Moments.HeatFlux[Z,:,:,species,frame]
                        T = getTime(fileh5.root.Moments.Time)[frame,:]
                        cmap = cm.jet
                        #YlOrRd
   else : print "No such Data ", data
                      

   X, Y, data = getRealFromXky(fileh5, data, modes, interpolate)
   print "Shape of Data : " , shape(Y), " " , shape(data)
   if norm == True:
            v_minmax =  max(abs(data.min()), abs(data.max()))
            norm = normalize(vmin = -v_minmax, vmax = v_minmax)
   else : 
            norm = normalize(vmin = data.min(), vmax = data.max())
   if clearPlot == True  : clf()
   contourf(X,Y,real(data), 100, cmap=cmap, norm=norm, interpolation=interpolation)
   colorbar(norm=norm)
   xlim((X.min(), X.max()))
   ylim((Y.min(), Y.max()))
   xlabel("Radial Direction")
   ylabel("Poloidal Direction")
   if printTitle == True : title("TimeStep : %i   Time : %.3f " % (T[0], T[1]))



def getData(Var, fileh5, Z=0, frame=-1, species=0):

    if   Var == "Phi" : 
                        data = fileh5.root.Visualization.Phi[Z,:,:,frame]
                        T = getTime(fileh5.root.Visualization.Time)[frame,:]
                        cmap = cm.jet
    elif Var == "Tp"  : 
                        data = fileh5.root.Moments.Temperature_v[Z,:,:,species,frame]
                        T = getTime(fileh5.root.Moments.Time)[frame,:]
                        cmap = cm.gist_heat
    elif Var == "HeatFlux"  : 
                        data = fileh5.root.Moments.HeatFlux[Z,:,:,species,frame]
                        T = getTime(fileh5.root.Moments.Time)[frame,:]
                        #cmap = cm.YlOrRd
                        cmap = cm.jet
    else : print "No such Data ", data


    return data

def plotCorrelateValues(A="Phi",B="Tp",frame=-1, Z=0, species=1, fileh5=fileh[0]):
    D = getDomain(fileh5)

    X, Y, A = getRealFromXky(fileh5, getData(A, fileh5, Z, frame))
    X, Y, B = getRealFromXky(fileh5, getData(B, fileh5, Z, frame))

    #B = getData(B, fileh5, Z, frame)
    #    A = fileh5.root.Potential.phi[:,:,:,timeStep]
    #    B = fileh5.root.Potential.phi[:,:,:,timeStep]


    # first correlate then average
    clf()
    Corr = scipy.signal.correlate(A, B)
    #Corr = (ifftn(fftn(A)*ifftn(B))).real
  
    #Corr = sum(Corr, axis=0)
    #    x_n  = list(-(D['kx']))
    #    x_n.reverse()
    #    x = x_n + list((D['kx'][:]))
    #    y = [ 0 ] + list((D['ky']))
    #    print "Shape Ap : ", shape(Corr), " x : " , len(x), "   y : " , len(y)
    contourf(Corr[:,:].T,100, cmap=cm.hot)#, locator=ticker.LogLocator())
    colorbar()

def plotCrossCorrelateValues(A="Phi",B="Tp",frame=-1, Z=0, species=1, fileh5=fileh[0]):
    D = getDomain(fileh5)

    X, Y, F = getRealFromXky(fileh5, getData(B, fileh5, Z, frame))
    A = getData(A, fileh5, Z, frame)
    B = getData(B, fileh5, Z, frame)
    
    #    A = fileh5.root.Potential.phi[:,:,:,timeStep]
    #    B = fileh5.root.Potential.phi[:,:,:,timeStep]


    # first correlate then average
    C = []
    clf()
    print "Real : shape : ", shape(B)
    #for nky in range(len(A[:,0])): C.append(abs(ifftn(fft(A[nky,:])*ifftn(B[nky,:])).imag))
    for nky in range(len(A[:,0])): C.append(abs(scipy.signal.correlate(A[nky,:], B[nky,:].imag)))
    Corr = array(C)
  
    #Corr = sum(Corr, axis=0)
    #    x_n  = list(-(D['kx']))
    #    x_n.reverse()
    #    x = x_n + list((D['kx'][:]))
    #    y = [ 0 ] + list((D['ky']))
    #    print "Shape Ap : ", shape(Corr), " x : " , len(x), "   y : " , len(y)
    X = linspace(-pi, pi, len(Corr[0,:]))
    Y = (D['ky'])

    print " Y : ", len(Y)
    print " X : ", len(X)
    print " Shap : " , shape(Corr[:len(Y),:])

    ax = subplot(111)

    ax.contourf(X,Y,Corr[:len(Y),:],100, cmap=cm.hot)
    #colorbar()
    ax.set_yscale("log") 
    # set ticks and tick labels
    ax.set_xlim((-pi, pi))
    ax.set_xticks([-pi,0,pi])
    pichr = unichr(0x03C0)
    ax.set_xticklabels(['$\\pi$','0', '$\\pi$'])
    ax.plot(linspace(0., 0., 101), linspace(Y.min(), Y.max(), 101), 'r-')
    ylim((Y.min(), Y.max()))
    xlabel("Phase")
    ylabel("$k_y$")
    #xlim((X.min(), X.max()))
    # plot zero phase line

