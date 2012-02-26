from scipy.signal import *
from pylab import *
from scipy import *
from scipy.special.orthogonal import *
from PyPPL import *
from scipy.optimize import *
from numpy import *
from scipy.signal import *
from pylab import *

import mpmath as mp

import HeliosPlot




def getComplexZero(Func, S, Init=complex(0.01,-0.01), ztol=1.e-6, ftol=1.e-6, dx=0.05, maxfev=9999):
   

    # need to hack a bit around as input parameters are different

    def NewtonFunc(vec_z):
        w = Func(complex(vec_z[0],vec_z[1]), S)
        return [real(w), imag(w)]
   
    w     =  zermuller(Func, S, Init, ztol=ztol, ftol=ftol, dx=dx)[0]
    w_vec = fsolve(NewtonFunc, [real(w), imag(w)], maxfev=maxfev)
    w = complex(w_vec[0], w_vec[1])
    
    return w


def fitDampedOscillationOptmize(T, Y):

    fitfunc = lambda p, t: p[0]*sin(p[1]*t+p[2]) * exp(+p[3]*t) # Target function
    errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
    
    p, success = optimize.leastsq(errfunc, [1.0e-5, -1.0, 0.0, 0.00001], args=(T, Y))

    return (abs(p[1]), p[3])



def fitDampedOscillationFourier(T, Y):

       # get frequency
       print " Y : " , Y
       hatY = rfft(Y[20:])
       print " hatY : " , hatY
       freq = fftfreq(len(Y[20:]), d=T[21]-T[20])[:len(T[20:])/2+1]

       idx = hatY.argmax()
       w0  =  freq[idx]


       print "freq " , freq, " idx : " , idx, " w0 : ", w0
       #we could caluclated now gamme but this is not accurate enough
       # so we fit the decline

       a = abs(Y)
       a_idx = numpy.r_[True, a[1:] > a[:-1]] & numpy.r_[a[:-1] > a[1:], True]

       fitfunc = lambda p, x: p[0]*x + p[1] 
       errfunc = lambda p, x, y: fitfunc(p, x) - y
       
       if len(log(abs(Y[a_idx]))) > 1: p, success = optimize.leastsq(errfunc, [-1.0e-5, 0.0001], args=(T[a_idx], log(abs(Y[a_idx]))))
       else : p = [1.e-6]
       return (w0, p[0])



def fitExpGrowthOptimize(T,Y,offset=0):

       fitfunc = lambda p, x: p[0]*x + p[1] 
       errfunc = lambda p, x, y: fitfunc(p, x) - y
       
       p, success = optimize.leastsq(errfunc, [-1.0e-5, 0.0001], args=(T[offset:], log(abs(Y[offset:]))))
       return p[0]




def getGrowth2(fileh5, pos=(-10,-1), dir='Y'):
    
  D = getDomain(fileh5)

  dPhi = log(fileh5.root.Analysis.PowerSpectrum.Y[0,1:,pos[1]])-log(fileh[0].root.Analysis.PowerSpectrum.Y[0,1:,pos[0]])
  dT   = fileh5.root.Analysis.PowerSpectrum.Time[pos[1]][1]-fileh5.root.Analysis.PowerSpectrum.Time[pos[0]][1]

  print "Fitting from T = ", fileh5.root.Analysis.PowerSpectrum.Time[pos[0]][1], \
                 " to T = " , fileh5.root.Analysis.PowerSpectrum.Time[pos[1]][1]

  return (D['ky'], dPhi/dT)

def getGrowth(fileh5, pos=(-10,-1), dir='Y'):
  #import FitFunction  
  if   dir == 'Y': Phi = fileh5.root.Analysis.PowerSpectrum.Y[0,1:,pos[0]:pos[1]]
  elif dir == 'X': Phi = fileh5.root.Analysis.PowerSpectrum.X[0,1:,pos[0]:pos[1]]
  else : print "No Suck direction"
  #T   = HeliosPlot.getTime(fileh5.root.Analysis.PowerSpectrum.Time[pos[0]:pos[1]])[:,1]
  T   = HeliosPlot.getTime(fileh5.root.Analysis.PowerSpectrum.Time[pos[0]:pos[1]])[:,1]

  D = HeliosPlot.getDomain(fileh5)
  
  growth = []
  for nky in range(len(D['ky'])): growth.append(fitExpGrowthOptimize(T,Phi[nky,:]))

  print "Fitting from T = ", fileh5.root.Analysis.PowerSpectrum.Time[pos[0]][1], \
                 " to T = " , fileh5.root.Analysis.PowerSpectrum.Time[pos[1]][1]

  return (D['ky'], growth)

def getGrowthOld(fileh5, pos=(-10,-1), dir='Y'):
  #import FitFunction  

  Phi = fileh5.root.Analysis.PowerSpectrum.Y[0,1:,pos[0]:pos[1]]
  #T   = HeliosPlot.getTime(fileh5.root.Analysis.PowerSpectrum.Time[pos[0]:pos[1]])[:,1]
  T   = HeliosPlot.getTime(fileh5.root.Analysis.PowerSpectrum.Timing[pos[0]:pos[1]])[:,1]

  D = HeliosPlot.getDomain(fileh5)
  
  growth = []
  for nky in range(len(D['ky'])): growth.append(fitExpGrowthOptimize(T,Phi[nky,:]))

  print "Fitting from T = ", fileh5.root.Analysis.PowerSpectrum.Timing[pos[0]][1], \
                 " to T = " , fileh5.root.Analysis.PowerSpectrum.Timing[pos[1]][1]


  return (D['ky'], growth)


def getGrowthEnergy(fileh5, pos=(-1,-10)):
    
  D = getDomain(fileh5)

  dPhi = log(fileh5.root.Analysis.scalarValues.cols.phiEnergy[pos[0]]) - log(fileh5.root.Analysis.scalarValues.cols.phiEnergy[pos[1]])
  dT   = fileh5.root.Analysis.scalarValues.cols.Time[pos[0]] - fileh5.root.Analysis.scalarValues.cols.Time[pos[1]]
     
  # factor 1/2 because we are looking for phi=exp[gamma*t] but calculating from phi^2
  return 0.5*dPhi/dT



