import gkcData
import gkcStyle

import pylab
import numpy as np
import scipy.optimize


def fitDampedOscillationOptmize(T, Y):

    fitfunc = lambda p, t: p[0]*sin(p[1]*t+p[2]) * exp(+p[3]*t) # Target function
    errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
    
    p, success = scipy.optimize.leastsq(errfunc, [1.0e-5, -1.0, 0.0, 0.00001], args=(T, Y))

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
       
       if len(log(abs(Y[a_idx]))) > 1: p, success = scipy.optimize.leastsq(errfunc, [-1.0e-5, 0.0001], args=(T[a_idx], log(abs(Y[a_idx]))))
       else : p = [1.e-6]
       return (w0, p[0])


def getGrowthrate(T, D, start, stop, dir='Y'):

  growth = []

  print "Fitting from T : ", T[start], " - ", T[stop]

  def fitExpGrowthOptimize(T,Y):

       fitfunc = lambda p, x: p[0]*x + p[1] 
       errfunc = lambda p, x, y: fitfunc(p, x) - y
       
       p, success = scipy.optimize.leastsq(errfunc, [-1.0e-5, 0.0001], args=(T, np.log(abs(Y))))
       return p[0]
  
  
  if D.ndim == 1:
        growth.append(fitExpGrowthOptimize(T[start:stop],D[start:stop]))
  elif D.ndim == 2:      
    N      = len(D[:,0]) 
    for n in range(N): growth.append(fitExpGrowthOptimize(T[start:stop],D[n,start:stop]))
  else : TypeError("Dimension of Array should be either one or two")

  print "Getting Mode Growth from T = ", T[start], " to T = " , T[stop]

  return np.array(growth)

def getFrequency(T, D, start, stop, dir='Y'):
  import scipy.ndimage
  
  freq_list   = []
  N      = len(D[:,0]) 
    
  print "Fitting from T : ", T[start], " - ", T[stop]
 
  # Note We assume constant time-steps !
  for n in range(N): 
     
    time_series = D[n,start:stop]
    FS = np.fft.rfft(time_series) #np.sin(time_series))
    # Get Maximum Frequency
    m     = np.argmax(abs(FS))
    fftfreq = np.fft.fftfreq(len(time_series), d = (T[-10]-T[-11])) 
    
    abs_freq =  2.*np.pi*fftfreq[m]

    # Needs sqrt(2.) from velocity normalization
    
    # Get sign of frequency by taking gradient of phase shift (how to deal with jump?)
    time_series = scipy.ndimage.gaussian_filter(time_series, 0.01)
    grad = np.gradient(time_series, T[-10]-T[-11])
    # remove jump values
    np.putmask(grad, abs(grad) > 1.05*abs_freq, 0.)
    np.putmask(grad, abs(grad) < 0.95*abs_freq, 0.)
    sig = -np.sign(sum(grad))

    freq_list.append(sig * abs_freq)

  print "Getting Frequency from T = ", T[start], " to T = " , T[stop]

  return np.array(freq_list)







def getGrowthEnergy(fileh5, pos=(-1,-10)):
    
  D = getDomain(fileh5)

  dPhi = np.log(fileh5.root.Analysis.scalarValues.cols.phiEnergy[pos[0]]) - np.log(fileh5.root.Analysis.scalarValues.cols.phiEnergy[pos[1]])
  dT   = fileh5.root.Analysis.scalarValues.cols.Time[pos[0]] - fileh5.root.Analysis.scalarValues.cols.Time[pos[1]]
     
  # factor 1/2 because we are looking for phi=exp[gamma*t] but calculating from phi^2
  return 0.5*dPhi/dT


def plotFrequencyGrowthrates(fileh5, which='b', markline="-", **kwargs):
    """
        Plots time evolution of mode power.


        Optional keyword arguments:

        Keyword           Description
        ===============   ==============================================
         *dir*             Direction 'X' (for radial) or 'Y' for poloidal
         *modes*           List of modes (default plot modes). 
                           e.g. modes = [1,4,5]         - to plot all modes
                                modes = range(Nky)[::2] - to plot every second mode
         *field*           'phi' electric potential
                           'A' parallel magnetic vector potential
                           'B' parallel magnetic field
         *dir*             Direction 'X' (for radial) or 'Y' for poloidal
         *doCFL*           clear previous figure
         *label*           'ky' or 'm'
         *offset*          Offset due to zeroset to 2 .

    """
    D = gkcData.getDomain(fileh5)
    
    doCFL    = kwargs.pop('doCFL', True)
    dir      = kwargs.pop('dir', 'Y')
    modes    = kwargs.pop('modes' , range(D['Nky']))
    field    = kwargs.pop('field', 'phi')  
    n_offset = kwargs.pop('offset', 2)  
    label    = kwargs.pop('label', 'ky')  
    leg_loc  = kwargs.pop('loc', 'best')  
    start    = kwargs.pop('start', 1)  
    stop     = kwargs.pop('stop', -1)  
    

    if doCFL == True : pylab.clf()
    
    T = gkcData.getTime(fileh5.root.Analysis.PowerSpectrum.Time)[:,1]
        
    print "Fitting from T : ", T[start], " - ", T[stop]
    
    if   field == 'phi' : n_field = 0
    elif field == 'A'   : n_field = 1
    elif field == 'B'   : n_field = 2
    else : raise TypeError("Wrong argument for field : " + str(field))

    if(dir == 'X'):
      pl = semilogy(T, fileh5.root.Analysis.PowerSpectrum.X[n_field,:numModes,2:].T)
      legend_list = []
      for i in range(len(fileh5.root.Analysis.PowerSpectrum.X[n_field, :numModes,0])):
        legend_list.append("kx = %i" % i)
      leg = pylab.legend(legend_list, loc='lower right', ncol=2)
      leg.draw_frame(0)
    
    elif(dir == 'Y'):

      scale = fileh5.root.Grid._v_attrs.Ly/(2. * np.pi)
      
      legend_list = []
      if   which=='i' or which=='b':
        power = fileh5.root.Analysis.PowerSpectrum.Y[n_field, :,:]
        growthrates = getGrowthrate(T,power, start,stop, dir='Y')
        pl = pylab.semilogx(D['ky'], growthrates, "s" + markline, label='$\\gamma$')
      if which=='r' or which=='b':
        shift = fileh5.root.Analysis.PhaseShift.Y   [n_field, :,:]
        frequency   = getFrequency(T,shift, start,stop, dir='Y')
        pl = pylab.semilogx(D['ky'], frequency, "v" + markline, label='$\\omega_r$')
      #if which !='r' or which != 'i' or which !='b':
      #      raise TypeError("Wrong argument for which (r/i/b) : " + str(dir))
    
      #pylab.twinx()
      pylab.xlim((0.8*min(D['ky']), 1.2*max(D['ky'])))

    
      
    else : raise TypeError("Wrong argument for dir : " + str(dir))
     
    pylab.xlabel("$k_y$")
   
    leg = pylab.legend(loc=leg_loc, ncol=1, mode="expand").draw_frame(0)

    gkcStyle.plotZeroLine(0.8*min(D['ky']), 1.2*max(D['ky']))

    pylab.ylabel("Growthrate $\\gamma(k_y)$ / Frequency $\\omega_r(k_y)$")
    
    #return pl, leg

