import gkcData
import gkcStyle
import pylab
import numpy as np
import scipy.optimize


def fitDampedOscillationOptmize(T, Y):

    fitfunc = lambda p, t: p[0]*np.sin(p[1]*t+p[2]) * np.exp(+p[3]*t) # Target function
    errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
    
    p, success = scipy.optimize.leastsq(errfunc, [1.0e-5, -1.0, 0.0, 0.00001], args=(T, Y))

    return abs(p[1]) + 1.j* p[3]

def fitDampedOscillationFourier(T, Y):

       # get frequency
       hatY = np.fft.rfft(Y)
       freq = np.fft.fftfreq(len(Y), d=T[21]-T[20])

       idx = hatY.argmax()
       w0  =  freq[idx] / (2. * np.pi)
    
       #fftfreq = np.fft.fftfreq(len(time_series), d = (T[-10]-T[-11])) 
       #abs_freq =  2.*np.pi*fftfreq[m]


       print "freq " , freq, " idx : " , idx, " w0 : ", w0
       #we could caluclated now gamme but this is not accurate enough
       # so we fit the decline

       # Try simple smoothening
       Y = scipy.ndimage.gaussian_filter(Y,0.5)
       a = abs(Y)
       a_idx = np.r_[True, a[1:] > a[:-1]] & np.r_[a[:-1] > a[1:], True]
       
       # we ignore first and last value due to errors from initial condition and stop
       a_idx = a_idx[1:-2]

       fitfunc = lambda p, x: p[0]*x + p[1] 
       errfunc = lambda p, x, y: fitfunc(p, x) - y
       
       if len(Y[a_idx]) > 1 : p, success = scipy.optimize.leastsq(errfunc, [-1.0e-5, 0.0], args=(T[a_idx], np.log(a[a_idx])))
       else : p = [1.e-6, 1.e-1]
       return w0 + 1.j * p[0], np.exp(p[1])


def getGrowthrate(T, D, start, stop):

  print "Fitting growthrates from T : ", T[start], " - ", T[stop]
  

  def fitExpGrowthOptimize(T,Y):

       fitfunc = lambda p, x: p[0]*x + p[1] 
       errfunc = lambda p, x, y: fitfunc(p, x) - y
       
       #p, success = scipy.optimize.leastsq(errfunc, [0.001, min(np.log(abs(Y)))], args=(T, np.log(abs(Y))))
       p, success = scipy.optimize.leastsq(errfunc, [0.05, 0.], args=(T, np.log(abs(Y))))
       print p, success
       return p[0]
  
  if D.ndim == 1:
    return fitExpGrowthOptimize(T[start:stop],D[start:stop])
  elif D.ndim == 2:      
    growth = []
    N      = len(D[:,0]) 
    for n in range(N): growth.append(fitExpGrowthOptimize(T[start:stop],D[n,start:stop]))
    return np.array(growth)
  else : raise TypeError("Dimension of Array should be either one or two")



def getFrequency(T, D, start=0, stop=-1):
  import scipy.ndimage
  
  freq_list   = []
  
  if   np.ndim(D) == 1 :  N = 1
  elif np.ndim(D) == 2 :  N = len(D[:,0]) 
  else : raise TypeError("Dimension of Array should be either one or two")
    
  # Note We assume constant time-steps !
  for n in range(N): 
    time_series = 0 
    if   N == 1 : time_series = D[start:stop]
    elif N >  1 : time_series = D[n,start:stop]
    else : raise TypeError("Dimension of Array should be either one or two")
    FS = np.fft.rfft(time_series) #np.sin(time_series))
    # Get Maximum Frequency (ignore DC component)
    m     = np.argmax(abs(FS[1:]))+1
    fftfreq = np.fft.fftfreq(len(time_series), d = (T[-10]-T[-11])) 
   
    # crap should use rfftfreq and skip the two
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

  if   N == 1 : return freq_list[0]
  elif N >  1 : return np.array(freq_list)







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
    
    doCLF    = kwargs.pop('doCLF', True)
    dir      = kwargs.pop('dir', 'Y')
    modes    = kwargs.pop('modes' , range(D['Nky']))
    field    = kwargs.pop('field', 'phi')  
    n_offset = kwargs.pop('offset', 2)  
    label    = kwargs.pop('label', 'ky')  
    leg_loc  = kwargs.pop('loc', 'best')  
    start    = kwargs.pop('start', 1)  
    stop     = kwargs.pop('stop', -1)  
    m        = kwargs.pop('m', 0)  
    useLog   = kwargs.pop('useLog', True ) 
   
    if useLog == True : pf = pylab.semilogx
    else              : pf = pylab.plot

    if doCLF == True : pylab.clf()
    
    T = gkcData.getTime(fileh5.root.Analysis.PowerSpectrum.Time)[:,1]
        
    if   field == 'phi' : n_field = 0
    elif field == 'A'   : n_field = 1
    elif field == 'B'   : n_field = 2
    else : raise TypeError("Wrong argument for field : " + str(field))

    if(dir == 'X'):
      pl = pf(T, fileh5.root.Analysis.PowerSpectrum.X[n_field,:numModes,2:].T)
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
        growthrates = getGrowthrate(T,power, start,stop)
        pl = pf(D['ky'], growthrates, "s" + markline, label='$\\gamma$')
      if which=='r' or which=='b':
        shift = fileh5.root.Analysis.PhaseShift.Y   [n_field, :,:]
        frequency   = getFrequency(T,shift, start,stop)
        print "ky : ", np.shape(D['ky']), " freq : ", np.shape(frequency)
        pl = pf(D['ky'], frequency, "v" + markline, label='$\\omega_r$')
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


def plotEigenvalues(fileh5, **kwargs):
    """
        Plot the eigenvalues of the phase space function


        Optional keyword arguments:

        Keyword           Description
        ===============   ==============================================
         *dir*             Direction 'X' (for radial) or 'Y' for poloidal
         *offset*          Offset due to zeroset to 2 .

    """
    import gkcStyle


    D = gkcData.getDomain(fileh5)

    try:
        eigv = fileh5.root.Eigenvalue.EigenValues.cols.Eigenvalue[:] 
    except:
       print "Problems openning eigenvalues. Not included ?"
       return


    eigv_r = np.real(eigv)
    eigv_i = np.imag(eigv)

    pylab.plot(eigv_r, eigv_i, '.')


    gkcStyle.plotZeroLine(1.1*min(eigv_r), 1.1*max(eigv_r), direction='horizontal', color="#666666", lw=0.8)
    gkcStyle.plotZeroLine(1.1*min(eigv_i), 1.1*max(eigv_i), direction='vertical'  , color="#666666", lw=0.8)

    pylab.xlim((1.1*min(eigv_r), 1.1*max(eigv_r)))
    pylab.ylim((1.1*min(eigv_i), 1.1*max(eigv_i)))


    pylab.xlabel("$\\omega_r$")
    pylab.ylabel("$\\omega_i$")


def plotEigenfunctions(fileh5, mode=1, n_max=5, **kwargs):
    """
        Plot the n largerst eigenfunctions

        Optional keyword arguments:

        Keyword           Description
        ===============   ==============================================
         *dir*             Direction 'X' (for radial) or 'Y' for poloidal
         *offset*          Offset due to zeroset to 2 .

    """
    import gkcStyle
    import gkcAnalysis


    D = gkcData.getDomain(fileh5)

    try:
        eigv = fileh5.root.Eigenvalue.EigenValues.cols.Eigenvalue[:] 
    except:
       print "Problems openning eigenvalues. Not included ?"
       return


    # Sort eigenvalues (note automatically sorts real values first)
    eigv     = fileh5.root.Eigenvalue.EigenValues.cols.Eigenvalue[:]
    # absteigend sort
    idx_sort = np.argsort(eigv)[::-1]
            
    pylab.subplot(121)
    plotEigenvalues(fileh5, **kwargs)

    # Plot points
    for n in range(n_max):
            w = eigv[idx_sort[n]]
            pylab.plot(np.real(w), np.imag(w), 'o', markersize=7.5, color=gkcStyle.markers_C[n])

    pylab.subplot(122)

    labels = []
    for n in range(n_max):
            w = eigv[idx_sort[n]]
            gkcAnalysis.plotModeStructure(fileh5=fileh5, fied='phi', mode=mode, frame=idx_sort[n], part='a', m=n)
            labels.append("%.3f+%.3f i" % (np.real(w), np.imag(w)))

    pylab.legend(labels).draw_frame(0)

def getFrequencyGrowthrates(fileh5, start=1, stop=-1, q=0):
          
   T = gkcData.getTime(fileh5.root.Analysis.PowerSpectrum.Time)[:,1]

   power = fileh5.root.Analysis.PowerSpectrum.Y[q, :,:]
   phase = fileh5.root.Analysis.PhaseShift.Y   [q, :,:]
        
   frequency   = getFrequency (T,phase, start,stop)
   growthrates = getGrowthrate(T,power, start, stop)

   return np.array(frequency) + 1.j * np.array(growthrates)


