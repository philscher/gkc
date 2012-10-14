import numpy as np

def getTime(timing):
  timing_new = []
  for t in timing:
    timing_new.append((t[0], t[1]))
  return np.array(timing_new)


def getDomain(fileh5):
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
           'Tscale' : np.sqrt(1837.) / (2. *np.sqrt(2.) * 1.e-3), \
           'X'     : fileh5.root.Grid._v_attrs.X, \
           'Z'     : fileh5.root.Grid._v_attrs.Z, \
           'V'     : fileh5.root.Grid._v_attrs.V, \
           'M'     : fileh5.root.Grid._v_attrs.M, \
           'kx'   : np.arange(0, Nx/2+1) * 2 * np.pi / Lx,\
           'ky'   : np.arange(0, Ny    ) * 2 * np.pi / Ly,\
           'kp'   : np.arange(0, Nz/2+1) * 2 * np.pi / Lz,\
           'species': species,\
           'Debye2' : fileh5.root.Plasma._v_attrs.Debye2[0],\
           }
  D['Y'] = np.linspace(0., D['Ly'], D['Ny'])
  D['TScale'] =  D['Lv']/D['Lz']

  return D


def getData(Var, fileh5, Z=0, frame=-1, species=0):

    if   Var == "2DPhi" : 
                        data = fileh5.root.Visualization.Phi[Z,:,:,frame]
                        T = getTime(fileh5.root.Visualization.Time)[frame,:]
    elif Var == "2DAp"  : 
                        data = fileh5.root.Visualization.Ap[Z,:,:,frame]
                        T    = getTime(fileh5.root.Visualization.Time)[frame,:]
    elif Var == "2DTp"  : 
                        data = fileh5.root.Moments.Temperature_v[Z,:,:,species,frame]
                        T = getTime(fileh5.root.Moments.Time)[frame,:]
    elif Var == "2Dn"  : 
                        data = fileh5.root.Moments.Density[Z,:,:,species,frame]
                        T = getTime(fileh5.root.Moments.Time)[frame,:]
    elif Var == "2DHeatFlux"  : 
                        data = fileh5.root.Moments.HeatFlux[Z,:,:,species,frame]
                        T = getTime(fileh5.root.Moments.Time)[frame,:]
    else                : raise TypeError("No such var")
    
    ## One dimensional scalar Data
    """
    elif Var == "1DPhiTMode" : 
            T   = getTime(fileh5.root.Analysis.PowerSpectrum.Time)[:,1]
           data =   fileh5.root.Analysis.PowerSpectrum.Y[field,m,:]
    else : print "No such Data ", data
    """
    D = getDomain(fileh5)

    return D, T, data



def getRealFromXky(fileh5, data, modes=[], min_modes=129, interpolate=False):

   D = getDomain(fileh5)
   
   # Fourier back transform c2r, renormaliza, not it stored in fftw shape
   Nky = len(data[:,0])
   Nx  = len(data[0,:])
   
   if modes != [] :
     for m in range(Nky):
        # remove unwanted frequencies
        if m not in modes :
          data[m,:] = 0.
   
   Ny_mod = 2 * max(len(data[:,0]), min_modes)-2
   data = np.fft.irfft(data, n=Ny_mod, axis=0)
   
   Y = np.linspace(D['Y'].min(), D['Y'].max(), Ny_mod)

   return D['X'], Y, data


def plotInfo(fileh5):
  
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


