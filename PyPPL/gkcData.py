

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


def getData(VarID, fileh5, Z=0, frame=-1, species=0):

    if   Var == "2DPhi" : 
                        data = fileh5.root.Visualization.Phi[Z,:,:,frame]
                        T = getTime(fileh5.root.Visualization.Time)[frame,:]
    elif Var == "2DTp"  : 
                        data = fileh5.root.Moments.Temperature_v[Z,:,:,species,frame]
                        T = getTime(fileh5.root.Moments.Time)[frame,:]
    elif Var == "2DHeatFlux"  : 
                        data = fileh5.root.Moments.HeatFlux[Z,:,:,species,frame]
                        T = getTime(fileh5.root.Moments.Time)[frame,:]
    
    ## One dimensional scalar Data
    """
    elif Var == "1DPhiTMode" : 
            T   = getTime(fileh5.root.Analysis.PowerSpectrum.Time)[:,1]
           data =   fileh5.root.Analysis.PowerSpectrum.Y[field,m,:]
    else : print "No such Data ", data
    """
    D = getDomain(fileh5)

    return D, T, data



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
     for m in range(Nky/2+1):
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


