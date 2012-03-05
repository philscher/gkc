import gkcData


def plotTimeEvolutionModePower(fileh5=fileh[0], **kwargs):
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
    D = getDomain(fileh5)
    
    doCFL    = kwargs.pop('doCFL', True)
    dir      = kwargs.pop('dir', 'Y')
    modes    = kwargs.pop('modes' , range(D['Nky']))
    field    = kwargs.pop('field', 'phi')  
    n_offset = kwargs.pop('offset', 2)  
    label    = kwargs.pop('label', 'ky')  
    leg_loc  = kwargs.pop('loc', 'best')  
    
    if doCFL == True : pylab.clf()
    
    T = getTime(fileh5.root.Analysis.PowerSpectrum.Time)[2:,1]
    
    if   field == 'phi' : n_field = 0
    elif field == 'A'   : n_field = 1
    elif field == 'B'   : n_field = 2
    else : raise TypeError("Wrong argument for field : " + str(field))

    if(dir == 'X'):
      pl = semilogy(T, fileh5.root.Analysis.PowerSpectrum.X[n_field,:numModes,2:].T)
      legend_list = []
      for i in range(len(fileh5.root.Analysis.PowerSpectrum.X[n_field, :numModes,0])):
        legend_list.append("kx = %i" % i)
      leg = legend(legend_list, loc='lower right', ncol=2)
      leg.draw_frame(0)
    
    elif(dir == 'Y'):

      scale = fileh5.root.Grid._v_attrs.Ly/(2. * pi)
      
      legend_list = []
      
      for m in modes:
            pl = pylab.semilogy(T, (scale*fileh5.root.Analysis.PowerSpectrum.Y[n_field, m,n_offset:]).T)
            if  (label == 'm' ) : legend_list.append("m = %i" % m)
            elif(label == 'ky') : legend_list.append("ky = %.1f" % (m / scale)) 
            else                : print "Name Error"
    
      leg = legend(legend_list, loc=leg_loc, ncol=4, mode="expand").draw_frame(0)
    else : raise TypeError("Wrong argument for dir : " + str(dir))
     
    xlabel("Time")
    xlim((0.,max(T)))


    if    field == "phi" : ylabel("Mode Power $|\\phi|^2$")
    elif  field == "A"   : ylabel("Mode Power $|A_\\parallel|^2$")
    elif  field == "B"   : ylabel("Mode Power $|B_\\parallel|^2$")
    else : raise TypeError("Wrong argument for field : " + str(field))
    
    return pl, leg


def plotInstantGrowthrates(numModes=10, dir='Y', fileh5=fileh[0], order=11, filterType='hanning', field=0):
    """
        Plots instant growthrates of mode power

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
    
    doCFL    = kwargs.pop('doCFL', True)
    dir      = kwargs.pop('dir', 'Y')
    modes    = kwargs.pop('modes' , range(D['Nky']))
    field    = kwargs.pop('field', 'phi')  
    n_offset = kwargs.pop('offset', 2)  
    label    = kwargs.pop('label', 'ky')  
    leg_loc  = kwargs.pop('loc', 'best')  
    

    if(dir == 'X'):

      plot(T, grad)
      legend_list = []
      for i in range(len(fileh5.root.Analysis.PowerSpectrum.X[field,:numModes,0])):
        legend_list.append("kx = %i" % i)
      
      leg = legend(legend_list, loc='best', ncol=4)
      leg.draw_frame(0)
    
    if(dir == 'Y'):
      for m in numModes:
        Var =   fileh5.root.Analysis.PowerSpectrum.Y[field,m,:]
        gamma = array(gradient(log10(Var), T))
        gamma = smoothFilter(gamma, filterType, window_len=order)[0:-order+1]
        plot(T, gamma, label='ky = %.2f' % (D['ky'][m]), color='#'+colors_216[(16*m+1) % 215])
      
      leg = legend(loc='best', ncol=4)
      leg.draw_frame(0)
    
    
    xlabel("Time")
    xlim((0.,max(T)))

    ylabel("Mode Growth $\gamma$")
        

def plotCrossCorrelateValues(A="Phi",B="Tp",frame=-1, Z=0, species=1, fileh5=fileh[0]):
    D = getDomain(fileh5)

    X, Y, F = getRealFromXky(fileh5, gkcData.gkcData.getData(B, fileh5, Z, frame))
    A = gkcData.getData(A, fileh5, Z, frame)
    B = gkcData.getData(B, fileh5, Z, frame)
    
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


################### Ploshed functions


def plotVisualizeFFT(fileh5=fileh[0], data="Phi"):
    """
        Plots Scalar Data (2D) in XY  coordinates


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
    norm     = kwargs.pop('Normalize', True)
    Z        = kwargs.pop('Z', 0)
    modes    = kwargs.pop('modes' , range(D['Nky']))
    doCFL    = kwargs.pop('doCFL' , True)
    interpolation = kwargs.pop('interpolation' , 'bilinear')
    printTitle    = kwargs.pop('printTitle' , True)
    orientation   = kwargs.pop('orientation' , 'vertical')
    frame         = kwargs.pop('frame' , -1)
    
    
  ###
   
   D, T, data = gkcData.getData(Field, fileh5, Z, frame, species)
  
   X, Y, data = getRealFromXky(fileh5, data, modes, interpolate)
   
   HeliosStyle.plotContourWithColorbar(X,Y, data, cmap, norm, interpolation, orientation)
   
   
   xlabel("Radial Direction")
   ylabel("Poloidal Direction")
   if printTitle == True : title("TimeStep : %i   Time : %.3f " % (T[0], T[1]))




def plotScalarDataTimeEvolution(fileh5=fileh[0], var = "phi", log=0, opt="", legend=1):
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
    D = getDomain(fileh5)
    
    doCFL    = kwargs.pop('doCFL', True)
    log      = kwargs.pop('log', True)

    if(log == 1) : 
      plotf = semilogy
      af    = abs
    else         : 
      plotf = plot
      af    = lambda x : x
        
    T = fileh5.root.Analysis.scalarValues.cols.Time[3:]

  for s in range(D['Ns']):


    # Plot Field Energy"
    if var == "phi":
        plotf(T, af(fileh5.root.Analysis.scalarValues.cols.phiEnergy[3:]), opt)
        ylabel("Field Energy")

    if var == "kinEnergy" or var == "KinPhi":
        for s in sp:
            plotf(T, af(fileh5.root.Analysis.scalarValues.cols.KineticEnergy[3:][:,s]) , label="Kinetic Energy (" + fileh5.root.Species.cols.Name[s]+ ")" )
        plotf(T, af(af(fileh5.root.Analysis.scalarValues.cols.phiEnergy[3:])+fak*af(sum(fileh5.root.Analysis.scalarValues.cols.KineticEnergy[3:], axis=1))) , label="Total Energy")
        ylabel("Energy")
    
    if var == "kinEnergy" or var == "KinPhi":
        for s in sp:
            plotf(T, af(fileh5.root.Analysis.scalarValues.cols.KineticEnergy[3:][:,s]) , label="Kinetic Energy (" + fileh5.root.Species.cols.Name[s]+ ")" )
        plotf(T, af(af(fileh5.root.Analysis.scalarValues.cols.phiEnergy[3:])+fak*af(sum(fileh5.root.Analysis.scalarValues.cols.KineticEnergy[3:], axis=1))) , label="Total Energy")
        plotf(T, af(fileh5.root.Analysis.scalarValues.cols.phiEnergy[3:]), opt, label="Field Energy")
        # Plot Toal nergy
        ylabel("Energy")

  if var == "heatFlux":
  
    for s in sp:
      plot(T, fileh5.root.Analysis.scalarValues.cols.HeatFlux[3:][:,s], label=fileh5.root.Species.cols.Name[s]) 
    ylabel("Heat Flux")
 
  if var == "particleNumber":
    for s in sp:
        plotf(fileh5.root.Analysis.scalarValues.cols.Time[3:], abs(fileh5.root.Analysis.scalarValues.cols.ParticleNumber[3:][:,s]) , label=fileh5.root.Species.cols.Name[s])
    if D['Ns'] > 1:
        charge = 0.
        for s in range(D['Ns']):
            charge = charge + fileh5.root.Species.cols.Charge[s]*fileh5.root.Analysis.scalarValues.cols.ParticleNumber[3:][:,s]
        
        plotf(fileh5.root.Analysis.scalarValues.cols.Time[3:], abs(charge) , label="Charge")
  
  if var == "particleFlux":
    for s in sp:
        plotf(T, fileh5.root.Analysis.scalarValues.cols.ParticleFlux[3:][:,s], label=fileh5.root.Species.cols.Name[s]) 
    ylabel("Particle Flux")
  
  leg = legend(loc='lower right', ncol=2)
  leg.draw_frame(0)
  xlabel("Time")


def plotTurbulenceTime(dir='Y', pos=(-2,-1), fileh5=fileh[0], doFit='False', posT=(1,-1)):
  
  D = getDomain(fileh5)
  data = fileh5.root.Analysis.PowerSpectrum.Y[1:,:]
  T = getTime(fileh5.root.Analysis.PowerSpectrum.Timing)[:,1]
  timeEvolution = []
  for step in range(len(data[0,:])):
    timeEvolution.append(data[:,step]/sum(data[:,step]))

  contourf(log10(D['ky']), T,timeEvolution, 250, locator=ticker.LogLocator(), cmap=cm.jet)
  colorbar()


def plotTurbulenceSpectra(dir='Y', pos=(-2,-1), fileh5=fileh[0], doFit='False', posT=(1,-1), field=0):
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
    D = getDomain(fileh5)
    
    doCFL    = kwargs.pop('doCFL', True)
    doFit    = kwargs.pop('doFit', True)

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


