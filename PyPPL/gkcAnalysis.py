import gkcData
import gkcStyle

import numpy as np
import pylab
import Filtering

def plotTimeEvolutionModePower(fileh5, **kwargs):
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
    ncol     = kwargs.pop('ncol', 3)  
    
    if doCFL == True : pylab.clf()
    
    T = gkcData.getTime(fileh5.root.Analysis.PowerSpectrum.Time)[2:,1]
    
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
      
      for m in modes:
            pl = pylab.semilogy(T, (scale*fileh5.root.Analysis.PowerSpectrum.Y[n_field, m,n_offset:]).T)
            if  (label == 'm' ) : legend_list.append("m = %i" % m)
            elif(label == 'ky') : legend_list.append("ky = %.1f" % (m / scale)) 
            else                : print "Name Error"
    
      leg = pylab.legend(legend_list, loc=leg_loc, ncol=ncol, mode="expand").draw_frame(0)
 
    else : raise TypeError("Wrong argument for dir : " + str(dir))
     
    pylab.xlabel("Time")
    pylab.xlim((0.,max(T)))


    if    field == "phi" : pylab.ylabel("Mode Power $|\\phi|^2$")
    elif  field == "A"   : pylab.ylabel("Mode Power $|A_\\parallel|^2$")
    elif  field == "B"   : pylab.ylabel("Mode Power $|B_\\parallel|^2$")
    else : raise TypeError("Wrong argument for field : " + str(field))
    
    #return pl, leg

def plotTimeEvolutionModePhase(fileh5, **kwargs):
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
    leg_loc  = kwargs.pop('loc', 'best')  
    dir      = kwargs.pop('dir', 'Y')
    modes    = kwargs.pop('modes' , range(D['Nky']))
    field    = kwargs.pop('field', 'phi')  
    n_offset = kwargs.pop('offset', 2)  
    label    = kwargs.pop('label', 'ky')  
    leg_loc  = kwargs.pop('loc', 'best')  
    ncol     = kwargs.pop('ncol', 2)  
    
    if doCFL == True : pylab.clf()
    
    T = gkcData.getTime(fileh5.root.Analysis.PhaseShift.Time)[2:,1]
    
    if   field == 'phi' : n_field = 0
    elif field == 'A'   : n_field = 1
    elif field == 'B'   : n_field = 2
    else : raise TypeError("Wrong argument for field : " + str(field))

    if(dir == 'X'):
      pl = plot(T, fileh5.root.Analysis.PhaseShift.X[n_field,:numModes,2:].T)
      legend_list = []
      for i in range(len(fileh5.root.Analysis.PhaseShift.X[n_field, :numModes,0])):
        legend_list.append("kx = %i" % i)
      leg = pylab.legend(legend_list, loc='lower right', ncol=2)
      leg.draw_frame(0)
    
    elif(dir == 'Y'):

      scale = fileh5.root.Grid._v_attrs.Ly/(2. * np.pi)
      
      legend_list = []
      
      for m in modes:

            data  = fileh5.root.Analysis.PhaseShift.Y[n_field, m,n_offset:]
            
            # set jump value to nan so it is not plotted, (can we speed up using ma ?)
            data_m = []
            for n in range(len(data[:-1])):
                if abs(data[n] - data[n+1]) < 1.: data_m.append(data[n])
                else                            : data_m.append(float('nan'))
            data_m.append(data[-1])
            data_m = np.array(data_m)

            print np.shape(T), np.shape(data_m)

            pl = pylab.plot(T, data_m)
            if  (label == 'm' ) : legend_list.append("m = %i" % m)
            elif(label == 'ky') : legend_list.append("ky = %.1f" % (m / scale)) 
            else                : print "Name Error"
    
      leg = pylab.legend(legend_list, loc=leg_loc, ncol=ncol, mode="expand").draw_frame(0)
 
    else : raise TypeError("Wrong argument for dir : " + str(dir))
    


    pylab.xlabel("Time")
    pylab.xlim((0.,max(T)))
    
    ax = pylab.gca()
    ax.set_yticks([-np.pi, -np.pi/2.,0.,np.pi/2., np.pi])
    ax.set_yticklabels(['$-\\pi$','$-\\pi/2$','$0$', '$\\pi/2.$', '$\\pi$'])
    pylab.ylim((-3.5,3.5))



    if    field == "phi" : pylab.ylabel("Mode Phase $|\\phi|^2$")
    elif  field == "A"   : pylab.ylabel("Mode Phase $|A_\\parallel|^2$")
    elif  field == "B"   : pylab.ylabel("Mode Phase $|B_\\parallel|^2$")
    else : raise TypeError("Wrong argument for field : " + str(field))
    
    return pl, leg



def plotInstantGrowthrates(fileh5, **kwargs):
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
    import scipy.ndimage 
    import scipy.interpolate
    
    D = gkcData.getDomain(fileh5)
    
    doCFL    = kwargs.pop('doCFL', True)
    dir      = kwargs.pop('dir', 'Y')
    modes    = kwargs.pop('modes' , range(D['Nky']))
    field    = kwargs.pop('field', 'phi')  
    n_offset = kwargs.pop('offset', 2)  
    label    = kwargs.pop('label', 'ky')  
    leg_loc  = kwargs.pop('loc', 'best')  
    off      = kwargs.pop('off', 2)  
    
    sigma      = kwargs.pop('sigma', 10)  
    #filterType = kwargs.pop('filterType', 'hanning')  
    
    if doCFL == True : pylab.clf()

    T = gkcData.getTime(fileh5.root.Analysis.PowerSpectrum.Time)[:,1]

    # TimeStep is irregular thus needs to cast into regular 
    def cast2Equidistant(T, Var):
        f = scipy.interpolate.interp1d(T, Var)
        Tnew = np.linspace(min(T), max(T), 1000)
        return Tnew, f(Tnew)

    if   field == 'phi' : n_field = 0
    elif field == 'A'   : n_field = 1
    elif field == 'B'   : n_field = 2
    else : raise TypeError("Wrong argument for field : " + str(field))
    
    

    if(dir == 'X'):

      plot(T, grad)
      legend_list = []
      for i in range(len(fileh5.root.Analysis.PowerSpectrum.X[n_field,:numModes,0])):
        legend_list.append("kx = %i" % i)
      
      leg = legend(legend_list, loc='best', ncol=4)
      leg.draw_frame(0)
    
    if(dir == 'Y'):
      for m in modes:
        Var =   fileh5.root.Analysis.PowerSpectrum.Y[n_field,m,off:]
        Tn, V = cast2Equidistant(T[off:], np.log10(Var))
        # Use NdImage for line-smoothening (convolution with gaussian kernel)
        V   = scipy.ndimage.gaussian_filter(V, sigma=sigma, mode='nearest')
        gamma = np.gradient(V, Tn[1]-Tn[0])
        pylab.plot(Tn, gamma, "-", label='ky = %.2f' % (D['ky'][m]))#, color=gkcStyle.markers_C[m])
      
      leg = pylab.legend(loc='best', ncol=4)
      leg.draw_frame(0)
    
    
    pylab.xlabel("Time")
    pylab.xlim((0.,max(T)))

    pylab.ylabel("Instant Mode Growth $\gamma$")
        

def plotCrossCorrelateValues(fileh5, A="Phi",B="Tp",frame=-1, Z=0, species=1):
    D = gkcData.getDomain(fileh5)

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
    X = np.linspace(-pi, pi, len(Corr[0,:]))
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
    ax.plot(np.linspace(0., 0., 101), np.linspace(Y.min(), Y.max(), 101), 'r-')
    ylim((Y.min(), Y.max()))
    xlabel("Phase")
    ylabel("$k_y$")
    #xlim((X.min(), X.max()))
    # plot zero phase line


################### Ploshed functions


def plotContour(fileh5, var="2DPhi", **kwargs):
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
    D = gkcData.getDomain(fileh5)
    
    norm     = kwargs.pop('Normalize', True)
    Z        = kwargs.pop('Z', 0)
    modes    = kwargs.pop('modes' , range(D['Nky']))
    doCFL    = kwargs.pop('doCFL' , True)
    interpolation = kwargs.pop('interpolation' , 'bilinear')
    printTitle    = kwargs.pop('printTitle' , True)
    orientation   = kwargs.pop('orientation' , 'vertical')
    frame         = kwargs.pop('frame' , -1)
    
    
    D, T, data = gkcData.getData(var, fileh5, Z, frame, species=0)
  
    X, Y, data = gkcData.getRealFromXky(fileh5, data, modes)
   
    gkcStyle.plotContourWithColorbar(X,Y, data, norm=norm, interpolation=interpolation, orientation=orientation)
   
   
    pylab.xlabel("Radial Direction")
    pylab.ylabel("Poloidal Direction")
    if printTitle == True : pylab.title("TimeStep : %i   Time : %.3f " % (T[0], T[1]))




def plotScalarDataTimeEvolution(fileh5, var = "KinPhi", **kwargs):
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
    log      = kwargs.pop('log', True)
    
    leg_loc  = kwargs.pop('leg_loc' , 'best')
    leg_ncol = kwargs.pop('leg_ncol', 1)

    if(log == 1) : 
      plotf = pylab.semilogy
      af    = abs
    else         : 
      plotf = pylab.plot
      af    = lambda x : x
        
    T = fileh5.root.Analysis.scalarValues.cols.Time[5:]
    n_sp = len(fileh5.root.Species.cols.Name[1:]) 



    # Plot Field Energy"
    if var == "phi":
        plotf(T, af(fileh5.root.Analysis.scalarValues.cols.phiEnergy[3:]), label="Field Energy")
        pylab.ylabel("Field Energy")

    if var == "kinEnergy":
        for s in range(n_sp):
            plotf(T, af(fileh5.root.Analysis.scalarValues.cols.KineticEnergy[3:][:,s]) , label="Kinetic Energy (" + fileh5.root.Species.cols.Name[s]+ ")" )

        #plotf(T, af(fileh5.root.Analysis.scalarValues.cols.phiEnergy[3:]+1.*(np.sum(fileh5.root.Analysis.scalarValues.cols.KineticEnergy[3:], axis=1))) , label="Total Energy")
        pylab.ylabel("Energy")
    
    if var == "KinPhi":
        pylab.plot(T, fileh5.root.Analysis.scalarValues.cols.phiEnergy[5:], label="Field Energy")
        E_kin_total = fileh5.root.Analysis.scalarValues.cols.phiEnergy[5:]
        for s in range(n_sp):
          E_kin = fileh5.root.Analysis.scalarValues.cols.KineticEnergy[5:][:,s]-fileh5.root.Analysis.scalarValues.cols.KineticEnergy[5][s]
          E_kin = -abs(E_kin)
          E_kin_total = E_kin_total + E_kin
          pylab.plot(T, E_kin , label="Kinetic Energy (" + fileh5.root.Species.cols.Name[s]+ ")" )
        gkcStyle.plotZeroLine(min(T), max(T))
        pylab.plot(T, -abs(E_kin_total)/100. , label="Total Energy")
        pylab.yscale('symlog', linthreshy=1.e-7)
        # Plot Toal nergy
        pylab.ylabel("Energy")

    if var == "heatFlux":
  
       for s in sp:
        for s in range(n_sp):
          plot(T, fileh5.root.Analysis.scalarValues.cols.HeatFlux[3:][:,s], label=fileh5.root.Species.cols.Name[s]) 
       pylab.ylabel("Heat Flux")
 
    if var == "Charge":
       charge = np.zeros(len(fileh5.root.Analysis.scalarValues.cols.ParticleNumber[3:][:,1]))
       for s in range(n_sp):
          dn = (fileh5.root.Analysis.scalarValues.cols.ParticleNumber[3:][:,s] - fileh5.root.Analysis.scalarValues.cols.ParticleNumber[3][s])
          n  = fileh5.root.Analysis.scalarValues.cols.ParticleNumber[3:][:,s] 
          q = fileh5.root.Species.cols.Charge[s]
          charge = charge + q * n

          pylab.plot(fileh5.root.Analysis.scalarValues.cols.Time[3:], q*n , label=fileh5.root.Species.cols.Name[s])
          pylab.plot(fileh5.root.Analysis.scalarValues.cols.Time[3:], q*dn , label=fileh5.root.Species.cols.Name[s])
       if D['Ns'] > 1:
          plotf(fileh5.root.Analysis.scalarValues.cols.Time[3:], charge , label="Charge")
       
       gkcStyle.plotZeroLine(min(T), max(T))
       pylab.yscale('symlog', linthreshy=1.e-7)
       pylab.ylabel("Charge")
  
    if var == "particleFlux":
      for s in range(n_sp):
        plotf(T, fileh5.root.Analysis.scalarValues.cols.ParticleFlux[3:][:,s], label=fileh5.root.Species.cols.Name[s]) 
      pylab.ylabel("Particle Flux")
 
    pylab.xlim((min(T), max(T)))
    leg = pylab.legend(loc=leg_loc, ncol=leg_ncol).draw_frame(0)
    pylab.xlabel("Time")


def plotTurbulenceTime(fileh5, dir='Y', pos=(-2,-1), doFit='False', posT=(1,-1)):
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
    data = fileh5.root.Analysis.PowerSpectrum.Y[1:,:]
    T = gkcData.getTime(fileh5.root.Analysis.PowerSpectrum.Timing)[:,1]
    timeEvolution = []
    for step in range(len(data[0,:])):
       timeEvolution.append(data[:,step]/sum(data[:,step]))

    contourf(np.log10(D['ky']), T,timeEvolution, 250, locator=ticker.LogLocator(), cmap=cm.jet)
    colorbar()


def plotTurbulenceSpectra(fileh5, dir='Y', start=1, end=-1, posT=(1,-1), field=0, **kwargs):
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
    doFit    = kwargs.pop('doFit', True)
    m        = kwargs.pop('m', 0)

    if(dir == 'X'):
      data = sum(fileh5.root.Analysis.PowerSpectrum.X[field,0:,start:end], axis=1)/abs(end-start)
      loglog(D['kx'], data/sum(data), '.-')
      xlabel("$k_x$")
    elif(dir == 'Y'):
      data = np.mean(fileh5.root.Analysis.PowerSpectrum.Y[field, 0:,start:end], axis=1)
      mypl = pylab.loglog(D['ky'][1:], data[1:], gkcStyle.markers_D[m] +  '-', color=gkcStyle.markers_C[m], markersize=8.)
      
      # Plot Zonal flow seperately
      pylab.loglog(0.9*D['ky'][1], data[0], gkcStyle.markers_D[m], markersize=8., color=gkcStyle.markers_C[m])

      pylab.xlim((0.75 * D['ky'][1], 1.25*D['ky'][-1]))
      pylab.xlabel("$k_y$")
     

      # Draw vertical line
      if m == 0:
        min_x = min(data)
        max_x = max(data)
        pylab.loglog(np.linspace(0.9*D['ky'][1], 0.9*D['ky'][1], 201), np.logspace(np.log10(0.8*min_x), np.log10(1.2*max_x), 201), "-", linewidth=1.5, color="#666666")
        #pylab.text(D['ky'][1], 1.05*np.sqrt(max_x), "Zonal Flow", weight='bold', color="#666666", rotation='vertical')
        pylab.text(0.95*D['ky'][1], 5*min_x, "Zonal Flow", weight='bold', color="#666666", rotation='vertical')

    elif(dir == 'Z'):
      data = sum(fileh5.root.Analysis.PowerSpectrum.Z[field, 0:,start:end], axis=1)/abs(start-end)
      pylab.loglog(D['kp'], data[1:], '.-')
      pylab.loglog(D['kp'][0], data[0], 'o')
      xlabel("$k_z$")
      
    pylab.ylabel("$|\\phi_k(k_y)|^2$")
        
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
        text(log(D['ky'][5]), log(data[5]), "$\\propto k_y^{%2.1f}$" % p1[0], ha="center", size=14)
    
    return mypl


def plotXPropTempDensity(fileh5, time=0, Z=4, Y=6):

  markers = ['v-c', '^-y', '<-r', '>-m', '*-b', 'd-g', 'p-b', '1-r', '2-m']
  pl = []
  pl_name = []
  D = gkcData.getDomain(fileh5)
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

    

  leg = pylab.legend(pl, pl_name, loc='lower center')
  leg.draw_frame(0)


def plotModeStructure(fileh5, mode, part = "r",  **kwargs):
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
    field    = kwargs.pop('field', 'phi')  
    Z        = kwargs.pop('Z', 0)  
    frame    = kwargs.pop('frame', -1)  
    label    = kwargs.pop('label', "")  
    m        = kwargs.pop('m', 0)  

    
    if   field == 'phi' : data_X = fileh5.root.Visualization.Phi[Z,mode,:,frame]
    elif field == 'A'   : n_field = 1
    elif field == 'B'   : n_field = 2
    else : raise TypeError("Wrong argument for field : " + str(field))

    # Normaliza
    if   part == "a" : data_X = abs(data_X) / np.sum(abs(data_X))
    elif part == "r" : data_X = np.real(data_X) / np.sum(abs(np.real(data_X)))
    elif part == "i" : data_X = np.imag(data_X) / np.sum(abs(np.imag(data_X)))
    else : raise TypeError("Wrong argument for part : " + str(part))

    pylab.plot(D['X'], data_X, gkcStyle.markers_C[m], label=label)

    if    field == "phi" : pylab.ylabel("Mode Power $|\\phi|^2$")
    elif  field == "A"   : pylab.ylabel("Mode Power $|A_\\parallel|^2$")
    elif  field == "B"   : pylab.ylabel("Mode Power $|B_\\parallel|^2$")
    else : raise TypeError("Wrong argument for field : " + str(field))

    pylab.xlim((min(D['X']), max(D['X'])))
    pylab.xlabel("X")

    
