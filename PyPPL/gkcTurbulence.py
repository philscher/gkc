
def plotFrequencySpectra(fileh5, which='b', markline="-", **kwargs):
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
    import gkcData
    import gkcStyle

    import pylab
    import numpy as np
    
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
    
    if   field == 'phi' : n_field = 0
    elif field == 'A'   : n_field = 1
    elif field == 'B'   : n_field = 2
    else : raise TypeError("Wrong argument for field : " + str(field))

    if(dir == 'X'):
        raise TypeError("Not implemented for X-direction")
    elif(dir == 'Y'):
        
      shift = fileh5.root.Analysis.PhaseShift.Y   [n_field, :,start:stop]
      print "Using data from T=", T[start], " to T = ", T[stop]
      freq = [] 
      for nky in np.arange(1,len(D['ky'])):
        time_series = np.sin(shift[nky,:])
        FS = np.fft.rfft(time_series)
      
        # Not only if we have constant time step !
        freq.append(abs(FS))
        fftfreq = np.fft.fftfreq(len(abs(np.fft.fftshift(FS))), d = (T[-10]-T[-11])) 
    
      freq = np.array(freq)
      print np.shape(fftfreq), " 2 : ", np.shape(D['ky'][1:]), np.shape(freq)
      pylab.contourf(D['ky'][1:], np.fft.fftshift(fftfreq), freq.T, 100, cmap=pylab.cm.jet) 
      pylab.xlim((D['ky'][1], D['ky'][-1]))
      pylab.ylim((min(fftfreq), max(fftfreq)))
      pylab.colorbar()
    else : raise TypeError("Wrong argument for dir : " + str(dir))
    pylab.gca().set_xscale("log") 
     
    pylab.xlabel("$k_y$")
   
    gkcStyle.plotZeroLine(D['ky'][1], D['ky'][-1], color='r')

    pylab.ylabel("Frequency $\\omega_r(k_y)$")
    
    #return pl, leg




def plotAveragedFluxes(fileh5, var='H', **kwargs):
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
    import gkcData
    import gkcStyle

    import pylab
    import numpy as np
    
    D = gkcData.getDomain(fileh5)

    doCFL    = kwargs.pop('doCFL', True)
    start    = kwargs.pop('start', 1)  
    stop     = kwargs.pop('stop', -1)  
    scale    = kwargs.pop('scale', 'ky')  
    field    = kwargs.pop('field', 'phi')  
    

    if doCFL == True : pylab.clf()
    
    if   field == 'phi' : n_field = 0
    elif field == 'A'   : n_field = 1
    elif field == 'B'   : n_field = 2
    else : raise TypeError("Wrong argument for field : " + str(field))
    
    T = gkcData.getTime(fileh5.root.Analysis.PowerSpectrum.Time)[:,1]

    # Averaged over Z and start stop
    if   var == 'P' : 
      data = np.mean(fileh5.root.Analysis.Flux.Density[n_field,:,:,start:stop], axis=2)
      ylabel = "Particle Flux $k_y \\Gamma/\\Gamma_{gB}\,(ky)$"
    elif var == 'H' : 
      data = np.mean(fileh5.root.Analysis.Flux.Heat   [n_field,:,:,start:stop], axis=2)
      #ylabel = "Heat Flux $Q/Q_\\textrm{gB}(k_y)$"
      ylabel = "Heat Flux $k_y Q/Q_{gB}\,(k_y)$"
    else : raise TypeError("No such variable")
     
    if scale == 'ky' : gor_ky = D['ky'][1:]
    else             : gor_ky = 1.

    for s in range(D['Ns']):
      species_name = fileh5.root.Species.cols.Name[s+1]
      pylab.semilogx(D['ky'][1:], gor_ky * data[1:,s], gkcStyle.markers_D[s]+'-', label=species_name, color=gkcStyle.markers_C[s], markersize=7.)
    
    pylab.xlabel("$k_y \\rho_i$")
    pylab.ylabel(ylabel)
  
    pylab.legend(ncol=2).draw_frame(0)
    pylab.xlim((0.8*D['ky'][1], 1.2*D['ky'][-1]))
    




