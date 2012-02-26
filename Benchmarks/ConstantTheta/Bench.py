"""
    Simple Benchmark for Drift kinetic, Dorland, Gyro-kinetic comparioson.

Take care of nomralization


"""

markers_D = [ '<r', 'sk', 'vb', '^c', '>g', '^m', 'dy', 'kx', 'kx', 'kx']
markers_L = [ 'r-', 'k-', 'b-', 'c-', 'g-', 'm-', 'y-', 'k-', 'k-', 'k-' ]

from Execute_Thin import *
from sets import Set
import string 

N = 5


def linspace(start, stop, count):
	x = []
	for v in range(count) : x.append(start+v/float(count-1)*(stop-start))
        return x

shear = 0.025


# sysargv is input paraeters file
if sys.argv[1] == "Exe" : 

     decomposition = "1:1:1:1:16"
     helios_options = []
     
     for theta in [ 0.033, 0.05, 0.133, 0.15, 0.2]: 
        helios_options.append("\"DataOutput.OutputFileName=Nakata11_Theta_%.2f.h5;Geometry.Theta=%.2f\"" % (theta, theta))
     executeThin(sys.argv[2], helios_options, decomposition, queue='lh10247', system_shared=False)
     

if sys.argv[1] == "Ana":
    from tables import *
    from FitFunction import *
    import FitFunction
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import interpolate
    import Dispersion_ConstTheta 

    n_start = int(sys.argv[2])
    n_stop  = int(sys.argv[3])	

    #determine Type

    result_A = [] ;
    result_Theory = [] ;
    ky = 0.
    ky_list = 0.
    for fileh_name in sys.argv[4:]:
          try:
            fileh5  = openFile(fileh_name)

            ky, gamma = FitFunction.getGrowth(fileh5, pos=(n_start,n_stop))           
            gamma = array(gamma) #fileh5.root.Species[1]          
            result_A.append((gamma, fileh5.root.Geometry._v_attrs.Theta[0]))
  
            #Setup = { 'eta_e' : 6., 'kx' : 0.0, 'v_te' : sqrt(2.), 'rho_te2' : 1., 'tau' : 1., 'theta' : 0. }
            Setup = { 'eta' : 6., 'kx' : 0.0, 'v_te' : mp.sqrt(2.), 'rho_te2' : 1., 'tau' : 1., 'theta' : 0. , 'm_ie' : 1837., 'lambda_D2' : 0.}
            Setup['theta'] = fileh5.root.Geometry._v_attrs.Theta[0]
            ky_list = linspace(min(ky), min(max(ky),9999.), 151)

            _disp="Gyro1st"
            if(fileh5.root.Grid._v_attrs.Nm[0] > 1) : _disp="Gyro"
            print _disp
            print Setup['theta'], " " , _disp
            kye, y, e = Dispersion_ConstTheta.getGrowth(ky_list, Setup, disp=_disp)
            result_Theory.append((kye, y, e, Setup['theta']))
  
            fileh5.close()
    
          except : #  BaseException as Error :
               print "Cannot open : ", fileh_name, "   (skipping)"
               #print "Cannot open : ", fileh_name, "   (skipping)", Error.message, " " , Error.args
    pl = []
    n = 0
    #fig = subplot(111).set_size_inches((8,8))
    for nw in result_A: 
      #pl.append(plt.semilogx(ky, nw[0], markers_D[n], markersize=12., label=nw[1]))
      pl.append(plt.semilogx(ky, nw[0], markers_D[n], markersize=9., label="%.3f" % nw[1]))
    
  
    # Get Setting
  #pylab.plot(results[:,0], results[:,1], 'r', label='real')
      """
      try:
            # Interpolate result_A A
            tck = interpolate.splrep(ky,nw[0],s=0)
            xnew = linspace(ky.min(), ky.max(), 201)
            ynew = interpolate.splev(xnew,tck,der=0)
            #plt.semilogx(xnew, ynew, markers_L[n])
            plt.plot(xnew, ynew, markers_L[n])
      except: print "Interpolation Failed"
      """ 
      n = n + 1
    n = 0
    for nw in result_Theory: 
      plt.semilogx(nw[0], imag(nw[1]), markers_L[n])
      #plt.plot(nw[0], real(nw[1]), markers_L[n])
      n = n + 1
    
    plt.legend(ncol=5).draw_frame(0)
    plt.xlabel("Poloidal Wavenumber $k_y/\\rho_{te}$")
    plt.ylabel("Growthrate $\\gamma_L [\\nu_{te}/L_T]$")
    
    #plt.xlim((0., 2.0))
    plt.ylim((-0.05, 0.25))

    #plt.show() 
    plt.savefig("ModeGrowth.png", bbox_inches='tight')
    plt.savefig("ModeGrowth.pdf", bbox_inches='tight')

    clf()
    # Plot here Theory Erro
    n = 0
    for nw in result_Theory:
      print nw[2]
      plt.semilogy(nw[0], nw[2], markers_L[n], label="%.3f" % nw[3])
      n = n + 1
    
    #plt.legend(ncol=5).draw_frame(0)
    plt.xlabel("Poloidal Wavenumber $k_y/\\rho_{te}$")
    plt.ylabel("Residuum Growthrate $\\gamma_L [\\nu_{te}/L_T]$")
    
    #plt.xlim((0., 2.0))

    #plt.show() 
    plt.savefig("ModeGrowth_Res.png", bbox_inches='tight')
    plt.savefig("ModeGrowth_Res.pdf", bbox_inches='tight')


    # Plot here Theory Erro










