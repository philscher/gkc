from pylab import *
import sys
import os
import streamlines
from StringIO import StringIO   # StringIO behaves like a file object
from tables import *
import string 

def getXYFromDataFile(fileh5, frame, addPsi0=True):
    # Write in Data
    Psi_1k = fileh5.root.Field.Psi[:,:,frame]

    nky =  len(Psi_1k[:,0])
    nx =  len(Psi_1k[0,:])
    #print "Nky : " , nky, " Nx : " , nx

    # reshape data into (ky, x)
    #Psi_1k = complex(1.,0.) * data[:,4].reshape(nky,nx) + complex(0.,1.) * data[:,5].reshape(nky,nx)

    # Resize to have higher resolution in y
    Psi_1k = resize(Psi_1k,(129, nx))
    Psi_1k[nky:,:] = 0.

    # Fourier back transform c2r, renormaliza
    Psi = 128.*irfft(Psi_1k, axis=1)

    Psi0 =  fileh5.root.Init.Psi0[:,1]
    if(addPsi0 == True): Psi = Psi + Psi0
    
    # shift array by half in y !
    PsiA = Psi[128:,:]
    Psi [128:,:] = Psi[:128,:]
    Psi[:128,:]  = PsiA[::-1,:]

    x = fileh5.root.Init._v_attrs.X[:]
    #print x
    #return (x[::8], Psi0[::8], Psi[:,::8])
    return (x[:], Psi0[:], Psi[:,:])



if sys.argv[1] == 'Plot':
 
  fileh5 = openFile(sys.argv[2])

  n = 1
  for frame in range(len(fileh5.root.Field.Psi[0,0,:]))[:]:
    clf()
    T = fileh5.root.Field.Time[frame]
    x, Psi0, Psi = getXYFromDataFile(fileh5, frame, True)
    y = linspace(  0, 5., len(Psi[:,0]))

    #contourf(x,y,Psi,100)
    contourf(linspace(-5,5, len(Psi[0,:])),y,Psi,100)

    #dA_dy, dA_dx = gradient(Psi[:,::8])
    #streamlines.Streamlines(x[::8],y, -dA_dy, dA_dx, spacing=5., res=2.).plotArrows(color="#888888", size=40., lw=1.1)
    colorbar()
    xlim(min(x), max(x))
    ylim(min(y), max(y))
    xlabel("X")
    ylabel("Y")
    title("Step : %7.0f Time : %7.2f" % ( T[0], T[1])) 
    savefig(os.path.splitext(sys.argv[2])[0] + '_' + string.zfill(n,5) + '.png', bbox_inches='tight')
    print "Plotted contour : " , sys.argv[2] , " (", n, "/", str(len(fileh5.root.Field.Psi[0,0,:])), ")"
    n += 1



if sys.argv[1] == 'Convert':

  x, Psi = getXYFromDataFile(sys.argv[3], sys.argv[2], False)
  yf     = linspace( 0., 20.*pi, len(Psi[:,0])) 
  x      = linspace(-10.*pi, 10.*pi, len(x)) 


  # write file
  outputfilename = "MagneticIsland.dat"
  
  ofile = open("MagneticIsland.dat", 'w')



  def getString1D(A):
    str_A = ''
    for nx in range(len(A[:])): str_A = str_A + '%.3e ' % A[nx]
    return str_A + '\n'

  ofile.write(getString1D(x))
  ofile.write(getString1D(yf))

  for z in range(len(Psi[:,0])): ofile.write(getString1D(Psi[z,:]))


  ofile.close()



def getIslandWidth(x,y,Psi):
    from scipy import interpolate
    import mpmath as mp      
    Nx = len(x)
    Ny = len(y)
    XPointValue_1 = Psi[0,Nx/2]
    XPointValue_2 = Psi[-1,Nx/2]

    
    # Interpolate Psi
    Psi_inter = interpolate.splrep(x,Psi[Ny/2,:],s=0)

    # Find zero
    def W(x) : 
      return XPointValue_1 - interpolate.splev(float(x),Psi_inter,der=0)
    try:
        w1 = float(mp.findroot(W,  1))
        w2 = float(mp.findroot(W, -1))
    except:
        w1 = 0.
        w2 = 0.
    return (w1, w2) 


if sys.argv[1] == 'PlotWidth':

    width = []
    filename = sys.argv[3]
    x, Psi = getXYFromDataFile(fileh5, -1, True)
    y = linspace(  0, 5., len(Psi[:,0]))

    w1, w2  = getIslandWidth(x,y,Psi)
    print "Island Width is : ",  (w1-w2)/2., " +- ", abs(w1+w2)/2. 
    
    fig = plt.figure(1, figsize=(12,12))
    fig.clf()
    ax = fig.add_subplot(111)
    contourf(x[::1],y,Psi[:,::1],100)

    T = loadtxt(StringIO(open(filename).readline()))[:]
    dA_dy, dA_dx = gradient(Psi)
    streamlines.Streamlines(x,y, -dA_dy, dA_dx, spacing=20., res=1.).plotArrows(color="#888888", size=100., lw=1.1)
    colorbar()
    xlim(min(x), max(x))
    ylim(min(y), max(y))
    xlabel("X")
    ylabel("Y")
    title("Step : %7.0f Time : %7.2f" % ( T[0], T[1]))

    # Plot Island width Annotae

    #ax.annotate("", xy=(w2, 2.5), xycoords='data', xytext=(w1, 2.5), textcoords='Island', \
    #            arrowprops=dict(arrowstyle="->", connectionstyle="arc3"),)

    ax.annotate('', xy=(w1, 2.5), xytext=(0., 2.5), size=140., arrowprops=dict(facecolor='white', shrink=0.05),)
    ax.annotate('', xy=(w2,2.5), xytext=(0., 2.5), size=140., arrowprops=dict(facecolor='white', shrink=0.05),)
    #ax.annotate('', xy=(0., 2.5), xytext=(w2, 2.5), arrowprops=dict(facecolor='white', shrink=0.05, arrowstyle='-|'),)

    fig.savefig('Island.png', bbox_inches='tight')



if sys.argv[1] == 'WidthScan':


 
    width = []
    for filename in sys.argv[3:]:
      try:     
        fileh5 = openFile(filename)
   
        x, Psi0, Psi = getXYFromDataFile(fileh5, -1, True)
        y = linspace(  0, 5., len(Psi[:,0]))
        Ly = fileh5.root.Init._v_attrs.Ly[0]
        Lx = fileh5.root.Init._v_attrs.Lx[0]
        w1, w2  = getIslandWidth(x,y,Psi)
        width.append(((w1-w2)/2., Ly, Lx) )
      except BaseException as Error : print " Error : " , Error.message
    print width
    width = array(width)


    fig = plt.figure(1, figsize=(12,12))
    fig.clf()
    ax = fig.add_subplot(111)
    
    plot(width[:,2], width[:,0], 's-')
    
    #ylim(min(y), max(y))
    xlabel("Lx")
    ylabel("Island Width")

    # Plot Island width Annotae

    #ax.annotate("", xy=(w2, 2.5), xycoords='data', xytext=(w1, 2.5), textcoords='Island', \
    #            arrowprops=dict(arrowstyle="->", connectionstyle="arc3"),)

    fig.savefig('Island.png', bbox_inches='tight')



