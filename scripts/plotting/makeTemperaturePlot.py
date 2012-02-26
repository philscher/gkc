
from enthought.mayavi import mlab
import tables
import string
import sys
from numpy import *


fileh = tables.openFile(sys.argv[1])



fmax = 0.0
#//for frame in range(len(fileh.root.Time[:])):
#    fmax = max([fmax, fileh.root.Phi[:,2:-2,2:-2, frame].max() ])


X = linspace(0.0,32., len(fileh.root.Phi[:,0,0,0]))
Y = X/2
Z = 8000./32. * X



ext = [ min(X), max(X), min(Y), max(Y), 0.0, fmax ]



fig3d = mlab.figure(bgcolor=(1.0, 1.0, 1.0), fgcolor=(0.0, 0.0, 0.0), size=(1280, 1280))

angle = 0.0

def Temperature(f, v):
  temperature = zeros( shape(f[0,:,:,:]))

  print "Shape Temp: ", shape(temperature)
  print "Shape v   : ", shape(v), v


  for i in range(len(f[0,:,0,0])):
    for j in range(len(f[0,0,:,0])):
      for k in range(len(f[0,0,0,:])):
        # Calculate kinetic energy and norm per particle
        temperature[i,j,k] = (f[:,i,j,k]* v[2:-2]**2).sum() / f.sum(axis=0)[i,j,k]


  return temperature


for frame in range(len(fileh.root.Time[:])):
  # 3D Movie
  frame = 4*frame
  mlab.clf(figure=fig3d)
  scene = mlab.contour3d(Temperature(fileh.root.Phasespace[:,:,:,:,frame], fileh.root.Phasespace.attrs.v), contours=50)
  # Density profile
  #scene = mlab.contour3d(fileh.root.Phasespace[:,:,:,:,frame].sum(axis=0), contours=50)
  #source = mlab.pipeline.scalar_field(Temperature(fileh.root.Phasespace[:,:,:,:,frame], fileh.root.Phasespace.attrs.v))
  #vol = mlab.pipeline.volume(source)
  angle = 50.
  mlab.view(azimuth=angle%360, distance=72.0)
  #mlab.axes(color=(0.0,0.0,0.),extent=ext,xlabel='X',ylabel='Y',zlabel='Z')
  #a = array(gradient(fileh.root.Phi[:,:,:,frame]))**2
  #mlab.quiver3d(a[0],a[1],a[2])
  mlab.colorbar(orientation='vertical')#, label_fmt="%.2f", nb_labels=3)
  mlab.outline()
  mlab.savefig("Rho_" + string.zfill(frame, 4) + ".jpg")#, size=(600,600))
  print "[", string.zfill(frame, 4),"/", len(fileh.root.Time[:]), "]"





#mencoder "mf://*.jpg" -mf fps=10 -ovc lavc -o mymovie.av
