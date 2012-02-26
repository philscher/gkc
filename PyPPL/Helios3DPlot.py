from HeliosPlot import *


from enthought.mayavi import mlab
import tables
import string
import sys
import os
from numpy import *



def plotPhi3D(



fileh = tables.openFile(sys.argv[1])

frame_skip = 1
if len(sys.argv) > 2:
    frame_skip = int(sys.argv[2])

filename = os.path.splitext(os.path.basename(sys.argv[1]))[0]

fmax = 0.0
fmin = 0.0

for frame in range(len(fileh.root.Time[:])):
  fmax = max([fmax, (fileh.root.Phi[2:-2,2:-2, :, frame]**2).max() ])
  fmin = min([fmin, ((fileh.root.Phi[2:-2,2:-2, :, frame]**2)).min() ])
# Set to zero because we take the square
fmin = 0.0


print "Determined Min/Max Values ....  min : ", fmin, "  max : ", fmax


X = linspace(-1.0,1.0, len(fileh.root.Phasespace.attrs.x[2:-2]))
Y = X
Z = 1.5 * X


ext = [ min(X), max(X), min(Y), max(Y), min(Z), max(Z) ]



fig3d = mlab.figure(bgcolor=(1.0, 1.0, 1.0), fgcolor=(0.0, 0.0, 0.0), size=(1280, 1280))

angle = 1.0*25


for frame in range(len(fileh.root.Time[1:])):
  frame=frame_skip*frame
  # 3D Movie
  
  mlab.clf(figure=fig3d)
  #Simple contour
  scene = mlab.contour3d((fileh.root.Phi[2:-2,2:-2,2:-2,frame]**2), contours=50, extent=ext, vmax=fmax, vmin=fmin)#, colormap='ylgn')

  # Cutplane
  #E  = gradient(fileh.root.Phi[:,:,:,frame])
  #mlab.quiver3d(E[0], E[1], E[2])
  #src = mlab.pipeline.vector_field(E[2], E[1], E[0])
  #mlab.pipeline.vectors(src, mask_points=20, scale_factor=3.)
  #mlab.pipeline.vector_cut_plane(src, mask_points=1, scale_factor=5)
  #src = mlab.pipeline.vector_field(E[0], E[1], E[2])
  #magnitude = mlab.pipeline.extract_vector_norm(src)
  #flow = mlab.pipeline.streamline(magnitude, seedtype='plane', seed_visible=False, seed_scale=1.0, seed_resolution=100, linetype='line',integration_direction='both')
  #$mlab.axes(color=(0.0,0.0,0.),extent=ext,xlabel='X',ylabel='Y',zlabel='Z')
  #source = mlab.pipeline.scalar_field(fileh.root.Phi[:,:,:,frame]**2)
  #vol = mlab.pipeline.volume(source)
  #a = array(gradient(fileh.root.Phi[:,:,:,frame]))**2
  #mlab.quiver3d(a[0],a[1],a[2])
  #  mlab.view(azimuth=angle%360, distance=128.0)
  mlab.view(azimuth=angle%360, distance=10.0)
  angle += 1.0
  mlab.colorbar(orientation='vertical', label_fmt="%.2e", nb_labels=5)
  mlab.outline()
  mlab.savefig("Phi3D_" + filename + "_" + string.zfill(frame, 4) + ".jpg")#, size=(600,600))
  print "[", string.zfill(frame, 4),"/", len(fileh.root.Time[:]), "]"


#mencoder "mf://*.jpg" -mf fps=10 -ovc lavc -o mymovie.av

from enthought.mayavi import mlab
import tables
import string
import sys
from numpy import *


fileh = tables.openFile(sys.argv[1])



fmax = 0.0
for frame in range(len(fileh.root.Time[:])):
  fmax = max([fmax, fileh.root.Phi[2:-2,2:-2, :, frame].max() ])


X = linspace(0.0,32., len(fileh.root.Phi[:,0,0,0]))
Y = X/2
Z = 8000./32. * X



ext = [ min(X), max(X), min(Y), max(Y), 0.0, fmax ]



fig3d = mlab.figure(bgcolor=(1.0, 1.0, 1.0), fgcolor=(0.0, 0.0, 0.0), size=(1280, 1280))

angle = 0.0

def Temperature(f, v):
  temperature = zeros( shape(f[0,:,:,:]))

  print "Shape Temp: ", shape(temperature)
  print "Shape v   : ", shape(v)


  for i in range(len(f[0,:,0,0])):
    for j in range(len(f[0,0,:,0])):
      for k in range(len(f[0,0,0,:])):
        # Calculate kinetic energy and norm per particle
        temperature[i,j,k] =  f.sum(axis=3)[i,j,k]


  return temperature


for frame in range(len(fileh.root.Time[:])):
  # 3D Movie
  #frame = 25*frame+150
  mlab.clf(figure=fig3d)
  #scene = mlab.contour3d(Temperature(fileh.root.Phasespace[:,:,:,:,frame], fileh.root.Phasespace.attrs.v), contours=50)
  # Density profile
  scene = mlab.contour3d(fileh.root.Phasespace[:,:,:,:,frame].sum(axis=0), contours=50)
  #source = mlab.pipeline.scalar_field(Temperature(fileh.root.Phasespace[:,:,:,:,frame], fileh.root.Phasespace.attrs.v))
  #vol = mlab.pipeline.volume(source)
  angle = 50.
  mlab.view(azimuth=angle%360)
  #mlab.axes(color=(0.0,0.0,0.),extent=ext,xlabel='X',ylabel='Y',zlabel='Z')
  #a = array(gradient(fileh.root.Phi[:,:,:,frame]))**2
  #mlab.quiver3d(a[0],a[1],a[2])
  mlab.colorbar(orientation='vertical')#, label_fmt="%.2f", nb_labels=3)
  mlab.outline()
  mlab.savefig("Rho_" + string.zfill(frame, 4) + ".jpg")#, size=(600,600))
  print "[", string.zfill(frame, 4),"/", len(fileh.root.Time[:]), "]"





#mencoder "mf://*.jpg" -mf fps=10 -ovc lavc -o mymovie.av

from enthought.mayavi import mlab
import tables
import string
import sys
import os
from numpy import *


fileh = tables.openFile(sys.argv[1])

frame_skip = 1
if len(sys.argv) > 2:
    frame_skip = int(sys.argv[2])

filename = os.path.splitext(os.path.basename(sys.argv[1]))[0]

fmax = 0.0
fmin = 0.0
for frame in range(len(fileh.root.Time[:])):
  fmax = max([fmax, (fileh.root.Phi[2:-2,2:-2, :, frame]**2).max() ])
  fmin = min([fmin, ((fileh.root.Phi[2:-2,2:-2, :, frame]**2)).min() ])
# Set to zero because we take the square
fmin = 0.0


print "Determined Min/Max Values ....  min : ", fmin, "  max : ", fmax


X = linspace(-1.0,1.0, len(fileh.root.Phasespace.attrs.x[2:-2]))
Y = X
Z = 1.5 * X


ext = [ min(X), max(X), min(Y), max(Y), min(Z), max(Z) ]



fig3d = mlab.figure(bgcolor=(1.0, 1.0, 1.0), fgcolor=(0.0, 0.0, 0.0), size=(1280, 1280))

angle = 1.0*25


for frame in range(len(fileh.root.Time[:])):
  frame=frame_skip*frame
  # 3D Movie
  
  mlab.clf(figure=fig3d)
  #Simple contour
  scene = mlab.contour3d((fileh.root.Phi[2:-2,2:-2,2:-2,frame]**2), contours=50, extent=ext, vmax=fmax, vmin=fmin)#, colormap='ylgn')

  # Cutplane
  #E  = gradient(fileh.root.Phi[:,:,:,frame])
  #mlab.quiver3d(E[0], E[1], E[2])
  #src = mlab.pipeline.vector_field(E[2], E[1], E[0])
  #mlab.pipeline.vectors(src, mask_points=20, scale_factor=3.)
  #mlab.pipeline.vector_cut_plane(src, mask_points=1, scale_factor=5)
  #src = mlab.pipeline.vector_field(E[0], E[1], E[2])
  #magnitude = mlab.pipeline.extract_vector_norm(src)
  #flow = mlab.pipeline.streamline(magnitude, seedtype='plane', seed_visible=False, seed_scale=1.0, seed_resolution=100, linetype='line',integration_direction='both')
  #$mlab.axes(color=(0.0,0.0,0.),extent=ext,xlabel='X',ylabel='Y',zlabel='Z')
  #source = mlab.pipeline.scalar_field(fileh.root.Phi[:,:,:,frame]**2)
  #vol = mlab.pipeline.volume(source)
  #a = array(gradient(fileh.root.Phi[:,:,:,frame]))**2
  #mlab.quiver3d(a[0],a[1],a[2])
  #  mlab.view(azimuth=angle%360, distance=128.0)
  mlab.view(azimuth=angle%360, distance=10.0)
  angle += 1.0
  mlab.colorbar(orientation='vertical', label_fmt="%.2e", nb_labels=5)
  mlab.outline()
  mlab.savefig("Phi3D_" + filename + "_" + string.zfill(frame, 4) + ".jpg")#, size=(600,600))
  print "[", string.zfill(frame, 4),"/", len(fileh.root.Time[:]), "]"


#mencoder "mf://*.jpg" -mf fps=10 -ovc lavc -o mymovie.av

from enthought.mayavi import mlab
import tables
import string
import sys
from numpy import *


fileh = tables.openFile(sys.argv[1])

fmax = 0.0
for frame in range(len(fileh.root.Time[:])):
  fmax = max([fmax, fileh.root.Phi[2:-2,2:-2, :, frame].max() ])


X = fileh.root.Values.attrs.x#[2:-2]
Y = X/2
Z = 8000./32. * X



ext = [ min(X), max(X), min(Y), max(Y), 0.0, fmax ]


def EnergyFlux(f,v):
  energyflux = zeros( shape(f[0,:,:,:]))


  for i in range(len(f[0,:,0,0])):
    for j in range(len(f[0,0,:,0])):
      for k in range(len(f[0,0,0,:])):
        # Calculate kinetic energy and norm per particle
        energyflux[i,j,k] = (f[:,i,j,k]* v**2).sum() * (f[:,i,j,k]* v).sum() / f.sum(axis=0)[i,j,k]

  return energyflux




angle = 0.0


fig3d = mlab.figure(bgcolor=(1.0, 1.0, 1.0), fgcolor=(0.0, 0.0, 0.0), size=(1280, 1280))

from enthought.tvtk.util.ctf import PiecewiseFunction

for frame in [290]: #range(len(fileh.root.Time[:])):
  frame = frame
  # 3D Movie
  mlab.clf(figure=fig3d)
  #E  = gradient(fileh.root.Phi[:,:,:,frame])
  #pV = ParallelVelocity(fileh.root.Values[:,:,:,:,frame],fileh.root.Values.attrs.v)
  #mlab.quiver3d(E[0], E[1], E[2])
  #src = mlab.pipeline.vector_field(E[2], E[1], E[0])
  #mlab.pipeline.vectors(src, mask_points=20, scale_factor=3.)
  #mlab.pipeline.vector_cut_plane(src, mask_points=1, scale_factor=5)
  #src = mlab.pipeline.vector_field(pV, E[1], E[0])
  #magnitude = mlab.pipeline.extract_vector_norm(src)
  #flow = mlab.pipeline.streamline(magnitude, seedtype='sphere', seed_visible=False, seed_scale=1.0, seed_resolution=12, linetype='line',integration_direction='both')


  angle += 5.0
  mlab.view(azimuth=angle%360)
  #$mlab.axes(color=(0.0,0.0,0.),extent=ext,xlabel='X',ylabel='Y',zlabel='Z')
  source = mlab.pipeline.scalar_field(EnergyFlux(fileh.root.Values[:,:,:,:,frame],fileh.root.Values.attrs.v))
  # Changing the otf:
  vol = mlab.pipeline.volume(source)
  otf = PiecewiseFunction()
  otf.add_point(0.0, 0.0)
  otf.add_point(0.33, 0.85)
  vol._otf = otf
  vol._volume_property.set_scalar_opacity(otf)



  #a = array(gradient(fileh.root.Phi[:,:,:,frame]))**2
  #mlab.quiver3d(a[0],a[1],a[2])
  mlab.colorbar(orientation='vertical')#, label_fmt="%.2f", nb_labels=3)
  mlab.outline()
  mlab.savefig("EnergyFlux_" + string.zfill(frame, 4) + ".jpg")#, size=(600,600))
  print "[", string.zfill(frame, 4),"/", len(fileh.root.Time[:]), "]"



#mencoder "mf://*.jpg" -mf fps=10 -ovc lavc -o mymovie.av
#from pylab import *
from enthought.mayavi import mlab
import tables
import string
import sys
from numpy import *


fileh = tables.openFile(sys.argv[1])



fmax = 0.0
for frame in range(len(fileh.root.Time[:])):
  fmax = max([fmax, fileh.root.Phi[2:-2,2:-2, :, frame].max() ])


X = fileh.root.Phasespace.attrs.x#[2:-2]
Y = X/2
Z = 8000./32. * X



ext = [ min(X), max(X), min(Y), max(Y), 0.0, fmax ]


def ParallelVelocity(f, v):
  velocity = zeros( shape(f[0,:,:,:]))



  for i in range(len(f[0,:,0,0])):
    for j in range(len(f[0,0,:,0])):
      for k in range(len(f[0,0,0,:])):
        # Calculate kinetic energy and norm per particle
        velocity[i,j,k] = (f[:,i,j,k]* v).sum() / f.sum(axis=0)[i,j,k]

  return velocity



fig3d = mlab.figure(bgcolor=(1.0, 1.0, 1.0), fgcolor=(0.0, 0.0, 0.0), size=(1280, 1280))

angle = 0.0


for frame in range(len(fileh.root.Time[:])):
  # 3D Movie
  frame = 10*frame+207
  mlab.clf(figure=fig3d)
  #Simple contour
  #scene = mlab.contour3d(fileh.root.Phi[:,:,:,frame]**2, contours=50)

  # Cutplane

  #E  = gradient(fileh.root.Phi[:,:,:,frame])
  pV = ParallelVelocity(fileh.root.Phasespace[:,:,:,:,frame],fileh.root.Phasespace.attrs.v[2:-2])
  #mlab.quiver3d(pv[0], pv[1], pv[2])
  #plot(pV)
  src = mlab.pipeline.vector_field(pv[2], pv[1], pv[0])
  #mlab.pipeline.vectors(src, mask_points=20, scale_factor=3.)
  #mlab.pipeline.vector_cut_plane(src, mask_points=1, scale_factor=5)
  #src = mlab.pipeline.vector_field(pV, E[1], E[0])
  #magnitude = mlab.pipeline.extract_vector_norm(src)
  #flow = mlab.pipeline.streamline(magnitude, seedtype='sphere', seed_visible=False, seed_scale=1.0, seed_resolution=12, linetype='line',integration_direction='both')






  angle += 5.0
  mlab.view(azimuth=angle%360)
  #$mlab.axes(color=(0.0,0.0,0.),extent=ext,xlabel='X',ylabel='Y',zlabel='Z')
  #source = mlab.pipeline.scalar_field(fileh.root.Phi[:,:,:,frame]**2)
  #vol = mlab.pipeline.volume(source)
  #a = array(gradient(fileh.root.Phi[:,:,:,frame]))**2
  #mlab.quiver3d(a[0],a[1],a[2])
  #mlab.colorbar(orientation='vertical')#, label_fmt="%.2f", nb_labels=3)
  mlab.outline()
  mlab.savefig("E_" + string.zfill(frame, 4) + ".jpg")#, size=(600,600))
  print "[", string.zfill(frame, 4),"/", len(fileh.root.Time[:]), "]"



#mencoder "mf://*.jpg" -mf fps=10 -ovc lavc -o mymovie.av

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

from enthought.mayavi import mlab
import tables
import string
import sys
from numpy import *


fileh = tables.openFile(sys.argv[1])

fmax = 0.0
for frame in range(len(fileh.root.Time[:])):
  fmax = max([fmax, fileh.root.Phi[2:-2,2:-2, :, frame].max() ])


X = fileh.root.Values.attrs.x#[2:-2]
Y = X/2
Z = 8000./32. * X



ext = [ min(X), max(X), min(Y), max(Y), 0.0, fmax ]


def EnergyFlux(f,v):
  energyflux = zeros( shape(f[0,:,:,:]))


  for i in range(len(f[0,:,0,0])):
    for j in range(len(f[0,0,:,0])):
      for k in range(len(f[0,0,0,:])):
        # Calculate kinetic energy and norm per particle
        energyflux[i,j,k] = (f[:,i,j,k]* v**2).sum() * (f[:,i,j,k]* v).sum() / f.sum(axis=0)[i,j,k]

  return energyflux




angle = 0.0


fig3d = mlab.figure(bgcolor=(1.0, 1.0, 1.0), fgcolor=(0.0, 0.0, 0.0), size=(1280, 1280))

from enthought.tvtk.util.ctf import PiecewiseFunction

for frame in [290]: #range(len(fileh.root.Time[:])):
  frame = frame
  # 3D Movie
  mlab.clf(figure=fig3d)
  #E  = gradient(fileh.root.Phi[:,:,:,frame])
  #pV = ParallelVelocity(fileh.root.Values[:,:,:,:,frame],fileh.root.Values.attrs.v)
  #mlab.quiver3d(E[0], E[1], E[2])
  #src = mlab.pipeline.vector_field(E[2], E[1], E[0])
  #mlab.pipeline.vectors(src, mask_points=20, scale_factor=3.)
  #mlab.pipeline.vector_cut_plane(src, mask_points=1, scale_factor=5)
  #src = mlab.pipeline.vector_field(pV, E[1], E[0])
  #magnitude = mlab.pipeline.extract_vector_norm(src)
  #flow = mlab.pipeline.streamline(magnitude, seedtype='sphere', seed_visible=False, seed_scale=1.0, seed_resolution=12, linetype='line',integration_direction='both')


  angle += 5.0
  mlab.view(azimuth=angle%360)
  #$mlab.axes(color=(0.0,0.0,0.),extent=ext,xlabel='X',ylabel='Y',zlabel='Z')
  source = mlab.pipeline.scalar_field(EnergyFlux(fileh.root.Values[:,:,:,:,frame],fileh.root.Values.attrs.v))
  # Changing the otf:
  vol = mlab.pipeline.volume(source)
  otf = PiecewiseFunction()
  otf.add_point(0.0, 0.0)
  otf.add_point(0.33, 0.85)
  vol._otf = otf
  vol._volume_property.set_scalar_opacity(otf)



  #a = array(gradient(fileh.root.Phi[:,:,:,frame]))**2
  #mlab.quiver3d(a[0],a[1],a[2])
  mlab.colorbar(orientation='vertical')#, label_fmt="%.2f", nb_labels=3)
  mlab.outline()
  mlab.savefig("EnergyFlux_" + string.zfill(frame, 4) + ".jpg")#, size=(600,600))
  print "[", string.zfill(frame, 4),"/", len(fileh.root.Time[:]), "]"



#mencoder "mf://*.jpg" -mf fps=10 -ovc lavc -o mymovie.av

from enthought.mayavi import mlab
import tables
import string
import sys
from numpy import *
import scipy

fileh = tables.openFile(sys.argv[1])



fmax = 0.0
for frame in range(len(fileh.root.Time[:])):
  fmax = max([fmax, fileh.root.Phi[2:-2,2:-2, :, frame].max() ])


X = fileh.root.Values.attrs.x#[2:-2]
Y = X/2
Z = 8000./32. * X



ext = [ min(X), max(X), min(Y), max(Y), 0.0, fmax ]


def ParallelVelocity(f, v):
  velocity = zeros( shape(f[0,:,:,:]))

  for i in range(len(f[0,:,0,0])):
    for j in range(len(f[0,0,:,0])):
      for k in range(len(f[0,0,0,:])):
        # Calculate kinetic energy and norm per particle
        velocity[i,j,k] = (f[:,i,j,k]* v).sum() / f.sum(axis=0)[i,j,k]

  return velocity



angle = 0.0

rotEnergy = []

for frame in range(len(fileh.root.Time[:])):

  E  = gradient(fileh.root.Phi[:,:,:,frame])
  pV = ParallelVelocity(fileh.root.Values[:,:,:,:,frame],fileh.root.Values.attrs.v)

  omega = rot(pV, E[1], E[0])
  rotEnergy.append(omega**2 * fileh.root.Values[:,:,:,:,frame].sum(axis=0)).sum().sum().sum()

plot(fileh.root.Time[:], rotEnergy)
savefig("RotationalEnergy")




def plotApStreamlines(ntime = -1, nz = 4, fileh5=fileh[0]):


    D0 = getDomain(fileh5)
    X,Y = meshgrid(D['x'],D['y'])

    Ap = fileh5.root.Potentials.Ap[nz,:,:,ntime]
    #dA_dy, dA_dx = gradient(Ap(X,Y))

    
    dA_dy, dA_dx = gradient(Ap)
    
    streamlines.Streamlines(x,y, dA_dy, -dA_dx).plotArrows()
    
    xlim((min(D['X']), max(D['X'])))
    ylim((min(D['Y']), max(D['Y'])))
    show()


