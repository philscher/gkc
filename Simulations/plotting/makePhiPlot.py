
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
