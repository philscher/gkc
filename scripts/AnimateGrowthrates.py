from scipy.io import read_array
from pylab import *
from tables import *

import string

# open HDF5 File

fileh = openFile(sys.argv[1])


files = []
for file in sys.argv[1:]:
    files.append(openFile(file))




def contourXY(time=0, Z = 8, fileh5=fileh):
  contourf(fileh5.root.Phi[Z,:,:,time])
  xlabel("Y")
  ylabel("X")

fig = figure()#frameon=False)
ax = fig.add_axes((0.35,0.2,0.55,0.7))
#DefaultSize = gcf().get_size_inches()
#gcf().set_size_inches( (DefaultSize[0]*2.2, DefaultSize[1]) )

critmin = 0


#fig = plt.figure()
#ax1 = fig.add_subplot(111)

print time, "Len ", len(time)


# Iterate through all pictures
for frame in range(len(time)):
  frame=40*frame+2
  clf()

  # create line
  pl1 = semilogy(time[:frame], ek_mode0[:frame], 'y')

  pl2 = semilogy(time[:frame], ek_mode1[:frame], 'r')
  pl3 = semilogy(time[:frame], ek_mode2[:frame], 'b')
  #twinx()
  pl4 = semilogy(time[:frame], ek_mode3[:frame], 'g')
  pl5 = semilogy(time[:frame], ek_mode4[:frame], 'm')
  pl6 = semilogy(time[:frame], ek_mode5[:frame], 'c')
  xlim((time[0], time[-1]))
  xlabel("Time [$\omega_0^{-1}$]")
  ylabel("Energy in mode [$\\log \\left[|E_k^2|\\right]$")
  #t = text(0.1,0.8,'T = %2.f' % time[frame], size='x-large')
  legend( (pl1, pl2, pl3, pl4, pl5, pl6) , ("m = 0", "m = 1", "m = 2", "m = 3", "m = 4", "m = 5"), loc='upper left')
  #title("Time evolution of first five modes : T = %2.f" % time[frame])

  twiny()
  # create big dot
  semilogy(timestep[frame], ek_mode0[frame], 'y.')
  semilogy(timestep[frame], ek_mode1[frame], 'r.')
  semilogy(timestep[frame], ek_mode2[frame], 'b.')
  semilogy(timestep[frame], ek_mode3[frame], 'g.')
  semilogy(timestep[frame], ek_mode4[frame], 'm.')
  semilogy(timestep[frame], ek_mode5[frame], 'c.')
  xlim((timestep[0], timestep[-1]))

  xlabel("Timestep")

  print "Time   ", timestep[frame], "    ", time[frame]

  # Frame scrolling
  ymin = 1.0e99
  ymax = 0.0
  if frame < 202:
    for a in modes:
      ymin = min(ymin, a[2:frame+1].min())
      ymax = max(ymax, a[2:frame+1].max())

  elif ((frame > 202) and (critmin  == 0)):
    for a in modes:
      ymin = min(ymin, a[frame-200:frame+200].min())
      ymax = max(ymax, a[frame-200:frame+200].max())
      if ymin*0.1 > 1.0e-6: critmin = 1
  else:
    ymin = 1.0e-6
    ymax = 1.0e3
  ylim((0.5*ymin,5*ymax))



  savefig("Growthrate_" + string.zfill(frame, 6) + ".png")#, size=(600,600))
  print "[", string.zfill(frame, 6),"/", len(fileh.root.Values[0,0,:]), "]"

  print ek_mode0[frame],  ek_mode1[frame],  ek_mode2[frame],  ek_mode3[frame],  ek_mode4[frame],  ek_mode5[frame]


