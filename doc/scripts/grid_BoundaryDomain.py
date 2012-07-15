from pylab import *
import random
import matplotlib.lines as lines
import matplotlib.patches as patches
#import matplotlib.text as text
import matplotlib.collections as collections
import matplotlib.units as units
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.collections import PatchCollection
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import pylab 

# Add MPI communication and information exchange

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_frame_on(False)
ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)

#ax.fax1 = plt.axes(frameon=False)

def plotMeasure(start_xy, eo_xy, text):
  l_y = linspace(start_xy[1], eo_xy[1]+3., 128)
  plot(ones(128)*start_xy[0], l_y, 'k--')
  plot(ones(128)*eo_xy[0]   , l_y, 'k--')
  
  
  l_x = linspace(start_xy[0], eo_xy[0], 128) 
  plot(l_x, ones(128)*eo_xy[1], 'k-')
  
  pylab.text(0.5*(start_xy[0]+eo_xy[0]), eo_xy[1]+3., text, ha='center')


 # Second Domain
rect_1 = patches.Rectangle( (-50, 0.), width=50, height=20, alpha=0.2, axes=ax, color='b')
ax.add_patch(rect_1)
plotMeasure((-50., 20.), (0., 30.), "RxLD")


rect_1 = mpatches.Rectangle( (-60., 0.), width=70, height=20, alpha=0.2, axes=ax)
rect_1.set_fc('g')
ax.add_patch(rect_1)
plotMeasure((-60., 20.), (10., 45.), "RxLB")




 # First Domain
rect_1 = patches.Rectangle( (-100, -20), width=50, height=20, alpha=0.2, axes=ax, color='b')
ax.add_patch(rect_1)
plotMeasure((0., -20.), (50., -35.), "RxLD")


rect_1 = mpatches.Rectangle( (-110., -20), width=70, height=20, alpha=0.2, axes=ax)
rect_1.set_fc('g')
ax.add_patch(rect_1)
plotMeasure((-10., -20.), (60., -50.), "RxLB")

# Third
rect_1 = patches.Rectangle( (0., -20), width=50, height=20, alpha=0.2, axes=ax, color='b')
ax.add_patch(rect_1)
rect_1 = mpatches.Rectangle( (-10., -20), width=70, height=20, alpha=0.2, axes=ax)
rect_1.set_fc('g')
ax.add_patch(rect_1)

# Fourth 
rect_1 = patches.Rectangle( (50, 0.), width=50, height=20, alpha=0.2, axes=ax, color='b')
ax.add_patch(rect_1)
rect_1 = mpatches.Rectangle( (40., 0.), width=70, height=20, alpha=0.2, axes=ax)
rect_1.set_fc('g')
ax.add_patch(rect_1)



plotMeasure((-100., 0.), (100., 65.), "RxGD")
plotMeasure((-110., 0.), (110., 85.), "RxGB")

# Manually set legend
rect_1 = patches.Rectangle( (-120., -40), width=12, height=8, alpha=0.2, axes=ax, color='b')
ax.add_patch(rect_1)
pylab.text(-100, -36, "Domain", ha='left', va='center')
rect_1 = mpatches.Rectangle( (-120., -54), width=12, height=8, alpha=0.2, axes=ax)
rect_1.set_fc('g')
ax.add_patch(rect_1)
pylab.text(-100, -50, "Boundary Domain", ha='left', va='center')




#rect_1 = patches.Rectangle( (-100, 20.), width=200, height=20, alpha=0.2, axes=ax, color='b')
#rect_1.set_fc('r')
#ax.add_patch(rect_1)
#rect_1 = mpatches.Rectangle( (-110., 20.), width=220, height=20, alpha=0.2, axes=ax)
#rect_1.set_fc('g')
#ax.add_patch(rect_1)




#rect_2 = mpatches.Rectangle( (10, 10), width=120, height=20, alpha=0.2, axes=ax)


#rect_2 = mpatches.Rectangle( (10, 10), width=120, height=20, alpha=0.2, axes=ax)


#rect_2 = mpatches.Rectangle( (10, 10), width=120, height=20, alpha=0.2, axes=ax)
#rect_2.set_fc('r')
#ax.add_artist(rect_2)
#rect_2 = mpatches.Circle( (0, 0), 0.5, alpha=0.2, axes=ax, color='r')
#rect_2.set_fc('r')
#ax.add_artist(rect_2)


#xlim((-150, 150))
#ylim((-100, 100))
#collection.set_array(np.array(colors))
#ax.add_collection(collection)
savefig("Grid_BoundaryDomain.png", bbox_inches='tight')
