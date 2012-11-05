import os
import string

from pylab import *

from gkcData      import *
from gkcAnalysis  import *
from gkcLinear    import *
from gkcStyle     import *
from gkcNonLinear import *

import tables


setPlotOutputThesis()


fileh = []
D     = []

# Open gkc files which where attached as arguments
# Use ipython !
for fileh5 in sys.argv[1:]:
     try:
        fileh.append(tables.openFile(fileh5))
     except:
        print "Cannot open : " , fileh5, " (skipping)"
        #     D.append(getDomain(fileh5 = fileh[-1]))


# reload data states, (the simple way)
def reload():
   files= []  
   for file in sys.argv[1:]:
         files.append(openFile(file))


