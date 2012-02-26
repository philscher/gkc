from scipy import *

from scipy.io import read_array
from pylab import *
import string
import sys

import time

import threading
import Queue

import os
import subprocess 

import GrowthrateFit
############ set user parameters
# define the number of threads
use_nthreads = 2
print " use_nthreads = ", use_nthreads 

# define parameters of the parameter grid
helios_filename = "ITG_Growthrates.raw"
parameter_list       = logspace(1.3,2.,6)
parameter_name  = '$dvalue1'



# create queue
task_queue = Queue.Queue()

for parameter in parameter_list:
        task_unit = (parameter)
        task_queue.put(task_unit)
        

class ProgramRun(threading.Thread):
    vlock = threading.Lock()

    def __init__(self):
      threading.Thread.__init__(self)

    def run(self):
      
      # Use infinite loop, not Queue.Full() because of threading
      while(True):

        try:
            # Get task unit, non-blocking
            task_unit = task_queue.get(False)
            print("Running task : eta_i %.3f" % task_unit +  " Jobs  Left : " + str(task_queue.qsize()))
        except Queue.Empty:
            # Nothing to do anymore, return
            print("Nothing to do anymore .... exiting")
            return

        ############################################################
        # Here is main program
        # do not run setup in parallel! 
        self.vlock.acquire()
      
        parameter = task_unit

        helios_basename = os.path.splitext(helios_filename)[0]
        helios_paramfile_name = helios_basename + "_run_%.2f" %  (parameter) + "_.helios" 
        helios_hdf5_filename = helios_basename + "_%.2f.h5" % parameter

         #create helios file
        file = open(helios_filename)
        filestr = file.read()
        file.close()
        filestr = filestr.replace(parameter_name, "%.2f" % parameter)
        filestr = filestr.replace("$hdf5_filename", helios_hdf5_filename)
        file = open(helios_paramfile_name, "w")
        file.write(filestr)
        file.close()

        self.vlock.release()
        # run process in parallel (no process has to wait):  

        retcode=os.system('./helios -c ' + helios_paramfile_name)
       
        print " end of run: return code = ", retcode 

        #############################################################3            
        # Data analysis, e.g. Error calculation

        self.vlock.acquire()
        # do not run data analysis in parallel!!! 
        growthrates = linspace(0.0,1.0,6)

        statfile = open('Growthrates.txt', 'a')
        line = "%f %f %f %f %f %f %f\n" % (parameter, GrowthrateFit.getGrowthrates(helios_hdf5_filename),\
                growthrates[1], growthrates[2], growthrates[3], growthrates[4], growthrates[5])
                                                 
                                                 
        statfile.write(line)
  
        self.vlock.release()

# Meister proper
"""
file = open("Growthrates.txt", "w")
file.write("")
file.close()
"""
# Start the processes
for i in range(use_nthreads):
    ProgramRun().start()
    




