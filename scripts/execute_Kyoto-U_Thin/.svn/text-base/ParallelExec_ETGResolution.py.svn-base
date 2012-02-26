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
use_nthreads = int(sys.argv[3])
print " use_nthreads = ", use_nthreads 

# define parameters of the parameter grid
helios_filename = sys.argv[2]
prog_path = sys.argv[1]

#parameter_list       = logspace(1.3,2.,6)

parameter_list = array(["(32, 16, 8, 16)", \
                        "(64, 16, 8, 16)", \
                        "(32, 32, 8, 16)", \
                        "(64, 32, 8, 16)", \
                        "(32, 16, 16, 16)",\
                        "(32, 32, 16, 16)",\
                        "(64, 32, 16, 16)",\
                        "(64, 32, 16, 16)",\
                        "(32, 16, 8, 32)", \
                        "(64, 16, 8, 32)", \
                        "(32, 32, 8, 32)", \
                        "(64, 32, 16, 32)"])


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
            print("Running task : eta_i " + task_unit +  " Jobs  Left : " + str(task_queue.qsize()))
        except Queue.Empty:
            # Nothing to do anymore, return
            print("Nothing to do anymore .... exiting")
            return

        ############################################################
        # Here is main program
        # do not run setup in parallel! 
        self.vlock.acquire()
      
        parameter = task_unit

        helios_basename = os.path.splitext(os.path.split(helios_filename)[1])[0]
        helios_paramfile_name = helios_basename + "_run_" +  str(where(parameter_list == parameter)[0][0]) + "_.helios" 
        helios_hdf5_filename = os.getcwd() + "/" + helios_basename + "_" + str(where(parameter_list == parameter)[0][0]) + ".h5"

         #create helios file
        file = open(helios_filename)
        filestr = file.read()
        file.close()
        filestr = filestr.replace("$value", parameter)
        filestr = filestr.replace("$hdf5_filename", helios_hdf5_filename)
        file = open(helios_paramfile_name, "w")
        file.write(filestr)
        file.close()

        self.vlock.release()
        # run process in parallel (no process has to wait):  

        retcode=os.system(prog_path + ' -c ' + helios_paramfile_name)
       
        print " end of run: return code = ", retcode 

        #############################################################3            
        # Data analysis, e.g. Error calculation

        self.vlock.acquire()
        # do not run data analysis in parallel!!! 
        growthrates = linspace(0.0,1.0,6)

        #statfile = open('Growthrates.txt', 'a')
        #line = parameter + " %f %f %f %f %f %f\n" % ( GrowthrateFit.getGrowthrates(helios_hdf5_filename),\
        #        growthrates[1], growthrates[2], growthrates[3], growthrates[4], growthrates[5])
        #                                         
        #                                         
        #statfile.write(line)
        #statfile.close()
  
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
    




