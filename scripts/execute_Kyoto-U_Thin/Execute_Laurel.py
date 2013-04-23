"""
   Script to execute parameter scans for the Kyoto university supercomputer

"""

import string
import os
import subprocess

"""
 wtime = set the max. wall clock time
"""
def executeThin(gkc_script, parameters, decomposition="auto", extra_param="", **kwargs):

  gkc_exec      = kwargs.pop('gkc_exec', '$HOME/gkc/src/gkc')
  system_shared = kwargs.pop('system_shared'  , True)
  queue         = kwargs.pop('queue'  , 'lh10250')
  #decomposition = kwargs.pop('decomposition', '1:1:1:1:1:1')
  extra_param   = kwargs.pop('extra_param', '') 
  memory        = kwargs.pop('memory'  , 2)
  wtime         = kwargs.pop('wtime', 220) 

  # get Queue
  qsub_cmd = ""
  if   queue == "gr10140b" : qsub_cmd = 'qsub -q gr10140b'
  elif queue == "tb"       : qsub_cmd = 'qsub -q tb'
  else                     : print "No such queue"
  
  # get decomposition
  decomp, cpu_number, thread_number = [], 1, 1
  for n, k in enumerate(decomposition.split(':')): 
    decomp.append(int(k)) 
    if n != 1 : cpu_number    = cpu_number * int(k)
    else      : thread_number = int(k)
  
  ################## set options and execute  #################3
  extra_id = ""
  if extra_param != "" : extra_id = " -x "

  if queue == "tb" :

    # tb is interactive so we should execute process
    qsub_go = "tssrun -c 20:00 -W 20:00 -A p=%i /opt/app/intel/impi/4.0.3.008/bin64/mpiexec.hydra /home/b/b30683/gkc/src/gkc -c %s -f -d %s " \
                      % (cpu_number, gkc_script, decomposition)
    if parameters[i] != "" : qsub_go = qsub_go + " -o \"%s\" " % parameters[i]
    subprocess.call(qsub_go, shell=True)
    
  else : 
    qsub_go = qsub_cmd + " -W %i:00 -A p=%i:t=%i:c=%i:m=%iG /opt/app/intel/impi/4.0.3.008/bin64/mpiexec.hydra %s -c %s -f -d %s " \
                      % (wtime, cpu_number, thread_number, thread_number, memory, gkc_exec, gkc_script, decomposition)
    
    qsub_go = qsub_go + " -o \"%s\" " % parameters
    os.system(qsub_go)
    
  print "Job Executed" 

