import string
import sys
import time
import os
import copy
import subprocess
from cStringIO import StringIO

qsub_gkc_script = \
"""
#! /bin/bash
. /usr/Modules/3.2.9/init/bash
#QSUB -W 840:00 
module load intel/12.1

mpiexec.hydra -n $cpu_num -machinefile ${QSUB_NODEINF} $gkc_exec -f -c $gkc_script -o $helios_opt -d $decomp -x $extra_param 

"""

def executeThin(gkc_script, parameters, extra_param="", **kwargs):


      gkc_exec      = kwargs.pop('gkc_exec', '$HOME/gkc/src/gkc')
      system_shared = kwargs.pop('system_shared'  , True)
      queue         = kwargs.pop('queue'  , 'lh10250')
      decomposition = kwargs.pop('decomposition', '1:1:1:1:1:1')
      extra_param   = kwargs.pop('extra_param', '') 
      memory        = kwargs.pop('memory'  , 2)


      # get Queue
      qsub_cmd = ""
      if   queue == "gr10140b" : qsub_cmd = 'qsub -q gr10140b'
      elif queue == "lr10140b" : qsub_cmd = 'qsub -q lr10140b'
      elif queue == "tb"       : qsub_cmd = 'qsub -q tb'
      elif queue == "eb"       : qsub_cmd = 'qsub -q eb'
      else                     : print "No such queue"
      # get decomposition
      decomp = []
      cpu_number = 1
      thread_number = 1
      print extra_param
      for n, k in enumerate(decomposition.split(':')): 
        decomp.append(int(k)) 
        if n != 1 : cpu_number    = cpu_number * int(k)
        else      : thread_number = int(k)
      # set options and execute 
      
      
      for i in range(len(parameters)):
        filestr = qsub_gkc_script

        filestr = filestr.replace("$gkc_exec", gkc_exec)
        filestr = filestr.replace("$cpu_num", str(cpu_number))
        filestr = filestr.replace("$decomp", decomposition)
        filestr = filestr.replace("$gkc_script", gkc_script)
        filestr = filestr.replace("$helios_opt", parameters[i])
        filestr = filestr.replace("$extra_param", extra_param)
        
        if system_shared == True : filestr = filestr.replace("$system_shared", '@$-oi') 
        else                     : filestr = filestr.replace("$system_shared", '#') 
   
        extra_id = ""
        if extra_param != "" : extra_id = " -x "

        if queue == "tb" :
          # tb is interactive so we should execute process
          qsub_go = "tssrun -c 20:00 -W 20:00 -A p=%i /opt/app/intel/impi/4.0.3.008/bin64/mpiexec.hydra /home/b/b30683/gkc/src/gkc -c %s -f -d %s " \
                      % (cpu_number, gkc_script, decomposition)
          if parameters[i] != "" : qsub_go = qsub_go + " -o \"%s\" " % parameters[i]
          subprocess.call(qsub_go, shell=True)
        else : 
          qsub_go = qsub_cmd + " -W 336:00 -A p=%i:t=%i:c=%i:m=%iG /opt/app/intel/impi/4.0.3.008/bin64/mpiexec.hydra %s -c %s -f -d %s " \
                      % (cpu_number, thread_number, thread_number, memory, gkc_exec, gkc_script, decomposition)
          if parameters[i] != "" : qsub_go = qsub_go + " -o \"%s\" " % parameters[i]
          os.system(qsub_go)
        #        os.system(qsub_cmd + " -W 336:00 -A p=%i:t=1:c=1:m=%iG /opt/app/intel/impi/4.0.3.008/bin64/mpiexec.hydra /home/b/b30683/gkc/src/gkc -c %s -f -d %s -o \"%s\" %s \"%s\"" % (cpu_number, memory, gkc_script, decomposition, parameters[i], extra_id, extra_param))
        print "Job Executed" 

   


