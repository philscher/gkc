import string
import sys
import time
import os
import copy
import subprocess
from cStringIO import StringIO

def linspace(start, stop, count):
	x = []
	for v in range(count) : x.append(start+v/float(count-1)*(stop-start))
        return x




qsub_helios_script = \
"""
#! /bin/bash
. /usr/Modules/3.2.9/init/bash
#QSUB -W 840:00 
module load intel/12.1

mpiexec.hydra -n $cpu_num -machinefile ${QSUB_NODEINF} $helios_exec -f -c $helios_script -o $helios_opt -d $decomp -x $extra_param 

"""


def executeThin(helios_script, parameters, decomposition="1:1:1:1:1:1:1", helios_exec="../../src/helios", queue="lh10250", extra_param="", system_shared=True, memory=2):

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
      print extra_param
      for k in decomposition.split(':'): 
        decomp.append(int(k)) 
        cpu_number = cpu_number * int(k)
      # set options and execute 
      
      for i in range(len(parameters)):
        filestr = qsub_helios_script

        filestr = filestr.replace("$helios_exec", helios_exec)
        filestr = filestr.replace("$cpu_num", str(cpu_number))
        filestr = filestr.replace("$decomp", decomposition)
        filestr = filestr.replace("$helios_script", helios_script)
        filestr = filestr.replace("$helios_opt", parameters[i])
        filestr = filestr.replace("$extra_param", extra_param)
        
        if system_shared == True : filestr = filestr.replace("$system_shared", '@$-oi') 
        else                     : filestr = filestr.replace("$system_shared", '#') 
   
        extra_id = ""
        if extra_param != "" : extra_id = " -x "
        #helios_cfg_filename= "helios" + '_%i' % i + '.sh'
        #file = open(helios_cfg_filename, "w")
        #file.write(filestr)
        #file.close()

        #subprocess.call(qsub_cmd + ' -A p=%i:t=1 ' % (max(cpu_number/16,1)), stdin=StringIO('one\ntwo\nthree\nfour\nfive\nsix\n'))
        #os.system(qsub_cmd + ' -A p=%i:t=1 ' % (max(cpu_number/16,1)) + helios_cfg_filename)
        #print (qsub_cmd + ' -A p=%i:t=1 ' % (max(cpu_number/16,1)) + helios_cfg_filename)
        #os.system(qsub_cmd + " -W 1:00 -A p=%i:t=1:c=1  /opt/app/intel/impi/4.0.3.008/bin64/mpiexec.hydra /home/b/b30683/gkc/src/gkc -c MI_shear_eta_width.helios -f -d %s -o " % (cpu_number, decomposition) + parameters[i])
        print "0------iii------------------> \"%s\" " % parameters[i]
        if queue == "tb" :
          # tb is interactive so we should execute process
          qsub_go = "tssrun -c 20:00 -W 20:00 -A p=%i /opt/app/intel/impi/4.0.3.008/bin64/mpiexec.hydra /home/b/b30683/gkc/src/gkc -c %s -f -d %s -o %s " \
                      % (cpu_number, helios_script, decomposition, parameters[i])
          print "---##############33-->", qsub_go
          subprocess.call(qsub_go, shell=True)
        else : 
          qsub_go = qsub_cmd + " -W 336:00 -A p=%i:t=1:c=1:m=%iG /opt/app/intel/impi/4.0.3.008/bin64/mpiexec.hydra /home/b/b30683/gkc/src/gkc -c %s -f -d %s " \
                      % (cpu_number, memory, helios_script, decomposition)
          if parameters[i] != "" : qsub_go = qsub_go + " -o \"%s\" " % parameters[i]
          os.system(qsub_go)
        #        os.system(qsub_cmd + " -W 336:00 -A p=%i:t=1:c=1:m=%iG /opt/app/intel/impi/4.0.3.008/bin64/mpiexec.hydra /home/b/b30683/gkc/src/gkc -c %s -f -d %s -o \"%s\" %s \"%s\"" % (cpu_number, memory, helios_script, decomposition, parameters[i], extra_id, extra_param))
        print "Job Executed" 

#os.system('rm -f ' + helios_cfg_filename)
   


