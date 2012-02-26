import string
import sys
import time
import os
import copy


def linspace(start, stop, count):
	x = []
	for v in range(count) : x.append(start+v/float(count-1)*(stop-start))
        return x




qsub_helios_script = \
"""
#! /bin/bash 
# @$-q vh10140 
# @$-g vh10140
$system_shared 
# @$-eo 
# @$-lP 2 
# @$-lp 16 
# @$-lm 28gb 
# @$-llm unlimited 
# @$-Pvn abs_unpack -Pvs unpack -Pvc unpack 
#
set -x

# for openmpi
. /thin/local/etc/setprofile/openmpi-1.3.3+gcc-4.1.sh
. /thin/local/etc/setprofile/intel-11.1.sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/b/b30683/blitz-gcc/lib/

#. /thin/local/etc/setprofile/openmpi-1.3.3+intel-11.0.sh
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/b/b30683/blitz-intel/lib/

cd $QSUB_WORKDIR 

mpiexec -n $cpu_num -machinefile ${QSUB_NODEINF} $helios_exec -f -c $helios_script -o $helios_opt -d $decomp 

"""


def executeThin(helios_script, parameters, decomposition="1:1:1:1:1:1:1", helios_exec="../../src/helios", queue="lh10250", system_shared=True):

      # get Queue
      qsub_cmd = ""
      if   queue == "vh10140" : qsub_cmd = 'qsub -g vh10140 -q vh10140'
      elif queue == "lh10250" : qsub_cmd = 'qsub -g lh10250 -q lh10250'
      elif queue == "lh10247" : qsub_cmd = 'qsub -g lh10247 -q lh10247'
      elif queue == "eh"      : qsub_cmd = 'qsub -g vh10140      -q eh'
      else                    : print "No such queue"
      # get decomposition
      decomp = []
      cpu_number = 1
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
        
        if system_shared == True : filestr = filestr.replace("$system_shared", '@$-oi') 
	else                     : filestr = filestr.replace("$system_shared", '#') 
       
        helios_cfg_filename= "helios" + '_%i' % i + '.sh'
        file = open(helios_cfg_filename, "w")
        file.write(filestr)
        file.close()

        os.system(qsub_cmd + ' -lP %i' % (max(cpu_number/16,1)) + ' -lp %i ' % (min(cpu_number,16)) + helios_cfg_filename)
        print(qsub_cmd + ' -lP %i' % (max(cpu_number/16,1)) + ' -lp %i ' % (min(cpu_number,16)) + helios_cfg_filename)
        os.system('rm -f ' + helios_cfg_filename)
   


