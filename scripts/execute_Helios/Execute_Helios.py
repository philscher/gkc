"""
 File        : execute_Helios.py
 Author      : Paul P. Hilscher (2013)
 License     : GPLv3+
 Description : Script to perform parameter scans on
               helios (SLURM)

"""

import subprocess

sbatch_script_raw = \
"""#!/bin/bash
#SBATCH --job-name="gkc++" 
#SBATCH --account=$ACCOUNT 
#SBATCH --nodes=$NUM_NODES
#SBATCH --ntasks-per-node=$NUM_TASK_PER_NODE
#SBATCH --cpus-per-task=$NUM_CPU_PER_TASK
#SBATCH --partition=ALL
#SBATCH --time=1-00:00:00
#SBATCH --output=$JOB_NAME.stdout 
#SBATCH --error=$JOB_NAME.stderr 
#======START=============================== 

export OMP_STACKSIZE=1024M
mpirun -n $NUM_TASK $GKC_EXEC -c $GKC_SCRIPT -f -d $GKC_DECOMP $GKC_PARAM
"""


def execute(gkc_script, param, decomp, job_name, account, helios_exec="../../src/helios", extra_param="", system_shared=True, memory=2):

  # get Queue
  slurm_cmd = "srun --job-name=gkc++ -A %s" % account 
      
  #### Set decomposition options ########
  num_task, num_thread = 1, 1
  for n, k in enumerate(decomp.split(':')): 
    if n != 1 : num_task   = num_task * int(k)
    else      : num_thread = int(k)
 
  num_node = (num_thread*num_task)/16

  # Set batch options
  sbatch_script = sbatch_script_raw

  sbatch_script = sbatch_script.replace("$ACCOUNT", account)
  
  sbatch_script = sbatch_script.replace("$NUM_TASK_PER_NODE", str(16/num_thread))
  sbatch_script = sbatch_script.replace("$NUM_CPU_PER_TASK" , str(num_thread))
  sbatch_script = sbatch_script.replace("$NUM_NODES"        , str(num_node))
  sbatch_script = sbatch_script.replace("$NUM_TASK"         , str(num_task))
  
  sbatch_script = sbatch_script.replace("$JOB_NAME" , job_name)

  # Set gkc commands
  sbatch_script = sbatch_script.replace("$GKC_EXEC"  , "/csc/home3/philscher/gkc-0.20/src/gkc")
  sbatch_script = sbatch_script.replace("$GKC_SCRIPT", gkc_script)
  sbatch_script = sbatch_script.replace("$GKC_DECOMP", decomp)
  sbatch_script = sbatch_script.replace("$GKC_PARAM" , "" if param == "" else " -o \"%s\"" % param)
 

  # Write temporary file & execute & delete
  sbatch_script_name = "%s.sbatch" % job_name
  file_sb = open(sbatch_script_name, 'w')
  file_sb.write(sbatch_script)
  file_sb.close()

  subprocess.call("sbatch " + sbatch_script_name, shell=True)
  
  # delete temporary file is done
  subprocess.call(["rm -rf", sbatch_script_name], shell=True)

  print "Job " + job_name + " executed !" 

   
