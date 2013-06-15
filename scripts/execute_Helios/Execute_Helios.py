"""
 File        : execute_Helios.py
 Author      : Paul P. Hilscher (2013)
 License     : GPLv3+
 Description : Script to perform parameter scans on helios (SLURM)

"""

import subprocess

sbatch_script_raw = \
"""#!/bin/bash
#SBATCH --job-name="gkc++" 
#SBATCH --account=$ACCOUNT 
#SBATCH --nodes=$NUM_NODES
#SBATCH --ntasks-per-node=$NUM_TASK_PER_NODE
#SBATCH --cpus-per-task=$NUM_THREADS
#SBATCH --partition=ALL
#SBATCH --time=1-00:00:00
#SBATCH --output=$JOB_NAME.stdout 
#SBATCH --error=$JOB_NAME.stderr 
#======START=============================== 

export OMP_STACKSIZE=1024M
export OMP_NUM_THREADS=$NUM_THREADS
mpirun -n $NUM_TASK $GKC_EXEC -c $GKC_SCRIPT -f -d $GKC_DECOMP $GKC_PARAM
"""

"""
 if notify is sent a mail is send to the user. 
 echo "your_email@provider.com" > $HOME/.forward will forward
 ALL you emails from helios to your account.


"""
def execute(gkc_script, param_list, decomp, job_name, account, **kwargs):

  gkc_exec    = kwargs.pop('gkc_exec', '$HOME/gkc/src/gkc')
  extra_param = kwargs.pop('extra_param', '')
  mail        = kwargs.pop('mail', False)

  if mail == True : mail_mod = ' --mail-type=END '
  else            : mail_mod = ''
 

  #### Bring param into the form which gkc accepts ###
  param = "\"" + param_list[0]
  for p in param_list[1:] : param += "+" + p
  param += "\""
   
 
  #### Set decomposition options ########
  num_task, num_thread = 1, 1
  for n, k in enumerate(decomp.split(':')): 
    if n != 1 : num_task   = num_task * int(k)
    else      : num_thread = int(k)
 
  num_node = (num_thread*num_task)/16

  ############# Set batch options #############
  sbatch_script = sbatch_script_raw

  sbatch_script = sbatch_script.replace("$ACCOUNT", account)
  sbatch_script = sbatch_script.replace("$JOB_NAME"         , job_name)
  
  sbatch_script = sbatch_script.replace("$NUM_TASK_PER_NODE", str(16/num_thread))
  sbatch_script = sbatch_script.replace("$NUM_THREADS"      , str(num_thread))
  sbatch_script = sbatch_script.replace("$NUM_NODES"        , str(num_node))
  sbatch_script = sbatch_script.replace("$NUM_TASK"         , str(num_task))

  ################ Set gkc commands ###################
  sbatch_script = sbatch_script.replace("$GKC_EXEC"  , gkc_exec)
  sbatch_script = sbatch_script.replace("$GKC_SCRIPT", gkc_script)
  sbatch_script = sbatch_script.replace("$GKC_DECOMP", decomp)
  sbatch_script = sbatch_script.replace("$GKC_PARAM" , "" if param == "" else " -o \"%s\"" % param)
 

  ##############  Write temporary file & execute & delete ##########
  sbatch_script_name = "%s.sbatch" % job_name
  file_sb = open(sbatch_script_name, 'w')
  file_sb.write(sbatch_script)
  file_sb.close()

  subprocess.call("sbatch --signal USR1 " + mail_mod  + sbatch_script_name, shell=True)
  
  # delete temporary file is done
  subprocess.call("rm -rf " + sbatch_script_name, shell=True)

  print "Job " + job_name + " executed !" 

   
