#! /bin/bash 
# @$-q vh10140 
# @$-g vh10140 
# @$-oi 
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

#mpiexec -n $cpu_num   -machinefile ${QSUB_NODEINF} gdb --command=helios.gdb --args ../../src/helios -f -c $cfg_file -o $cfg_options -d $decomp 
mpiexec -n $cpu_num   -machinefile ${QSUB_NODEINF} ../../src/helios -f -c $cfg_file -o $cfg_options -d $decomp 
