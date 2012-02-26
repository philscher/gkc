from sets import Set
#from numpy import *
#import subprocess
# set to execute
import os
from math import *
import sys

test = sys.argv[1]

sub_exec   = 'qsub -g lh10250 -q lh10250'

cpu_range  = [1, 8, 16, 32, 64 ]
cfg_script = sys.argv[1]
script     = 'helios_thin.sh'

os.system('rm -f *.h5')
os.system('rm -f N*')

m_list  = [1 ]
#, 8, 16 ]
v_list  = [1, 2 ]
s_list  = [1]
x_list =  [0.25, 0.5, 1,2, 4, 8 ]
y_list =  [ 1 ]
z_list =  [0.25, 0.5, 1,2, 4, 8  ]

decomposition = Set([])


# check for hdf5 files, need to be clean

for cpus in cpu_range:
    for m in m_list:
        for s in s_list:
            for v in v_list:
                for z in z_list:
                  for x in x_list:
                   for y in y_list:
                     
                     rest  = max(1,cpus / (m * s * v))
                     rest2 = int(2**ceil(log(rest)/log(2.) * 0.5))

                     z = int(max(1, rest2 / x))
                     x = int(max(1, rest  / z))

                     # check for reasoanbles decompositions
                     if(x*y*z*v*m*s != cpus) : continue
                     if(modf(log(x)/log(2.))[0]  != 0.) : continue
		     # check if already exists
                     decomposition.add('%i:%i:%i:%i:%i:%i' % (x,y,z,v,m,s))



decomposition = list(decomposition)

for i in range(len(decomposition)):

    # replace file and fire.
    run = decomposition[i]

    file = open(script)
    filestr = file.read()
    file.close()
   
    cpu_number = 1
    decomp = []
    for k in decomposition[i].split(':'): 
        decomp.append(int(k)) 
        cpu_number = cpu_number * int(k)
    
    filestr = filestr.replace("$cpu_num", str(cpu_number))
    filestr = filestr.replace("$decomp", decomposition[i])
    filestr = filestr.replace("$cfg_file", cfg_script)
    filestr = filestr.replace("$cfg_options", "\"DataOutput.OutputFileName=check_%iX%iY%iZ%iV%iM%iS.h5;\"" % tuple(decomp))
       
      
    helios_cfg_filename= script + '_%i' % i + '.sh'
    file = open(helios_cfg_filename, "w")
    file.write(filestr)
    file.close()

    os.system(sub_exec + ' -lP %i' % (max(cpu_number/16,1)) + ' -lp %i ' % (min(cpu_number,16)) + helios_cfg_filename)
    print(sub_exec + ' -lP %i' % (max(cpu_number/16,1)) + ' -lp %i ' % (min(cpu_number,16)) + helios_cfg_filename)
    os.system('rm -f ' + helios_cfg_filename)
    
