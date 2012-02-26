#from numpy import *
#import subprocess
# set to execute
import os
from math import *
import sys

test = sys.argv[1]

if test == 'ms':
        sub_exec = 'qsub -g vh10140 -q vh10140'
        cpu_range = [16, 32, 48, 64]
	output = 'scale_ms.txt'
        input  = 'ETG5_sc.helios' 

elif test == 'ss':
	sub_exec = 'qsub -q eh'
	cpu_range = [1,2,4,8,12,16]
	output = 'scale_ss.txt'
        input  = 'ETG5_sc.helios' 

elif test == 'ls':
	sub_exec = 'qsub -q lh10205 -g lh10205'
	cpu_range = [512, 1024]
	output = 'scale_ls.txt'
        input  = 'ETG5_ls.helios' 

else:
	throw()


script = 'heliosThinScaling.sh'
prefix = '../../src/helios -c /home/b/b30683/helios-0.8.30_2/scripts/scaling/'




m_list  = [1, 2, 4, 8]
v_list  = [1, 2]
s_list  = [1]
x_list =  [0.25, 0.5, 1,2, 4, 8, 16]
z_list =  [0.25, 0.5, 1,2, 4, 8, 16]

decomposition = []

for cpus in cpu_range:
    for m in m_list:
        for s in s_list:
            for v in v_list:
                for z in z_list:
                  for x in x_list:
                    rest  = max(1,cpus / (m * s * v))
		    rest2 = int(2**ceil(log(rest)/log(2.) * 0.5))
 	            y = int(max(1, rest2 / x))
                    x = int(max(1, rest  / y))

                    # check for reasoanbles decompositions
                    if(x*y*v*m*s != cpus) : continue
                    if(modf(log(x)/log(2.))[0]  != 0.) : continue
                    if(modf(log(y)/log(2.))[0]  != 0.) : continue


                    decomposition.append('%i:%i:1:%i:%i:%i' % (x,y,v,m,s))

# remove duplicate values
#decomposition = list(set(decomposition))
def f2(seq):  
    # order preserving 
    checked = [] 
    for e in seq: 
        if e not in checked: 
            checked.append(e) 
    return checked

decomposition = f2(decomposition)

# check various decomposition stategies


# cpu decomposition






### Medium Scale Simulations 6D ###

#decomposition = [ '1', '4', '2:2', '2:2:1:2' ]

print decomposition

for i in range(len(decomposition)):

    # replace file and fire.
    run = decomposition[i]

    file = open(script)
    filestr = file.read()
    file.close()
   
    cpu_number = 1
    for k in decomposition[i].split(':'): 
        cpu_number = cpu_number * int(k)
    filestr = filestr.replace("$cpu_num", str(cpu_number))
    filestr = filestr.replace("$decomp", decomposition[i])
    filestr = filestr.replace("$scale_output", output)
    filestr = filestr.replace("$input", prefix+input)
        
    helios_cfg_filename= script + '_%i' % i + '.sh'
    file = open(helios_cfg_filename, "w")
    file.write(filestr)
    file.close()

    #os.system(sub_exec + ' -lP %i' % (max(cpu_number/16,1)) + ' -lp %i ' % (min(cpu_number,16)) + helios_cfg_filename)
    print(sub_exec + ' -lP %i' % (max(cpu_number/16,1)) + ' -lp %i ' % (min(cpu_number,16)) + helios_cfg_filename)
    os.system('rm -f ' + helios_cfg_filename)
    
