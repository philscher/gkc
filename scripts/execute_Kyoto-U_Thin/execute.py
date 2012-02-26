#import subprocess
import os


pi = 3.141592653589793238462643

  # for local computer + supported decompositions
sub_exec = 'qsub -g lh10205 -q lh10205'
script = 'helios.sh'

parameter_file = 'ETG.helios'


# format is ( #number of cpus,  command line options )

    # get shear paramter
def get_cparameters_list() :

    cp_list = []
    for n in range(16):
        shear =  16./32. * 1./(2. * pi) * n
    
        cp_list.append( ( 1, ' -c ETG_Shear.helios -o "Geometry.shear = %.3f; Output=ETG_Shear_%.3f.h5"' % (shear, shear)))
    
    return cp_list

#"""
#cparameter_list =  [\
#    ( 1, '-d 16:16 -o "L=(64,64,8000,10);N=(128, 256,  64, 128); output=/home/....."'),\
#    ( 1, '-d 16:16 -o "L=(64,64,8000,10);N=(128, 512,  64, 128); output=/home/....."'),\
#    ( 1, '-d 16:16 -o "L=(64,64,8000,10);N=(256, 256,  64, 128); output=/home/....."'),\
#    ( 1, '-d 16:16 -o "L=(64,64,8000,10);N=(128, 256, 128, 128); output=/home/....."'),\
#    ( 1, '-d 16:16 -o "L=(64,64,8000,10);N=(128, 256,  64,  64); output=/home/....."'),\
#    (1, '-d 32:32 -o "L=(64,64,8000,10);N=(256, 512, 128, 128); output=/home/....."'),\
#    ]
#
#"""

cparameter_list = get_cparameters_list()

for cparameters in cparameter_list:
    # open raw job file and replace parameter file and command line parameters
    file = open(script)
    filestr = file.read()
    file.close()
   
    filestr = filestr.replace("$cfgfile",     parameter_file)
    filestr = filestr.replace("$cparameters", cparameters[1])
    filestr = filestr.replace("$cpu_number" , str(cparameters[0]))
        
    # write new job file
    helios_cfg_filename= "execScript_" + script 
    file = open(helios_cfg_filename, "w")
    file.write(filestr)
    file.close()

    cpu_number = cparameters[0]
    # execute job
    print (sub_exec + ' -lP %i' % (max(cpu_number/16,1)) + ' -lp %i ' % (min(cpu_number,16)) + helios_cfg_filename + " param : " + cparameters[1])
    #os.system(sub_exec + ' -lP %i' % (max(cpu_number/16,1)) + ' -lp %i ' % (min(cpu_number,16)) + helios_cfg_filename)
    #os.system('rm -f ' + helios_cfg_filename)

  
