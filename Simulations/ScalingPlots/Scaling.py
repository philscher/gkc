from  scipy import *
import pylab
from pylab import *
from numpy import *
import matplotlib.pyplot as plt
       

markers = ['v-c', '^-y', '<-r', '>-m', '*-b', 'd-g', 'p-b', '1-r', '2-m']

# set some basic parameters       
def log_label_ticks(x, pos):
          #"""The two args are the value and tick position.
          # Label ticks with the product of the exponentiation"""
          # see http://bytes.com/forum/thread462268.html
          return '%1i' % (x)

def setPlot(scale=1.0):
  clf()
  # Create frame for plots
  rcParams['legend.loc'] = 'best'
  #rect=[0.1, 0.1, 0.85, 0.85] # = [left, bottom, width, height]
  rect=[0.1, 0.1, 0.85, 0.80] # = [left, bottom, width, height]
  ax=axes(rect, axisbg='w')

  # increase font size of axis
  #params = {'backend': 'ps',
  params = { 
        #     'axes.labelsize' : 50,
             'text.fontsize': 12,
             'legend.fontsize': 12,
              'axes.labelsize': 'large',
             'xtick.labelsize': 26,
             'ytick.labelsize': 26,
             'text.usetex': True,
              'marker.size' : 8,
             }

   #Set size
  DefaultSize = gcf().get_size_inches()
  gcf().set_size_inches( (DefaultSize[0]*scale, DefaultSize[1]) )

def plotTail(cpu_min, cpu_max, text, version="", ylimits="", ylog=False):
  ax=gca()
  xformatter = FuncFormatter(log_label_ticks)
  ax.xaxis.set_major_formatter(xformatter)
  if(ylog==True) : ax.yaxis.set_major_formatter(xformatter)
  xlim((cpu_min*0.9,cpu_max*1.1))
  if ylimits == "" : ylim((cpu_min*0.9,cpu_max*1.1))
  else : ylim(ylimits)
  
  leg = gca().get_legend()
  leg.draw_frame(False)
  #if text != "" :  pylab.text(0.5,0.5, text, fontsize=12)
  title(text)
  import datetime
  now = datetime.datetime.now()
  print version
  if version != "" :  pylab.text(cpu_max,cpu_min/1.2, version + " (" +str( now.strftime("%Y-%m-%d"))+")", fontsize=8)


###########################################################################
# Standard Scaling Plot
def plotScaling(name, cpus_list, times_list, legend_names, text, version = "", P=0.0):
  # get highest cpu number
  cpu_max = 1
  cpu_min = 99999999
  for cpus in cpus_list:
    cpu_min = min(min(cpus), cpu_min)
    cpu_max = max(max(cpus), cpu_max)

  setPlot(scale=1.4142, ylog=True) 

  ######################### Ploting ############################
  speedup_ideal = logspace(max(cpu_min-1, 0), log2(cpu_max**1.01),101, base=2.)
  print speedup_ideal
  pylab.loglog(speedup_ideal, speedup_ideal, 'k--', basex=2.0000001)

  # Amdahl law
  if P != 0.0:
        cpu_theo     = logspace(1, log2(cpu_max)+1,101, base=2.0)
        speedup_theo = 1./(1-P + P/cpu_theo)
        pylab.loglog(cpu_theo, speedup_theo, 'g-', basex=2.0000001, basey=2.000000001)

  
  for i in range(len(cpus_list)):
       times = times_list[i]
       cpus  = cpus_list[i]

       speedup = 1./times*times[0]*cpus[0]

       pylab.loglog(cpus,speedup, markers[i], basex=2.0000001, basey=2.000000001)
  #################################################################
  pylab.legend(legend_names)
  
  plotTail(cpu_min, cpu_max*1.8, text, version)
  

  pylab.ylabel('Speedup')
  pylab.xlabel('Number of CPUs')

  savefig(name + ".png", dpi=300)
  savefig(name + ".eps")


###########################################################################
# Normalized Scaling Plot
def plotScalingNormalized(name, cpus_list, times_list, legend_names, text, version, title, P=0.0):
  # get highest cpu number
  cpu_max = 1
  cpu_min = 99999999
  for cpus in cpus_list:
    cpu_min = min(min(cpus), cpu_min)
    cpu_max = max(max(cpus), cpu_max)

  time0_min = 9999999
  time0_max = 0
  for times in times_list:
    time0_max = max(times[0], time0_max)
    time0_min = min(times[0], time0_min)
    
  setPlot(scale=1.4142) 


  ######################### Ploting
  speedup_ideal = logspace(max(log2(cpu_min)-1, 0), log2(cpu_max**1.01),101, base=2.)
  speedup_ideal = logspace(4, 7,10, base=2.)
  print speedup_ideal
  pylab.semilogx(speedup_ideal, speedup_ideal/speedup_ideal, 'k--', basex=2.0000001)

  # Amdahl law
  #if P != 0.0:
  #      cpu_theo     = logspace(log2(cpu_min), log2(cpu_max)+1,101, base=2.0)
  #      speedup_theo = 1./(1-P + P/cpu_theo)
  #      pylab.semilogx(cpu_theo, cpu_theo, 'g-', basex=2.0000001)

  
  for i in range(len(cpus_list)):
       times = times_list[i]
       cpus  = cpus_list[i]
        

       speedup = 1./((times/time0_min) * (cpus/cpu_min))

       pylab.semilogx(cpus,speedup, markers[i], basex=2.0000001)
  ###################################
  pylab.legend(legend_names, ncol=2)
  pylab.ylabel('$\\frac{1}{\\rm{Time} \\cdot \\rm{CPUs}}$')
  pylab.xlabel('Number of CPUs')
  
  plotTail(cpu_min, cpu_max,  text, version, (0.0,1.1))


  #pylab.rcParams.update(params)
    
  savefig(name + ".png", dpi=300)
  savefig(name + ".eps")


###########################################################################
# Normalized Scaling Plot
def plotScalingNormDecompLabel(name, cpus_list, times_list, decomp_list, legend_names, text, version, title, P=0.0):
  # get highest cpu number
  cpu_max = 1
  cpu_min = 99999999
    
  cpu_min = min(cpus_list)
  cpu_max = max(cpus_list)

  time0 = min(array(times_list)[where(array(cpus_list) == cpu_min)[0]])
  cpu0  = cpu_min
   
  

  setPlot(scale=1.4142) 


  ######################### Ploting
  speedup_ideal = logspace(max(log2(cpu_min)-1, 0), log2(cpu_max**1.2),101, base=2.)
  pl = []
#  speedup_ideal = logspace(4, 7,10, base=2.)
  print speedup_ideal
  pl.append(pylab.semilogx(speedup_ideal, speedup_ideal/speedup_ideal, 'k--', basex=2.0000001))

  # Amdahl law
  #if P != 0.0:
  #      cpu_theo     = logspace(log2(cpu_min), log2(cpu_max)+1,101, base=2.0)
  #      speedup_theo = 1./(1-P + P/cpu_theo)
  #      pylab.semilogx(cpu_theo, cpu_theo, 'g-', basex=2.0000001)

  for i in arange(1,len(cpus_list)):
       times = float(times_list[i])
       cpus  = float(cpus_list[i])
        
       print i, time0, cpu0, times, cpus
       speedup = 1./((times/time0) * (cpus/cpu0))
       pl.append(pylab.semilogx(cpus_list[i],speedup, markers[0], basex=2.0000001))
       pylab.text(cpus_list[i]*(1.05), speedup-0.01, decomp_list[i])

  # get the best performance per CPU number
  cpus_s = list(set(cpus_list))
  cpus_s.sort()
  time_b = []
  for cpu in cpus_s:
    time_bs =  min(array(times_list)[where(array(cpus_list) == cpu)[0]]) 
    time_b.append(1./((time_bs/time0) * (cpu/cpu0)))

  print cpus_s
  print time_b
  pylab.semilogx(cpus_s,time_b, 'r-', basex=2.0000001, linewidth=5.0)


  #pylab.semilogx(linspace(min(cpus_list), max(cpus_list), 11), linspace(0.9,0.9,11))
  # 0.9 % line
  #pylab.semilogx(linspace(min(cpus_list)/2, max(cpus_list)*1.5, 11), linspace(0.9,0.9,11),'k-', basex=2.0000001)
  ###################################
#  pylab.legend(pl[:1], legend_names, ncol=2)
  pylab.legend(legend_names, ncol=2)
  pylab.ylabel('Parallelization Efficiency $\\left[\\frac{1}{\\rm{Time} \\cdot \\rm{CPUs}}\\right]$')
  pylab.xlabel('Number of CPUs')
  
  plotTail(min(cpus_list[1:]), cpu_max*1.5,  text, version, ylimits=(0.0,1.1))

  #pylab.rcParams.update(params)
  gca().yaxis.grid(True)    
  #setp(gca().get_ygridlines(), visible=False)
  savefig(name + ".png", dpi=300)
  savefig(name + ".eps")


