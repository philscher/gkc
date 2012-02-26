from Scaling import *

name = "Helios_svn132"
version = "Helios 0.8.6 (svn90)"
title = "Strong Scaling of Vlasov solver $(N_x = 512, N_y = 128, N_z = 64, N_v = 112)$"
text = "Strong Scaling of ETG simulation with $N = \\left( N_x = 512, N_y = 128, N_z=64, N_{v_\parallel}=96 \\right)$"
legend_names = ['Ideal', 'Decomposition in $\\left(x,y,v_\parallel\\right)$']
cpu_list    = []
times_list    = []
decomp_list    = []

P = 0.9985
# we use data from Helios svn. v.40



# we produce 2 plots
# read in scaling results

# For 1 CPUS
a = loadtxt(sys.argv[1], delimiter=":")

for entry in a[:,:] :
  decomp_list.append("(%i:%i:%i:%i:%i:%i:%i)" % (entry[0] , entry[1], entry[2], entry[3], entry[4], entry[5], entry[6]))
  cpu_list   .append(int(entry[0] * entry[1] * entry[2] * entry[3] * entry[4] * entry[5] * entry[6]))
  times_list .append(entry[7])

plotScalingNormDecompLabel(name, cpu_list, times_list, decomp_list, legend_names, text, version, title)

