"""
    Benchmark for collisionless drift Alfven waves.

    Ref. : F. Jenko, Computer Physics Communication 125, (2000) 196-209


    Solved is the Dispersion Relation from Eq. (26)

Take care of nomralization

alpha^2 = \frac{2}{\hat{\epsilon}}
\hat{b} = \beta_e \hat{\epsilon}
\hat{mu} = \mu_e \hat{\epsilon}



"""


from pylab import *
from scipy import *
from scipy.special.orthogonal import *
from PyPPL import *
from scipy.optimize import *
from numpy import *
import os
import mpmath as mp
from Execute import *
from FitFunction import *
from tables import *

S = { 'kx' : 0., 'ky' : 0.5, 'kp' : 1.0, 'mu' : 10, 'eta_e' : 0., 'm_ie' : 1837.}
#beta_list = logspace(-2,1.,6)
beta_list = linspace(1.,10.,4)





def performSimulation(beta):
        
        helios_paramfile_name = "DriftAlfvenWave.helios" 
        #helios_options = "\"Plasma.Beta=%f;DataOutput.OutputFileName=DriftAlfvenWave_%f.h5;\"" % (beta*S['mu']/S['m_ie'], beta)
        helios_options = "\"Plasma.Beta=%f;DataOutput.OutputFileName=DriftAlfvenWave_%f.h5;\"" % (beta, beta)

        print ('./helios -c ' + helios_paramfile_name + " -o " + helios_options)
        retcode=os.system('../../src/helios -f -c ' + helios_paramfile_name + " -o " + helios_options + " -f ")

        # Now start Analysis
        helios_benchfile = "DriftAlfvenWave_%f.h5" % beta
        
        fileh5 = openFile(helios_benchfile)
        """
        phi2 = fileh5.root.Analysis.PowerSpectrum.Y[0,1,:]
        def getTime(timing):
            timing_new = []
            for t in timing: timing_new.append((t[0], t[1]))
            return array(timing_new)

        T    = getTime(fileh5.root.Analysis.PowerSpectrum.Timing[:])
        """
        Y = fileh5.root.Analysis.scalarValues.cols.phiEnergy[:]
        T = fileh5.root.Analysis.scalarValues.cols.Time[:]
        #os.remove(helios_benchfile)
   
        # Plots snapshot
        clf()
        plot(T, Y)
        savefig("DriftAlfven_Beta_%f.png" % beta)



        #fit 
        gamma = fitExpGrowthOptimize(T, Y, len(T)/2)
        os.remove(helios_benchfile)
        return beta, gamma



##############################################  Dispersion Relation ###########################

def DriftAlfvenWave(w, S):
  
  kp2 = S['kx']**2 + S['ky']**2
  mu = S['mu']
  beta = S['beta'] 
  
  alpha = sqrt(2./mu)

  kp  = S['kp']
  ws = S['ky']/(alpha*kp)
  wh = w / (alpha*kp)
  eta_e = S['eta_e']


  def Y(w) : return w + (w**2 - 0.5) * Z(w)

  return (kp2 *w + ((w-ws)-eta_e * w* ws * Y(w)) * ( 1. + w * Z(w)) * ( 1. - 2. * w**2 * beta/mu))

########################################### 






# Solve theoretical dispersion relation
theo_solution = []
for beta in linspace(min(beta_list), max(beta_list)+1.,301):
#for beta in logspace(log10(min(beta_list)), log10(max(beta_list)),301):
  S['beta'] = beta
  theo_solution.append((beta, sqrt(2./S['mu'])*S['kp']*getComplexZero(DriftAlfvenWave, S)))
theo_solution = array(theo_solution)


##########################################
# Solve numerically
num_solution = []
      

num_solution = ManyExecture(performSimulation, beta_list, 2)






#for beta in beta_list: 
#  w, w_err  =  performSimulation(beta)
#  num_solution.append((beta, w, w_err))
num_solution = array(num_solution)


print "Result :: " , num_solution
# Start the processes
#for i in range(use_nthreads): ProgramRun().start()
clf()
#plt.subplot(111, xscale="log", yscale="log") 
#plot(theo_solution[:,0], abs(real(theo_solution[:,1])), 'r')
plot(theo_solution[:,0], imag(theo_solution[:,1]), 'b')
plot( num_solution[:,0], num_solution[:,1]       ,'b^', label="LeastSq", markersize=14., markerfacecolor='None', markeredgecolor='b', markeredgewidth=3.)
xlim((min(beta_list), max(beta_list)))
xlabel("$\\hat{\\beta}$")
ylabel("$\\omega$")

savefig("DriftAlfvenDamping.png")






