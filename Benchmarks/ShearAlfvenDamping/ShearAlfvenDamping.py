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



# use settings from Lin & Chen, A fluid hybrid electrom model for electromagnetic simulations (2001), PoP Letters 
S  = { "kp" : 0.01    , "ky" : 0.4, "kx" : 0., "T_e" : 1., "T_i"  : 1., "m_ie" : 1837. , 'Ion': 'Adiabatic'}
beta_list = logspace(-5,-2,6)

def performSimulation(beta):
        helios_paramfile_name = "ShearAlfvenDamping.helios" 
        helios_options = "\"Plasma.Beta=%f;Benchmark.FileName=ShearAlfven_%f.txt;Helios.MaxTime=%f;\"" % (beta, beta, 70.+(log10(beta)+5)**(1.5)*70)
       # helios_options = "\"Plasma.Beta=%f;Benchmark.FileName=ShearAlfven_%f.txt;\"" % (beta, beta)#.e4+(log10(beta)+5)**(0.75)*1.e4)

        print ('./helios -c ' + helios_paramfile_name + " -o " + helios_options)
        retcode=os.system('../../src/helios -f -c ' + helios_paramfile_name + " -o " + helios_options + " -f ")

        # Now start Analysis
        helios_benchfile = "ShearAlfven_%f.txt" % beta
        A = loadtxt(helios_benchfile)
        os.remove(helios_benchfile)
   
        # Plots snapshot
        clf()
        plot(A[:,0], A[:,3])
        savefig("ShearAlfven_Beta_%f.png" % beta)

        #fit 
        w, g   = fitDampedOscillationOptmize(A[:,0], A[:,3])
        w2, g2 = fitDampedOscillationFourier(A[:,0], A[:,3])
        return beta, complex(w,g), complex(w2, g2)



##############################################  Dispersion Relation ###########################

def ShearAlfvenWave(w, S):
    kp = S['kp']
    k_perp2 = S['kx']**2. + S['ky']**2.
    T_i  = S['T_i']
    T_e  = S['T_e']
    m_ie = S['m_ie']
    beta = S['beta'] 
    #basis
    m_e = 1./m_ie
    m_i = 1.

    tau = T_i / T_e 
    cs2 =  T_e / m_i
    v_A = sqrt(2. * T_e/(m_i * beta))
    rho_s2 = 1. #T_e * m_ie 

    v_thi = sqrt(T_i/m_i)
    v_the = sqrt(T_e/m_e)

    zeta_i = w / (sqrt(2.) * kp * v_thi)
    zeta_e = w / (sqrt(2.) * kp * v_the)

    Z_zeta_i = Z(zeta_i) if zeta_i < 20. else CDF(zeta_i)
    Z_zeta_e = Z(zeta_e) if zeta_e < 20. else CDF(zeta_e)
      
    if(S['Ion'] == 'Adiabatic') : return (w**2/(kp**2*v_A**2) - 1.) * (1. + zeta_e*Z_zeta_e  )- k_perp2 * rho_s2
    if(S['Ion'] == 'Kinetic'  ) : return (w**2/(kp**2*v_A**2) - 1.) * (1. + zeta_e*Z_zeta_e + tau + tau*zeta_i*Z_zeta_i )- k_perp2 * rho_s2
        
        
    return   (w**2/(kp**2*v_A**2) - 1.) * (1. + zeta_e*Z_zeta_e + tau + tau*zeta_i*Z_zeta_i )- k_perp2 * rho_s2
   
    #return w, S['kp']*v_A*sqrt(1.+k_perp2), S['kp']*sqrt(cs2) 


########################################### 






# Solve theoretical dispersion relation
theo_solution = []
for beta in logspace(log10(min(beta_list))-1, log10(max(beta_list))+1,301):
  S['beta'] = beta
  theo_solution.append((beta, getComplexZero(ShearAlfvenWave, S)))
theo_solution = array(theo_solution)


##########################################
# Solve numerically
num_solution = []
      

num_solution = ManyExecture(performSimulation, beta_list, 3)






#for beta in beta_list: 
#  w, w_err  =  performSimulation(beta)
#  num_solution.append((beta, w, w_err))
num_solution = array(num_solution)


print "Result :: " , num_solution
# Start the processes
#for i in range(use_nthreads): ProgramRun().start()
clf()
#plt.subplot(111, xscale="log", yscale="log") 
loglog(theo_solution[:,0], abs(real(theo_solution[:,1])), 'r')
loglog(theo_solution[:,0], abs(imag(theo_solution[:,1])), 'b')
#loglog( num_solution[:,0], abs(real( num_solution[:,1])),'rs', yerr=real(num_solution[:,2]), markersize=12., markerfacecolor='None')
#loglog( num_solution[:,0], abs(imag( num_solution[:,1])),'bs', yerr=imag(num_solution[:,2]), markersize=12., markerfacecolor='None')
loglog( num_solution[:,0], abs(real( num_solution[:,1])),'rs', label="LeastSq", markersize=14., markerfacecolor='None', markeredgecolor='r', markeredgewidth=2.)
loglog( num_solution[:,0], abs(imag( num_solution[:,1])),'bs',                  markersize=14., markerfacecolor='None', markeredgecolor='b', markeredgewidth=2.)
loglog( num_solution[:,0], abs(real( num_solution[:,2])),'r^', label="FFT"    , markersize=14., markerfacecolor='None', markeredgecolor='r', markeredgewidth=2.)
loglog( num_solution[:,0], abs(imag( num_solution[:,2])),'g^',                  markersize=14., markerfacecolor='None', markeredgecolor='b', markeredgewidth=2.)

#xlim((min(theo_solution[:,1]), max(theo_solution[:,1])))
#xlim((1.e-5,1.e-1))
#ylim((1.e-4, 1.e1))
xlabel("$\\beta$")
ylabel("$\\omega$")

savefig("ShearAlfvenDamping.png")






