from scipy.optimize import *
from scipy.special import *
import numpy as np
import mpmath as mp

 
def getDispersion(ky, parameters):  
  
  from PyPPL import PlasmaDispersion, getZero, Gamma0,Gamma1
  from numpy import exp

  mode     = parameters['mode']
  m_ie     = parameters['m_ie']
  tau      = parameters['tau']
  kx       = parameters['kx']
  kp       = parameters['kp']
    
  
  beta     = parameters['beta']
  lambda_D = parameters['lambdaD']
  rhoLn    = parameters['rhoLn']
  rhoLT    = parameters['rhoLT']

  def DR_Smoly(w):

    k_ortho2 = kx**2 + ky**2
    mass_e   = 1./m_ie

    # ion  diamagnetic frequency
    w_ni    = - ky/kp * rhoLn * 0.5
    w_Ti    = - ky/kp * rhoLT * 0.5 
    s_i     = w 
    b_i     = (kx**2 + ky**2)/ 2.
    
    w_ni    = - ky/kp * rhoLn 
    w_Ti    = - ky/kp * rhoLT 
    s_i     = w / sqrt(2.) 
    b_i     = (kx**2 + ky**2)
        
    #electron
    w_ne   = ky/kp * rhoLn 
    w_Te   = ky/kp * rhoLT
    s_e    =  w / sqrt(2./mass_e)  
    b_e    = (kx**2 + ky**2) * mass_e
        
    # Plasma skin depth
    delta =  tau/(beta) / m_ie #sqrt(m_ie)
    
    ######################## Dispersion relation for gyrokinetic EM mode
    
    # settings for the ions
    
    Gamma_0i = Gamma0(b_i) 
    Gamma_1i = Gamma1(b_i)
    Z_i      = PlasmaDispersion(s_i)
   
    
    l_i = 1. - (1. - w_ni / w) * Gamma_0i - w_Ti/w * (Gamma_0i - Gamma_1i ) * b_i
    
    D_i = (1. - w_ni/w) * ( 1. + s_i * Z_i) * Gamma_0i  + w_Ti/w * s_i * \
           (0.5 * Z_i - s_i - s_i**2 * Z_i ) * Gamma_0i \
           + w_Ti/w *  ( 1. + s_i * Z_i) * ( Gamma_0i - Gamma_1i) * b_i

    # settings for the electrons
    Gamma_0e = Gamma0(b_e) 
    Gamma_1e = Gamma1(b_e) 
    Z_e      = PlasmaDispersion(s_e) 
    
    l_e = 1. - (1. - w_ne / w) * Gamma_0e - w_Te/w * (Gamma_0e - Gamma_1e ) * b_e
    
    D_e = (1. - w_ne/w) * ( 1. + s_e * Z_e) * Gamma_0e  + w_Te/w * s_e * \
         (0.5 * Z_e - s_e - s_e**2 * Z_e ) * Gamma_0e \
         + w_Te/w * ( 1. + s_e * Z_e) * ( Gamma_0e - Gamma_1e) * b_e


    #############################################################
    # Standard Dispersion relation
    if(mode == "EM"):
            d =    k_ortho2 * delta**2 * (l_i * tau + l_e + D_i * tau + D_e ) - 2 * s_e**2 * (D_i * tau + D_e) * (l_i * tau + l_e ) \
                 + k_ortho2 * lambda_D**2 * ( k_ortho2 * delta**2 - 2 * s_e**2 * (D_i * tau + D_e)) 
    
        # approximation for electrons
    elif(mode == "ETG"):
            d =  1. + l_e + D_e + lambda_D**2 * k_ortho2 
   
    # bad charge is wrong sign
    elif(mode == "ITG"):
            d =  1. + l_i + D_i + lambda_D**2 * k_ortho2 
    # simplified ITG dispersion relation
    # approximation for electrons
    elif(mode == "ITGe"):
            d =  l_e + D_e + tau * (l_i + D_i) 
    else: 
         raise NameError('No such Dispersion Relation')
   
    return d
  
  # Take care of extra pie due to normalization
  #init, solver =  (0.01 +.0026j, 0.02 + 0.03j, 0.015 + 0.01j), 'muller'
  #init, solver =  (0.01 +.0026j, 0.02 + 0.03j, 0.015 + 0.01j), 'muller'
  
  init, solver =  0.01 + .03j, 'muller'
  
  return getZero(DR_Smoly, init=init, solver=solver)





import sys
import tables
import gkcData
import gkcLinear
import pylab
import numpy as np


if (len(sys.argv[:]) == 1):
  print "Use as python Analyze.py File_1.h5 File_2.h5 ..."

results = []



##################### Call for adiabatic case ######################

# Open Files and extract frequency and growthrates

fileh5  = tables.openFile(sys.argv[1])


isKinetic = (len(fileh5.root.Species.cols[:]) == 3)

rhoLn   = fileh5.root.Species.cols.w_n[1]
rhoLT   = fileh5.root.Species.cols.w_T[1]

#kp      = fileh5.root.Species.cols.w_T[1]
D = gkcData.getDomain(fileh5)

omega_sim =  gkcLinear.getFrequencyGrowthrates(fileh5, start=50, stop=-1)

# get Theory
ky_list = np.logspace(np.log10(D['ky'][1]), np.log10(D['ky'][-1]), 48)

if (isKinetic) : 
     mode = 'EM'
     m_ie = fileh5.root.Species.cols.Mass[1]/fileh5.root.Species.cols.Mass[2]
else           :
     mode = 'ETG'
     m_ie = 1.

param = {  'disp' : 'GyroSlab', 'mode' : mode, 'beta'  : 1.0e-5,  'tau' : 1., 'lambdaD' : 1.e-5, 'rhoLn' : rhoLn, \
          'rhoLT' : rhoLT, 'kx' : np.sqrt(2.) * 1.e-1 , 'kp' : 2. * np.sqrt(2.) * 1.e-3, 'm_ie' : m_ie, 'adiab' : lambda x : 1 }


omega_th = []

for ky in ky_list :  
  omega_th.append(getDispersion(ky, param))
  #omega_th.append(getDispersion(sqrt(1837.) * ky, param))
omega_th = np.array(omega_th) * param['kp']#/sqrt(1837.)

fileh5.close()


pylab.figure(figsize=(10,8))

# Plot Growthrates
pylab.subplot(211, title="Benchmark of Smolykov et al. (PRL, 2002)  using gkc++ (rev. 173)")
     


pylab.semilogx(D['ky'][1:-1], np.imag(omega_sim[1:-1]), '^')
pylab.semilogx(ky_list, np.imag(omega_th), '-', color='#444444', alpha=0.8, linewidth=5.)

pylab.xlim((D['ky'][1], D['ky'][-1]))
pylab.ylim((-0.01, 0.02))
#pylab.legend(loc='best').draw_frame(0)
pylab.xlabel("$k_y$")
pylab.ylabel("Growthrate  $\\omega_i [L_n / v_{te}]$")


# Plot Frequency
pylab.subplot(212)


pylab.semilogx(D['ky'][1:-1], np.real(omega_sim[1:-1]), '^-')
pylab.semilogx(ky_list, np.real(omega_th), '-', color='#444444', alpha=0.8, linewidth=5.)

pylab.xlim((D['ky'][1], D['ky'][-1]))
#pylab.ylim((0., 2.))
pylab.xlabel("$k_y$")
pylab.ylabel("Frequency  $\\omega_r [L_n / v_{te}]$")

if(isKinetic) : pylab.savefig("Smolyakov_Kinetic.png"  , bbox_inches='tight')
else          : pylab.savefig("Smolyakov_Adiabatic.png", bbox_inches='tight')




