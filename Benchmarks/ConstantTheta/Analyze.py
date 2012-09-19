##########################################################################  
#
#   Compare simulation of local, constant-theta geometry
#   with theoretical dispersion.
#
#   Author  : Paul P. Hilscher (2012)
#   License : GNU Public License v3 (or at your opinion any later version)
#
##########################################################################

# get revision number !!

def getDispersionGyro1(theta, eta, kx, ky, doGyro1=False):

  from PyPPL import PlasmaDispersion, getZero, Gamma0,Gamma1
  from numpy import exp
  kp  = np.sqrt(2.) * theta * ky 
  b   = kx**2 + ky**2
  kyp = ky/kp
  lambda_D2 = 0. 
  G0    = Gamma0(b)
  G0mG1 = Gamma0(b) - Gamma1(b)

  def DR_Gyro(w):
        
      zeta = w / kp
      Z = PlasmaDispersion(zeta) 
     
      if (doGyro1):
        Lambda = lambda_D2 * b + exp(b) * (1. + 1. - G0) 
        R =   - (1. - eta/2. * (1. + b))*kyp * Z \
                - eta * kyp * (zeta + zeta**2 * Z) + zeta * Z + 1. +  Lambda 
      else :
        Lambda = lambda_D2 *b  + (1. - G0) + 1.
        R =   - (1. - eta/2.)*kyp * G0 * Z +   eta * kyp * b * G0mG1 * Z   \
                - eta * kyp * G0 * ( zeta +  zeta**2 * Z) +  zeta * Z * G0  + G0  +  Lambda 
      return R
        
        # Take care of extra pie due to normalization
  if   ky > 0.33 : init = [0.05+0.05j, 0.0+0.0j, 0.1+0.1j]
  elif ky < 0.15 : init = [0.01+0.01j, -0.01+0.01j, 0.0+0.0j]
  else           : init = [0.000+0.04j, 0.002+0.045j, 0.001+0.035j]
  
  return getZero(DR_Gyro, init=init)
  

import sys
import tables
import gkcData
import gkcLinear
import pylab
import numpy as np


if (len(sys.argv[:]) == 1):
  print "Use as python Analyze.py File_1.h5 File_2.h5 ..."

results = []
# Open Files and extract frequency and growthrates
for fileh_name in sys.argv[1:]:

            fileh5  = tables.openFile(fileh_name)

            theta = fileh5.root.Geometry._v_attrs.Theta[0]
            eta   = fileh5.root.Species.cols.w_T[1]

            D = gkcData.getDomain(fileh5)
            
            omega_sim =  gkcLinear.getFrequencyGrowthrates(fileh5, start=250, stop=-1)

            useGyro1 = (fileh5.root.Grid._v_attrs.Nm[0] == 1)

            # get Theory
            ky_list = np.logspace(np.log10(D['ky'][1]), np.log10(D['ky'][-1]), 128)

            omega_th = []
            for ky in ky_list :  omega_th.append(getDispersionGyro1(theta, eta, 0., ky, useGyro1))
            omega_th = np.array(omega_th)
           
            results.append( (theta, omega_sim, omega_th) )
            

            #if(fileh5.root.Grid._v_attrs.Nm[0] > 1) : _disp="Gyro"
            #kye, y, e = Dispersion_ConstTheta.getGrowth(ky_list, Setup, disp=_disp)
            #result_Theory.append((kye, y, e, Setup['theta']))
  
            fileh5.close()


pylab.figure(figsize=(10,8))

# Plot Growthrates
pylab.subplot(211, title="Benchmark in local, constant $\\Theta$-Geometry using gkc++ (rev. 150)")
for (theta, omega_sim, omega_th) in results:
      
            pylab.semilogx(D['ky'][1:-1], np.imag(omega_sim[1:-1]), '^', label='$\\Theta=%.1f$' % theta)
            pylab.semilogx(ky_list, np.imag(omega_th), '-', color='#444444', alpha=0.8, linewidth=5.)

pylab.xlim((D['ky'][1], D['ky'][-1]))
pylab.ylim((-0.1, 0.25))
pylab.legend(loc='best').draw_frame(0)
pylab.xlabel("$k_y$")
pylab.ylabel("Growthrate  $\\omega_i [L_n / v_{te}]$")


# Plot Frequency
pylab.subplot(212)
for (theta, omega_sim, omega_th) in results:

            pylab.semilogx(D['ky'][1:-1], np.real(omega_sim[1:-1]), '^-')
            pylab.semilogx(ky_list, np.real(omega_th), '-', color='#444444', alpha=0.8, linewidth=5.)
pylab.xlim((D['ky'][1], D['ky'][-1]))
pylab.ylim((0., 2.))
pylab.xlabel("$k_y$")
pylab.ylabel("Frequency  $\\omega_r [L_n / v_{te}]$")

pylab.savefig("ConstantTheta_Gyro.png", bbox_inches='tight')

