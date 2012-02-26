import pylab
import mpmath as mp
import scipy
from PyPPL import *
# Plasma Dispersion function

def getGrowth(ky_list, Setup, disp="Gyro", init = -0.02 -.006j):
  mp.dps=15

  results  = []
  residuum = []
  eta     = Setup['eta']
  kx        = Setup['kx']
  v_te      = Setup['v_te']
  rho_te2   = Setup['rho_te2']
  tau       = Setup['tau']
  theta     = Setup['theta']
  m_ie      = Setup['m_ie']
  lambda_D2 = Setup['lambda_D2']

  for ky in ky_list:
    kp  = theta * ky * v_te 
        
    ko2 = kx**2 + ky**2
    kyp = ky/kp
    b   = ko2 * rho_te2
    G0  = Gamma0(b)
    G0mG1 = Gamma0(b) - Gamma1(b)
   

    if disp == 'Gyro1st':
      def DispersionRelation(w):
        
        Lambda = lambda_D2 * b + mp.exp(b) * (1. + 1. - G0) 
        zeta = w / kp
        Z = Z_PDF(zeta) 
        # Take care of extra pie due to normalization
        return  - (1. - eta/2. * (1. + b))*kyp * Z \
                - eta * kyp * (zeta + zeta**2 * Z) + zeta * Z + 1. +  Lambda 
    if disp == 'Gyro1stKin':
      def DispersionRelation(w):
        
        def Disp(w, kp, b, eta):
            zeta = w / kp
            Z  = Z_PDF(zeta)
            kyp = ky/kp
            return    - (1. - eta/2. * (1. + b))*kyp * Z \
                - eta * kyp * (zeta + zeta**2 * Z) + zeta * Z + 1. 
       
        # proton
        sum1MG0 = lambda_D2 * b + (1. - G0) +  (1. - G0/m_ie) 

        #return -mp.exp(-b) * Disp(w, kp, b, eta) +  mp.exp(-b/m_ie) * Disp(w,kp * mp.sqrt(m_ie), b/m_ie, eta*mp.sqrt(m_ie)) + sum1MG0
        return mp.exp(-b) * Disp(w, kp, b, eta) -  mp.exp(-b/m_ie) * Disp(w,kp * mp.sqrt(m_ie), b/m_ie, eta) + sum1MG0


    elif disp == 'Fluid':
      def DispersionRelation(w):
        K    = 1. + eta
        K_12 = 1. + eta/2.
        Lambda = 0.#lambda_D2 * b + mp.exp(b) * (1. + 1. - G0) 
        return -ko2 + ( (ky + w)/(K_12 * ky - w) + kp**2 /w**2 * (0.5 * (K * ky - w))  /(K_12 * ky - w)) + Lambda
        return -ko2 + ( (ky+w)/(K_12 * ky-w) + kp**2 /w**2 * (0.5 * (K * ky - w))  /(K_12 * ky-w))
    
     


    elif disp == 'Gyro':
      def DispersionRelation(w):
        Lambda = lambda_D2 *b  + (1. - G0) + 1.
        zeta = w / kp
        Z  = Z_PDF(zeta)
        # Take care of extra pie due to normalization
        return  - (1. - eta/2.)*kyp * G0 * Z +   eta * kyp * b * G0mG1 * Z   \
                - eta * kyp * G0 * ( zeta +  zeta**2 * Z) +  zeta * Z * G0  + G0  +  Lambda 
    elif disp == 'GyroKin':
      def DispersionRelation(w):
        
        # Take care of extra pie due to normalization
        # what is Debye lengtth effect here ?
        def Disp(w, kp, b, eta):
            zeta = w / kp
            Z  = Z_PDF(zeta)
            kyp = ky/kp
            G0    = Gamma0(b)
            G0mG1 = Gamma0(b) - Gamma1(b)
            return  - (1. - eta/2.)*kyp * G0 * Z +   eta * kyp * b * G0mG1 * Z   \
                    - eta * kyp * G0 * ( zeta +  zeta**2 * Z) +  zeta * Z * G0  + G0 
        
        sum1MG0 = (1. - G0) + (1. - Gamma0(b/m_ie))

        return Disp(w, kp, b, eta)  - Disp(w,kp * mp.sqrt(m_ie), b/m_ie, eta/mp.sqrt(m_ie)              ) + sum1MG0
        return Disp(w, kp, b, eta) + (1.-G0) + 1. # - Disp(w,kp * mp.sqrt(m_ie), b/m_ie, eta              ) + sum1MG0
        return Disp(w, kp, b, eta) - Disp(w,kp * mp.sqrt(m_ie), b/m_ie, eta/mp.sqrt(m_ie)) + sum1MG0
    else : print  "No such Dispersion Relation : ", disp
 
    omega = init

    def checkAndImprove(omega, solver, tol=1.e-9):
       omega2 = omega
       try  : omega2 = complex(mp.findroot(DispersionRelation, omega , solver=solver, maxsteps=100, tol=tol))
       except : pass
       if (abs(DispersionRelation(omega2)) < abs(DisperisonRelation(omega))) : return omega2
       else : return omega

    #omega = checkAndImprove(omega, "newton", 1.e-6)
    try:
                        if (ky < 0.3): omega = complex(mp.findroot(DispersionRelation, omega, solver='anewton', maxsteps=10, verify=False))
                        omega = complex(mp.findroot(DispersionRelation, omega, solver='halley', maxsteps=300))
                             
                             #omega = complex(mp.findroot(DispersionRelation, omega, solver='newton', maxsteps=100 ))
    except : 
      omega = float('nan')
    print "ky : ", ky , " w : ", omega 
    results.append(float(pylab.real(omega))  + 1.j * pylab.imag(omega))
    res = DispersionRelation(results[-1])
    residuum.append(min(1.,float(abs(res))))
  
  return (pylab.array(ky_list), pylab.array(results), pylab.array(residuum))
      
