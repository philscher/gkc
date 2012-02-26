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
    #kp  = 1.e-1
        
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
        def Disp(w, kp, b, eta, m_ie):
            eta = eta #/ mp.sqrt(m_ie)
            b = b / m_ie
            kp = kp * mp.sqrt(m_ie)
            zeta = w / kp
            Z  = Z_PDF(zeta)
            kyp = ky*mp.sqrt(m_ie)/kp
            G0    = Gamma0(b)
            G0mG1 = Gamma0(b) - Gamma1(b)
            Lambda =  (1. - G0) 
            return  - (1. - eta/2.)*kyp * G0 * Z +   eta * kyp * b * G0mG1 * Z   \
                    - eta * kyp * G0 * ( zeta +  zeta**2 * Z) +  zeta * Z * G0  + G0 + Lambda
       
        Proton   = Disp(w, kp, b, eta, 1.)  - 1./(1.-G0)
        Electron = Disp(w, kp, b, eta, m_ie) 
        #Electron = mp.sqrt(m_ie)*Disp(w, kp, b, eta, m_ie) + (1. - Gamma0(b/m_ie))

        return Proton + Electron + b * lambda_D2
    elif disp == 'GyroKinI':
      def DispersionRelation(w):
        
        # Take care of extra pie due to normalization
        # what is Debye lengtth effect here ?
        def Disp(w, kp, b, eta, m_ie, q):
            eta = eta #/ mp.sqrt(m_ie)
            b = b / m_ie
            kp = kp * mp.sqrt(m_ie)
            zeta = w / kp
            Z  = Z_PDF(zeta)
            kyp = ky*mp.sqrt(m_ie)/kp
            G0    = Gamma0(b)
            G0mG1 = Gamma0(b) - Gamma1(b)
            return  q * (- (1. - eta/2.)*kyp * G0 * Z +   eta * kyp * b * G0mG1 * Z   \
                    - eta * kyp * G0 * ( zeta +  zeta**2 * Z) +  zeta * Z * G0  + G0) 
       
        Proton   = Disp(w, kp, b, eta, 1., 1.)  +    (1. - G0) 
        Electron = -Disp(w, kp, b, eta, m_ie, -1.) + (1. - G0/m_ie) 
        #Electron = mp.sqrt(m_ie)*Disp(w, kp, b, eta, m_ie) + (1. - Gamma0(b/m_ie))
        return Proton + Electron #+ b * lambda_D2
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


    def Disp(w, kp, b, eta, m_ie, q):
            eta = eta #/ mp.sqrt(m_ie)
            b = b / m_ie
            kp = kp * mp.sqrt(m_ie)
            zeta = w / kp
            Z  = Z_PDF(zeta)
            kyp = ky*mp.sqrt(m_ie)/kp
            G0    = Gamma0(b)
            G0mG1 = Gamma0(b) - Gamma1(b)
            Lambda =  (1. - G0) 
            return  q * (- (1. - eta/2.)*kyp * G0 * Z +   eta * kyp * b * G0mG1 * Z   \
                    - eta * kyp * G0 * ( zeta +  zeta**2 * Z) +  zeta * Z * G0  + G0) + Lambda
       
    results.append(float(pylab.real(omega))  + 1.j * pylab.imag(omega))
    res = DispersionRelation(results[-1])
    residuum.append(min(1.,float(abs(res))))

    ###
    """
    def Disp(w, kp, b, eta, mass):
            eta = eta / mp.sqrt(m_ie)
            b = b / m_ie
            kp = mp.sqrt(m_ie)
            zeta = w / kp
            Z  = Z_PDF(zeta)
            kyp = ky*mp.sqrt(m_ie)/kp
            G0    = Gamma0(b)
            G0mG1 = Gamma0(b) - Gamma1(b)
            return  - (1. - eta/2.)*kyp * G0 * Z +   eta * kyp * b * G0mG1 * Z   \
                    - eta * kyp * G0 * ( zeta +  zeta**2 * Z) +  zeta * Z * G0  + G0 
        
        
       
    print " proton   ky : " , ky , "  " , complex(Disp(omega,kp,b,eta, 1.)  )
    print " electron ky : " , ky , "  " , complex(mp.sqrt(m_ie)*Disp(omega,kp,b,eta,m_ie)              )  
    """
    #print " electron ky : " , ky , "  " , complex(Disp(omega,kp * mp.sqrt(m_ie), b/m_ie, eta/mp.sqrt(m_ie)              )  + 1. - Gamma0(b/m_ie))
  
  return (pylab.array(ky_list), pylab.array(results), pylab.array(residuum))
     


#for theta in [ 0.033, 0.050, 0.133, 0.15, 0.2 ]:
for theta in [ 0.10]:
            Setup = { 'eta' : 4., 'kx' : 0.0, 'v_te' : mp.sqrt(2.), 'rho_te2' : 1., 'tau' : 1., 'theta' : 0. , 'm_ie' : 100., 'lambda_D2' : 0.}
            ky_list = pylab.logspace(pylab.log10(0.2), pylab.log10(10.), 51)
            #ky_list = pylab.linspace(20., 60.,101)
            Setup['theta'] = theta

            ky, gamma, err = getGrowth(ky_list, Setup, disp="Gyro", init= 0.02 + 0.045j)
            pylab.semilogx(ky, pylab.imag(gamma), 'b-', label='Gyro')
            #pylab.plot(ky, pylab.imag(gamma), 'b-', label='Gyro')

            #ky, gamma, err = getGrowth(ky_list, Setup, disp="Gyro1st", init=0.02+0.045j)
            #pylab.semilogx(ky[:350], pylab.imag(gamma)[:350], 'g-', label='Gyro $\\mathcal{O}{(1)}$')

            #ky, gamma, err = getGrowth(ky_list, Setup, disp="Fluid", init = 0.02 + 0.045j)
            #pylab.semilogx(ky, pylab.imag(gamma), 'r-', label='Fluid ')

            #ky, gamma, err = getGrowth(ky_list, Setup, disp="GyroKin", init=-0.01 + 0.01j)
            #ky, gamma, err = getGrowth(ky_list, Setup, disp="GyroKin", init= 0.02 + 0.045j)
            #pylab.semilogx(ky, pylab.imag(gamma), 'm-', label='Kinetic')
            #pylab.semilogx(ky, pylab.real(gamma), 'r-', label='Kinetic')
            #pylab.twinx()
            ky, gamma, err = getGrowth(ky_list, Setup, disp="GyroKinI", init=-0.01 + 0.01j)
            #ky, gamma, err = getGrowth(ky_list, Setup, disp="GyroKin", init= 0.02 + 0.045j)
            pylab.semilogx(ky, pylab.imag(gamma), 'b-', label='Kinetic-I')
            #pylab.semilogx(ky, pylab.real(gamma), 'r-', label='Kinetic')
            #pylab.ylim((-0.01,1.))
            #pylab.plot(ky, pylab.imag(gamma), 'm-', label='Kinetic')
            
pylab.legend(ncol=3, loc='best').draw_frame(0)

#pylab.xlim((0.,10.))
pylab.xlabel("Poloidal Wavenumber $k_y/\\rho_{te}$")
pylab.ylabel("Growthrate $\\gamma_L [\\nu_{te}/L_n]$")
pylab.savefig("Dispersion_DCT.pdf", bbox_inches='tight')
pylab.savefig("Dispersion_DCT.png", bbox_inches='tight')

