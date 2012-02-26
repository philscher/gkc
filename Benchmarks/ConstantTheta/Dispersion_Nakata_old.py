import pylab
import mpmath as mp
import scipy




def Gamma0(b) : return mp.besseli(0,b)*mp.exp(-b)

# Plasma Dispersion function
def Z(x): return (1.j * mp.sqrt(mp.pi) * mp.exp(- x**2) * (1 + mp.erf(1.j * x)))



def getGrowthNakata(ky_list, Setup, init = -0.07 -.015j):
  
  results = []
  
  eta_e     = Setup['eta_e']
  kx        = Setup['kx']
  v_te      = Setup['v_te']
  rho_te2   = Setup['rho_te2']
  tau       = Setup['tau']
  theta     = Setup['theta']

  for ky in ky_list:
    
    kp  = mp.sqrt(2.) * theta * ky 
    # Dispersion Relation Equation (9)
    if ky > 2.5 : init = 0.01 - 0.01j
    def DispersionRelation(w):
        ko2 = kx**2 + ky**2
        def Lambda(b): return mp.exp(b) * (1. + tau - Gamma0(b))
        #def Lambda(b): return mp.exp(b) * (tau + b/(1.+b))
      
        #zeta = w / (kp * v_te)
        zeta = w / (kp * v_te)
        w_star_e = ky * pylab.sqrt(rho_te2) * v_te
 
        # Take care of extra pie due to normalization
        #return  1. + Lambda(ko2 * rho_te2) + zeta * Z(zeta) - ky/kp *  eta_e * zeta - ky/kp * ( eta_e * zeta**2 +\
        return  1. + Lambda(ko2 * rho_te2) + zeta * Z(zeta) - ky/kp *  eta_e * zeta - ky/kp * ( eta_e * zeta**2 +\
                    (1. - eta_e/2. * (1. + ko2 * rho_te2)))*Z(zeta) * mp.sqrt(mp.pi)
                    #(1. - eta_e/2. * (1. + ko2 * rho_te2)))*Z(zeta) 

    try:
        omega = complex(mp.findroot(DispersionRelation, init, solver='muller', maxsteps=1000))
        #omega = complex(PT.zermuller(DispersionRelation, 0., -0.2 - 0.05j, dx=.001)[0])

    except:
        omega = .0
        print "Not found : ", ky,  "  Theta : ", theta
    results.append(float(pylab.real(omega))  + 1.j * pylab.imag(omega))

  return (pylab.array(ky_list), pylab.array(results))
        
"""
for theta in [ 0.033, 0.050, 0.133, 0.15, 0.2 ]:
            Setup = { 'eta_e' : 6., 'kx' : 0.0, 'v_te' : 1., 'rho_te2' : 1., 'tau' : 1., 'theta' : 0. }
            ky_list = pylab.linspace(0.01, 2.0, 151)
            Setup['theta'] = theta*6.

            ky, gamma = getGrowthNakata(ky_list, Setup)
            pylab.plot(ky, pylab.imag(gamma), label='$\\theta=%.3f$' % theta)
            pylab.legend(ncol=1, loc='best').draw_frame(0)

            pylab.ylim((-0.01,0.045))
            pylab.xlim((0.,2.0))
            pylab.xlabel("Poloidal Wavenumber $k_y/\\rho_{te}$")
            pylab.ylabel("Growthrate $\\gamma_L [\\nu_{te}/L_T]$")
            pylab.savefig("Dispersion_Nakata.pdf", bbox_inches='tight')
            pylab.savefig("Dispersion_Nakata.png", bbox_inches='tight')
"""
