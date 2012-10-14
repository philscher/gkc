import numpy  as np
import mpmath as mp
import scipy.special
# ************************ Plasma Dispersion Functions *****************

"""
# Now Plasma Dispersion Function
# Calculated according to Newberger, Comp. Phys. Comm. (1986)
"""
def PlasmaDispersion(z, n=0):
  #if(hasattr(x,'__iter__')):
  # Values from Callen , Fundamentals of Plasma Physics

  Z =  1.j * mp.sqrt(mp.pi) * mp.exp(- z**2) * mp.erfc(-1.j * z)

  if   n == 0 : return  Z
  elif n == 1 : return  1. + z * Z 
  elif n == 2 : return  z + z**2 * Z
  elif n == 3 : return  0.5 * ( 1. + 2. * z**2 * ( 1. + z * Z))
  else        : raise TypeError("BUG : Plasma Dispersion function only for N <= 3 defined") 
  """  
  if  (x.real >= 0. and x.imag >= 0.)    : Z = W(x)
  elif(x.real > 0.  and x.imag < 0.)      : Z = conj(W(-z)) 
  elif(x.real < 0.  and abs(x.imag) > 0.) : Z = 2. * exp(-z**2) - W(-z) 
  else                                   : Z = 0.
  print x.imag, " : " , x.imag < 0.
  return Z
  """
PDF = np.vectorize(PlasmaDispersion, otypes=[np.complex])



"""
    Two Pole Approximation of the Plasma Dispersion Function
"""
def PlasmaDispersion_TwoPole(X):
    A1 = -1.2359 - 1.2150j    
    B1 = +0.5468 - 0.0372j
    A2 = -0.3786 - 1.3509j
    B2 = -1.0468 + 2.1018j
    A3 =  0.3786 - 1.3509j
    B3 = -1.0468 - 2.1018j
    A4 =  1.2359 - 1.2150j
    B4 =  0.5468 + 0.0372j
    
    Z=B1/(X-A1)+B2/(X-A2)+B3/(X-A3)+B4/(X-A4)
    return Z


# approximation for small values of x
def PlasmaDispersion_Small(z):
    return 1.j * mp.sqrt(mp.pi) * mp.exp(-z**2) - 2. * z * ( 1. - 2./3.*z**2 + 4./15.*z**4 - 8./105.*z**6)		

# approximation for large values of x
def PlasmaDispersion_Large(z):
    #sigma = []
    #if(hassattr(x,'__iter__')):
    #    for xv in x:
    #        r = real(xv)
    #        i = imag(xv)
    #        sigmav = 0.
    #        sigma.append(sigmav)
    #sigma = mp.where(imag(z) > 1./abs(real(z)), 0, where(abs(imag(z)) < 1./abs(real(z)), 1., 2.))
    
    x = z.real
    y = z.imag
    sigma = 0.;            
    if   (y      >  1./abs(x)) : sigma = 0.
    elif (abs(y) <  1./abs(x)) : sigma = 1.
    elif (y      < -1./abs(x)) : sigma = 2.
    else                       : print "Something wrong"

    try :
        #return 1j * mp.sqrt(mp.pi) * sigma * mp.exp(-z**2) - 1./z * ( 1. + 1./(2.*z**2) + 3./(4.*z**4) + 15./(8.*z**6) + 105./(16.*z**8))		
        return 1j * mp.sqrt(mp.pi) * sigma * mp.exp(-z**2) - 1./z * ( 1. + 1./(2.*z**2) + 3./(4.*z**4) + 15./(8.*z**6))		
    except OverflowError:
        return 0.e0



# *********************** Gamma Functions **************************

#def Gamma0(b) : return    scipy.special.i0e(float(b))  #mp.besseli(0,b)*mp.exp(-b)
#def Gamma1(b) : return    scipy.special.i1e(float(b))  #mp.besseli(1,b)*mp.exp(-b)
def Gamma0(b) : return    mp.besseli(0,b)*mp.exp(-b)
def Gamma1(b) : return    mp.besseli(1,b)*mp.exp(-b)
#def Gamma1(b) : return    #mp.besseli(1,b)*mp.exp(-b)


def getZero(Func, init=1.+1.j, solver='muller'):
  import mpmath as mp
  w0 = complex(float('nan'), float('nan'))
  try:
     w0 = mp.findroot(lambda w : Func(complex(w)), (init, init+0.1, init+0.1j), solver=solver, maxsteps=600, ftol=1.e-3)
  except:
   try:
     w0 = mp.findroot(lambda w : Func(complex(w)), init, solver='halley', maxsteps=600, ftol=1.e-3)
   except:   
     print "Could not find zero"

  return complex(w0)
# ***************** Other
"""
# Number density, charge, mass
def PlasmaFreq(n, mass = const.m_e, charge=const.e ):
    return sqrt(charge**2 * n / (mass * const.epsilon_0 ))

# Gyration radius
def GyroFreq(B, mass = const.m_e, charge=const.e):
  return mass * B / charge
"""
