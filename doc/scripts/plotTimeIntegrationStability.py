"""
/*
 * =====================================================================================
 *       
 *        Project: GKC - Gyro-kinetic solver
 *
 *       Filename: TimeIntegrationStability.py
 *
 *    Description: Main file
 *
 *         Author: Paul P. Hilscher (2011)
 *
 *        License: GPLv3+
 * =====================================================================================
 */

"""




from pylab import *
from scipy.optimize import *
from mpmath import *
import pylab

rcParams['lines.linewidth']=4
font = {'family' : 'normal', 'weight' : 'bold', 'size' : 22 }
rc('font', **font)
fig = figure(figsize=(12, 6))

mp.dps = 60;

theta = linspace(1.e-3,2*pi,301)
fig = figure()


### Stability diagrams for Runge-Kutta Codes

def D(theta,n=4, zinit=0.) :
  p = taylor(exp, 0.0, n)
  def RKGrowth(a) :  return polyval(p[::-1], a) - exp(1.j*theta)
  return complex(findroot(RKGrowth, zinit, solver='muller'))

def plotRK(n):
    #z = [0.-1.j]
    z = [1.-1.j]
    for t in theta : z.append(D(n*t,n, z[-1]))
    print "RK : " , n , "  max : ", max(abs(array(z)))
    z = array(z)
   
    clf()
    pylab.fill(real(z[1:]), imag(z[1:]), color="#333333", label="RK-%i" % n)

    title("Stability region of explicit RK-%i integration" %n)	
    xlabel("Re$(\\lambda)$")
    ylabel("Im$(\\lambda)$")
   
    pylab.plot(linspace(0.,0., 128), linspace(-3.8, 3.8, 128), '-', color="#888888", linewidth=1.)
    xlim((-3, 1.))
    ylim((-3.8,3.8)) 
    pylab.savefig("TimeIntegration_RK%i_Stability.png" % n, bbox_inches='tight')

plotRK(2)
plotRK(3)
plotRK(4)
plotRK(5)
plotRK(6)


"""
## Euler backward stability region

clf()

theta = linspace(0.,2*pi,301)
fig = figure()



def Di(theta,n=4, zinit=0.) :
  return complex(1. - exp(- 1.j * theta))

def plotRK(n):
    #z = [0.-1.j]
    z = [0.-0.j]
    for t in theta : z.append(Di(n*t,n))
    print "RK : " , n , "  max : ", max(abs(array(z)))
    z = array(z)
    print shape(z)
    pylab.fill(real(z[1:]), imag(z[1:]), alpha=1., color="#FF" + str(n*11)+"00", label="RK-"+str(n))

plotRK(1)
xlabel("Re$(\\lambda)$")
ylabel("Im$(\\lambda)$")


#legend(loc='best', ncol=2).draw_frame(0)
#pylab.text(0.25,0.0, "Stability Region", ha='center', va='center', rotation='vertical')

#pylab.plot(linspace(0.,0., 128), linspace(-3.25, 3.25, 128), 'k-', linewidth=1.)
#ylim((-3.25,3.25))
pylab.savefig("ImplicitStability.pdf", bbox_inches='tight')

"""


## Euler backward stability region




