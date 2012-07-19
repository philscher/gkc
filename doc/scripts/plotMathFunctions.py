from pylab import *
import scipy.special

rcParams['lines.linewidth']=4
font = {'family' : 'normal', 'weight' : 'bold', 'size' : 22 }
rc('font', **font)

fig = figure(figsize=(12, 6))

def gplot(x_start, x_stop, f, ylabel_string, title_string, savename, ptype='linlin'):
	clf()
	if     ptype == 'linlin' : x = linspace(x_start, x_stop, 512)
	elif   ptype == 'linlog' : x = linspace(x_start, x_stop, 512)
	elif   ptype == 'loglin' : x = logspace(x_start, x_stop, 512)
	elif   ptype == 'loglog' : x = logspace(x_start, x_stop, 512)
	else   : print "No such ptype"
	xlabel("$x$")
	ylabel(ylabel_string)
	if   ptype == 'linlin' : 
				 plot(x, f(x))
				 xlim((x_start, x_stop))
	elif ptype == 'linlog' : 
				 semilogy(x,f(x))
				 xlim((x_start, x_stop))
	elif ptype == 'loglin' : 
				 semilogx(x,f(x))
				 xlim((10**x_start, 10**x_stop))
	elif ptype == 'loglog' : 
				 loglog(x,f(x))
				 xlim((10**x_start, 10**x_stop))
	else : print "No such plot type"
	
	title(title_string)
	savefig(savename, bbox_inches='tight')
	print "Plotted : " , title_string


################ functions 	

# Plot Bessel function J0
gplot(0., 42., scipy.special.j0, "$J_0(x)$", "Bessel Function", "Bessel_J0.png")
# Plot Bessel function J1
gplot(0., 42., scipy.special.j1, "$J_1(x)$", "Bessel Function", "Bessel_J1.png")
# Plot Bessel function I0
gplot(0, 8,  scipy.special.i0, "$I_0(x)$", "Bessel Function", "Bessel_I0.png", 'linlog')
# Plot Bessel function I1
gplot(0, 8,  scipy.special.i1, "$I_1(x)$", "Bessel Function", "Bessel_I1.png", 'linlog')
# Plot gyro-kinetic Gamma_0
func = lambda x : scipy.special.i0(x)*exp(-x)
gplot(-3, 2,  func, "$\\Gamma_0(x)$", "Bessel Function", "GK_Gamma0.png", 'loglin')
# Plot gyro-kinetic Gamma_1
func = lambda x : scipy.special.i1(x)*exp(-x)
gplot(-3, 2,  func, "$\\Gamma_1(x)$", "Bessel Function", "GK_Gamma1.png", 'loglin')
# Plot gyro-kinetic GK_Delta
func = lambda x : (scipy.special.i0(x)-scipy.special.i1(x))*exp(-x)
gplot(-3, 2,  func, "$\\Delta(x)$", "Bessel Function", "GK_Delta.png", 'loglin')
func = lambda x : (scipy.special.i0(x)-scipy.special.i1(x))*exp(-x)
gplot(-3, 2,  func, "$\\Delta(x)$", "Bessel Function", "GK_G0mG1.png", 'loglin')
func = lambda x : x/(x+1.)
gplot(-3, 2.,  func, "$1-\\Gamma_0(x) = \\frac{x}{1+x}$", "Bessel Function", "GK_1mG_Pade.png", 'loglin')



# Plot derivative of error function
x = linspace(0., 10., 512)
derf_dx = lambda x : 2./sqrt(pi) * exp(-x**2)

gplot(0., 10.,  derf_dx, "$erf'(x)$", "Derivative of Error function", "Deriv_ErrorFunction.png")

## Collisionality functions
Coll_F1 = lambda x : x*derf_dx(x) + (2. * x*x -1.) * scipy.special.erf(x)
gplot(-1., 1.,  Coll_F1, "$F_1$", "Collision F1", "Collision_F1.png", 'loglin')
	
Coll_F2 = lambda x : (1.-2./3.*x*x) * scipy.special.erf(x) - x* derf_dx(x)
gplot(-1., 1.,  Coll_F2, "$F_2$", "Collision F2", "Collision_F2.png", 'loglin')

Coll_F3 = lambda x : Coll_F1(x) + 3. * Coll_F2(x)
gplot(-1, 1.,  Coll_F3, "$F_3$", "Collision F3", "Collision_F3.png", 'loglin')

# Plot Chandrasekhar function
clf()
kp = linspace(0., 10., 512)

Chandra = lambda x : (scipy.special.erf(x) - x * derf_dx(x))/(2.*x**2)
plot(x, Chandra(x))
xlabel("$x$")
ylabel("$Chandra(x)$")
savefig("Chandrasekhar.png", bbox_inches='tight')





#

 
