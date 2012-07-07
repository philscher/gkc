from pylab import *
import scipy.special


def gplot(x_start, x_stop, f, ylabel_string, title_string, savename, ptype='linlin'):
	clf()
	if   ptype == 'linlin' : x = linspace(x_start, x_stop, 512)
	else                   : x = logspace(x_start, x_stop, 512) 
	xlabel("$x$")
	ylabel(ylabel_string)
	xlim((x_start, x_stop))
	if   ptype == 'linlin' : plot(x, f(x))
	elif ptype == 'linlog' : semilogy(x,f(x))
	elif ptype == 'loglin' : semilogx(x,f(x))
	elif ptype == 'loglog' : loglog(x,f(x))
	else : print "No such plot type"
	title(title_string)
	savefig(savename, bbox_inches='tight')


################ functions 	

# Plot Bessel function J0
gplot(0., 42., scipy.special.j0, "$J_0(x)$", "Bessel Function", "Bessel_J0.png")
# Plot Bessel function J1
gplot(0., 42., scipy.special.j0, "$J_1(x)$", "Bessel Function", "Bessel_J1.png")
# Plot Bessel function I0
gplot(0, 100,  scipy.special.i0, "$I_0(x)$", "Bessel Function", "Bessel_I0.png", 'linlog')
# Plot Bessel function I1
gplot(0, 100,  scipy.special.i1, "$I_1(x)$", "Bessel Function", "Bessel_I1.png", 'linlog')
# Plot gyro-kinetic Gamma_0
func = lambda x : scipy.special.i0(x)*exp(-x)
gplot(-3, 2,  func, "$\Gamma_0(x)$", "Bessel Function", "GK_Gamma0.png", 'loglin')
# Plot gyro-kinetic Gamma_1
func = lambda x : scipy.special.i1(x)*exp(-x)
gplot(-3, 2,  func, "$\Gamma_1(x)$", "Bessel Function", "GK_Gamma1.png", 'loglin')
# Plot gyro-kinetic GK_Delta
func = lambda x : (scipy.special.i0(x)-scipt.special.i1)*exp(-x)
gplot(-3, 2,  func, "$\Delta(x)$", "Bessel Function", "GK_Delta.png", 'loglin')
func = lambda x : (scipy.special.i0(x)-scipt.special.i1)*exp(-x)
gplot(-3, 2,  func, "$\Delta(x)$", "Bessel Function", "GK_G0mG1.png", 'loglin')


# Plot Gamma Function
clf()
x = linspace(0., 10., 512)

semilogx(x, x/(1.+x))
xlabel("$x$")
ylabel("$J_1(x)$")
savefig("GK_1mG_Pade.png", bbox_inches='tight')





# Plot Modified Bessel function I1
clf()
x = linspace(0., 10., 512)

plot(x, scipy.special.i1(x))
xlabel("$x$")
ylabel("$J_1(x)$")
savefig("Bessel_I1.png", bbox_inches='tight')




# Plot derivative of error function
x = linspace(0., 10., 512)
derf_dx = lambda x : 2./sqrt(pi) * exp(-x**2)

clf()
plot(x, derf_dx(x))
xlabel("$x$")
ylabel("$erf'(x)$")
savefig("Deriv_ErrorFunction.png", bbox_inches='tight')


# Plot Chandrasekhar function
clf()
kp = linspace(0., 10., 512)

Chandra = lambda x : (scipy.special.erf(x) - x * derf_dx(x))/(2.*x**2)
plot(x, Chandra(x))
xlabel("$x$")
ylabel("$Chandra(x)$")
savefig("Chandrasekhar.png", bbox_inches='tight')


#

 
