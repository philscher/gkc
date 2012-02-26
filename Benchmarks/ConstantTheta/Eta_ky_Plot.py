from scipy.optimize import *
from scipy.special import *
from scipy.integrate import *
from numpy import *
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt




from Dispersion_Nakata import *

def PDF(x) : return (1.j * sqrt(pi) * exp(- x**2) * (1. - erf(-1.j * x)))


#for theta in [ 0.033, 0.050, 0.133, 0.15, 0.2 ]:
ky_list = pylab.linspace(0.1, 3., 80)
eta_list = linspace(2.,6., 80)
theta = 0.1
Z  = zeros((len(eta_list), len(ky_list)), dtype=complex)

for i in range(len(eta_list)):
            Setup = { 'eta_e' : eta_list[i], 'kx' : 0.0, 'v_te' : mp.sqrt(2.), 'rho_te2' : 1., 'tau' : 1., 'theta' : theta , 'm_ie' : 1837.}

            #ky, gamma, err = getGrowthNakata(ky_list, Setup, disp="Gyro")
            #pylab.semilogx(ky, pylab.imag(gamma), label='Adiabatic $\\theta=%.3f$' % theta)
            ky, gamma, err = getGrowthNakata(ky_list, Setup, disp="GyroKin")
            pylab.semilogx(ky, pylab.imag(gamma),  label='Kinetic   $\\theta=%.3f$' % theta)
            Z[i,:] = nan_to_num(gamma)
   #
#norm = normalize(vmin = -v_minmax, vmax = v_minmax)
#contourf(log10(ky_list), mie_list, solution_EM, 400, rect=rect, cmap=grow, norm=norm)
#contourf(log10(ky_list), mie_list, solution_EM, 400, rect=rect, cmap=grow, norm=norm)
#contourf(log10(ky_list), mie_list, solution_EM, 400, rect=rect, cmap=grow, norm=norm)
  
#xlim((log10(ky_list)[0], log10(ky_list)[-1]))

#

fig = plt.figure()
ax = Axes3D(fig)
X,Y = np.meshgrid(eta_list,ky_list)
print real(Z.T)
pl = ax.plot_surface(X, Y, real(Z.T), cstride=1, rstride=1, linewidth=1., antialiased=True, cmap=cm.jet)
ax.set_zlabel("$Z(\\zeta)$")


ax.set_xlabel("$\\eta$")
ax.set_ylabel("$k_y$")

fig.savefig('PlasmaDispersionFunction3D_real.pdf')
fig.savefig('PlasmaDispersionFunction3D_real.png')

fig = plt.figure()

ax = Axes3D(fig)
X,Y = np.meshgrid(eta_list,ky_list)

print imag(Z.T)
norm = normalize(vmin = -0.02, vmax=imag(Z.T).max())
pl = ax.plot_surface(X, Y, imag(Z.T), cstride=1, rstride=1, linewidth=1., antialiased=True, cmap=cm.jet, norm=norm)
ax.set_zlim3d((-0.02, imag(Z.T).max()*1.1))
ax.set_zlabel("$Z(\\zeta)$")


ax.set_xlabel("$\\eta$")
ax.set_ylabel("$k_y$")

fig.savefig('PlasmaDispersionFunction3D_imag.pdf')
fig.savefig('PlasmaDispersionFunction3D_imag.png')

