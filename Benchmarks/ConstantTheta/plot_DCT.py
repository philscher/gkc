import pylab
import mpmath as mp
import scipy
from PyPPL import *
from Dispersion_ConstTheta import *

# Plasma Dispersion function

#for theta in [ 0.033, 0.050, 0.133, 0.15, 0.2 ]:
for theta in [ 0.10]:
            Setup = { 'eta' : 4., 'kx' : 0.0, 'v_te' : mp.sqrt(2.), 'rho_te2' : 1., 'tau' : 1., 'theta' : 0. , 'm_ie' : 400., 'lambda_D2' : 0.}
            ky_list = pylab.logspace(pylab.log10(0.2), pylab.log10(40.), 101)
            #ky_list = pylab.linspace(20., 60.,101)
            Setup['theta'] = theta

            #ky, gamma, err = getGrowth(ky_list, Setup, disp="Gyro", init= 0.02 + 0.045j)
            #pylab.semilogx(ky, pylab.imag(gamma), 'b-', label='Gyro')
            #pylab.plot(ky, pylab.imag(gamma), 'b-', label='Gyro')

            #ky, gamma, err = getGrowth(ky_list, Setup, disp="Gyro1st", init=0.02+0.045j)
            #pylab.semilogx(ky[:350], pylab.imag(gamma)[:350], 'g-', label='Gyro $\\mathcal{O}{(1)}$')

            #ky, gamma, err = getGrowth(ky_list, Setup, disp="Fluid", init = 0.02 + 0.045j)
            #pylab.semilogx(ky, pylab.imag(gamma), 'r-', label='Fluid ')

            ky, gamma, err = getGrowth(ky_list, Setup, disp="GyroKin", init=-0.01 + 0.01j)
            #ky, gamma, err = getGrowth(ky_list, Setup, disp="GyroKin", init= 0.02 + 0.045j)
            pylab.semilogx(ky, pylab.imag(gamma), 'm-', label='Kinetic')
            #pylab.plot(ky, pylab.imag(gamma), 'm-', label='Kinetic')
            
pylab.legend(ncol=3, loc='best').draw_frame(0)

#pylab.ylim((-0.01,0.21))
#pylab.xlim((0.,10.))
pylab.xlabel("Poloidal Wavenumber $k_y/\\rho_{te}$")
pylab.ylabel("Growthrate $\\gamma_L [\\nu_{te}/L_n]$")
pylab.savefig("Dispersion_DCT.pdf", bbox_inches='tight')
pylab.savefig("Dispersion_DCT.png", bbox_inches='tight')

