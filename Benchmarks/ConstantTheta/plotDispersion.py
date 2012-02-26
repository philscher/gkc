import pylab
import mpmath as mp
import scipy
from Dispersion_ConstTheta import *

        

for m_ie in [ 0. ]:
            Setup = { 'eta' : 6., 'kx' : 0.0, 'v_te' : 1., 'rho_te2' : 1., 'tau' : 1., 'theta' : 0.01 , 'm_ie' : m_ie, 'lambda_D2' : 0.}
            ky_list = pylab.logspace(-1, 1.3, 101)

            disp="Gyro1stKin"
            if m_ie == 0. : disp ="Gyro1st"
            ky, gamma, res = getGrowth(ky_list, Setup, disp=disp)
            pylab.semilogx(ky, pylab.imag(gamma), label='$m_ie=%.3f$' % m_ie)
            
pylab.legend(ncol=1, loc='best').draw_frame(0)
pylab.ylim((-0.01,0.1))
#pylab.xlim((0.,2.0))
pylab.xlabel("Poloidal Wavenumber $k_y/\\rho_{te}$")
pylab.ylabel("Growthrate $\\gamma_L [\\nu_{te}/L_T]$")
pylab.savefig("Dispersion_Nakata.pdf", bbox_inches='tight')
pylab.savefig("Dispersion_Nakata.png", bbox_inches='tight')
