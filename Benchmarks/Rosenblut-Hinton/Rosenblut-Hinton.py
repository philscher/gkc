
#reference GKW paper A.G. Peeters, Y. Camenen, F.J. Casson, W.A. Hornsby, A.P. Snodin, D. Strintzi and G. Szepesi
#Computer Physics Communications, 180, 2650 (2009)


from pylab import *



q = 1.3

eps = linspace(0.,0.5,201)


def Rosenbluth_Hinton(eps) :
  Theta = 1.6 * eps**(3./2.)
  return 1./(1.+q**2 * Theta/eps**2)

def Xiao_Catto(eps) :
  Theta = 1.6 * eps**(3./2.) + 0.5 * eps**2 + 0.36 * eps**(5./2.)
  return 1./(1.+q**2 * Theta/eps**2)



plot(eps, Rosenbluth_Hinton(eps), label="Rosenbluth-Hinton")
plot(eps, Xiao_Catto       (eps), label="Xiao-Catto")
xlabel("$\\epsilon")
ylabel("Bla")
legend()
savefig("Rosenbluth-Hinton.png")
