from pylab import *

matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"] 
import matplotlib

x = linspace(0., 1., 64)

# 3rd order Hermite functions



plot(x, 2*x**3-3*x**2+1, '-', label='$H^3_{00}=\;2x^3-3x^2+1$' )
plot(x,-2*x**3+3*x**2  , '-', label='$H^3_{10}=- 2x^3+3x^2$'   )
plot(x,   x**3-2*x**2+x, '-', label='$H^3_{01}=\;\;x^3-2x^2+x$')
plot(x,   x**3-x**2    , '-', label='$H^3_{11}=\;\;x^3-x^2$'   )

xlabel("$x$")
ylabel("$H^3_{nm}(x)$")
legend(loc='best', ncol=1).draw_frame(0)

savefig("HermiteInterpolation_3rdOrderFunctions.png", bbox_inches='tight')

clf()

plot(x,-6*x**5 + 15*x**4 -10*x**3 + 1.   , label="$H^5_{00} = -6x^5 + 15 x^4 - 10 x^3 + 1$")
plot(x, 6*x**5 - 15*x**4 +10*x**3        , label="$H^5_{10} =  6x^5 - 15 x^4 + 10 x^3$"    )
plot(x, -3.*x**5 + 8.* x**4 - 6* x**3 + x, label="$H^5_{01} = -3x^5 + 8 x^4 - 6 x^3 + x  $")
plot(x, -3.*x**5 + 7 * x**4 - 4* x**3    , label="$H^5_{11} =  -3x^5 + 7 x^4 - 4 x^3$"    )

# tfrac is not working, why ? included in preemble
plot(x, - 0.5 * x**5 + 1.5 * x**4 - 1.5* x**3 + 0.5 * x**2, \
    label="$H^5_{02} =  -\\frac{1}{2}x^5 + \\frac{3}{2} x^4 - \\frac{3}{2} x^3 + \\frac{1}{2} x^2$")
plot(x,   0.5 * x**5 - x**4 + 0.5 * x**3, label="$ H^5_{12} =  -\\frac{1}{2}x^5 - x^4 + \\frac{1}{2} x^3$")

xlabel("$x$")
ylabel("$H^5_{nm}(x)$")
legend(loc='best', ncol=1).draw_frame(0)


savefig("HermiteInterpolation_5thOrderFunctions.png", bbox_inches='tight')
