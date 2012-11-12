from pylab import *
import sys
import os

figure(figsize=(20,10))

input_filename = sys.argv[1]

# read-in file from command line argument
fl = open(input_filename).readlines()

r = linspace(0., 1., 257)

def getArrayFromLines(start, stop):
   eq_lines =  fl[start:stop]

   # join lines list together
   eq_lines = " ".join(eq_lines)

   # replace 1.1e-02-1.2e-03 -> 1.1e-02 -1.2-e03
   eq_lines = eq_lines.replace('E','e')
   eq_lines = eq_lines.replace('e-','em').replace('-',' -').replace('em','e-')

   # read in array
   eq = [float(a) for a in eq_lines.split()]

   return array(eq)



a1 = getArrayFromLines(2, 58)
a2 = getArrayFromLines(58, 109)
a3 = getArrayFromLines(110, 161)
a4 = getArrayFromLines(162, 214)
q  = getArrayFromLines(13424, 13475)
a6 = getArrayFromLines(13476, 13599)

subplot(241)
print shape(a1)
plot(a1)
xlabel("r/a")
ylabel("Te[eV]")
subplot(242)
print shape(a2)
plot(a2)
subplot(243)
print shape(a3)
plot(a3)
subplot(244)
print shape(a4)
plot(a4)
subplot(245)
print shape(1)
plot(q)
xlabel("r/a")
ylabel("Safety factor q")
subplot(246)
r = linspace(0., 1., len(q))
plot(r/q*gradient(q,1./257.))
xlabel("r/a")
ylabel("Shear")

subplot(247)
plot(a6)
subplot(248)

############### Plot equilibrium equilibrium ##########

eq = getArrayFromLines(213, 13423)
eq = array(eq).reshape((257,257))


contour(eq, 50)


# Split e
filename = input_filename.split(os.extsep, 1)[0]

savefig(filename + ".pdf", bbox_inches='tight')


