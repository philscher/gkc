
def plotStreamLines(ntime=-1, nz=4, fileh5=fileh[0]):

    Ap = fileh5.root.Potential.Ap[nz,:,:,ntime]
    dA_dy, dA_dx = gradient(Ap)
  
    D = getDomain(fileh5)
    rect = [0.08, 0.1, 0.75, 0.8 ] 
    x = linspace(0.0, D['Lx'], D['Nx'])
    y = linspace(0.0, D['Ly'], D['Ny'])
    conax = axes(rect, axisbg='w')

    streamplot.streamplot(x,y, dA_dy, -dA_dx)

    xlim((x[0], y[-1]))
    ylim((x[0], y[-1]))
