from pylab import *
import pylab

# open HDF5 File
rcParams['axes.labelsize'] = 18.
rcParams['axes.labelsize'] = 18.
rcParams['xtick.labelsize'] = 'large'
rcParams['ytick.labelsize'] = 'large'

color_indigo='#4B0082'

markers_D = [ '<', 's', 'v', '^', '>', '^', 'd', 'p', '*', 'h',  '1', '2', '3', '4', '+', 'x', 'D', '_', '|' ]
markers_C = [ 'r', 'g', 'b', 'c', 'm', 'y', 'k', '#CC0066', '#A21F12', '#9632F1' , '#1BC243', 'FF9900', '#92CC44',\
              '#99FFCC', '#AAFF33', '#FF00CC', '#FF33AA', '#33FFAA', '33AAFF' ]

colors_216 = ['000000',  '000033',  '000066',  '000099',  '0000CC',  '0000FF',
              '003300',  '003333',  '003366',  '003399',  '0033CC',  '0033FF',
              '006600',  '006633',  '006666',  '006699',  '0066CC',  '0066FF',
              '009900',  '009933',  '009966',  '009999',  '0099CC',  '0099FF',
              '00CC00',  '00CC33',  '00CC66',  '00CC99',  '00CCCC',  '00CCFF',
              '00FF00',  '00FF33',  '00FF66',  '00FF99',  '00FFCC',  '00FFFF',
              '330000',  '330033',  '330066',  '330099',  '3300CC',  '3300FF',
              '333300',  '333333',  '333366',  '333399',  '3333CC',  '3333FF',
              '336600',  '336633',  '336666',  '336699',  '3366CC',  '3366FF',
              '339900',  '339933',  '339966',  '339999',  '3399CC',  '3399FF',
              '33CC00',  '33CC33',  '33CC66',  '33CC99',  '33CCCC',  '33CCFF',
              '33FF00',  '33FF33',  '33FF66',  '33FF99',  '33FFCC',  '33FFFF',
              '660000',  '660033',  '660066',  '660099',  '6600CC',  '6600FF',
              '663300',  '663333',  '663366',  '663399',  '6633CC',  '6633FF',
              '666600',  '666633',  '666666',  '666699',  '6666CC',  '6666FF',
              '669900',  '669933',  '669966',  '669999',  '6699CC',  '6699FF',
              '66CC00',  '66CC33',  '66CC66',  '66CC99',  '66CCCC',  '66CCFF',
              '66FF00',  '66FF33',  '66FF66',  '66FF99',  '66FFCC',  '66FFFF',
              '990000',  '990033',  '990066',  '990099',  '9900CC',  '9900FF',
              '993300',  '993333',  '993366',  '993399',  '9933CC',  '9933FF',
              '996600',  '996633',  '996666',  '996699',  '9966CC',  '9966FF',
              '999900',  '999933',  '999966',  '999999',  '9999CC',  '9999FF',
              '99CC00',  '99CC33',  '99CC66',  '99CC99',  '99CCCC',  '99CCFF',
              '99FF00',  '99FF33',  '99FF66',  '99FF99',  '99FFCC',  '99FFFF',
              'CC0000',  'CC0033',  'CC0066',  'CC0099',  'CC00CC',  'CC00FF',
              'CC3300',  'CC3333',  'CC3366',  'CC3399',  'CC33CC',  'CC33FF',
              'CC6600',  'CC6633',  'CC6666',  'CC6699',  'CC66CC',  'CC66FF',
              'CC9900',  'CC9933',  'CC9966',  'CC9999',  'CC99CC',  'CC99FF',
              'CCCC00',  'CCCC33',  'CCCC66',  'CCCC99',  'CCCCCC',  'CCCCFF',
              'CCFF00',  'CCFF33',  'CCFF66',  'CCFF99',  'CCFFCC',  'CCFFFF',
              'FF0000',  'FF0033',  'FF0066',  'FF0099',  'FF00CC',  'FF00FF',
              'FF3300',  'FF3333',  'FF3366',  'FF3399',  'FF33CC',  'FF33FF',
              'FF6600',  'FF6633',  'FF6666',  'FF6699',  'FF66CC',  'FF66FF',
              'FF9900',  'FF9933',  'FF9966',  'FF9999',  'FF99CC',  'FF99FF',
              'FFCC00',  'FFCC33',  'FFCC66',  'FFCC99',  'FFCCCC',  'FFCCFF',
              'FFFF00',  'FFFF33',  'FFFF66',  'FFFF99',  'FFFFCC',  'FFAA0F']

################################## Set Environments for Plots #######################################



# good for presentation, bad for viewing
def setPlotOutputLarge():
    # open HDF5 File
    rcParams['xtick.major.size'] = 12.
    rcParams['xtick.minor.size'] =  5.
    rcParams['ytick.major.size'] = 12.
    rcParams['ytick.minor.size'] =  5.
    rcParams['axes.labelsize'] = 26.
    rcParams['axes.labelsize'] = 26.
    rcParams['xtick.labelsize'] = 'xx-large'
    rcParams['ytick.labelsize'] = 'xx-large'
    rcParams['lines.markersize']  = 22.   
    rcParams['font.size'] =           12.0
    rcParams['font.weight'] =           'bold'
    rcParams['lines.linewidth'] = 7.
    rcParams['legend.fontsize'] = 'xx-large'
    #rcParams['axes.color_cycle'] = markers_C

def setPlotOutputPublication():
    # open HDF5 File
    rcParams['xtick.major.size'] = 12.
    rcParams['xtick.minor.size'] =  5.
    rcParams['ytick.major.size'] = 12.
    rcParams['ytick.minor.size'] =  5.
    rcParams['axes.labelsize'] = 35.
    rcParams['xtick.labelsize'] = 'xx-large'
    rcParams['ytick.labelsize'] = 'xx-large'
    rcParams['lines.markersize']  = 20.   
    rcParams['lines.markeredgewidth']  = 5.   
    rcParams['font.size'] =         18.0
    rcParams['font.weight'] =     'medium'
    rcParams['lines.linewidth'] = 8.
    rcParams['legend.fontsize'] = 'x-large'
    rcParams['legend.fontsize'] = 'x-large'
    #rcParams['axes.color_cycle'] = markers_C


# good for presentation, bad for viewing
def setPlotOutputThesis():
    # open HDF5 File
    rcParams['xtick.major.size'] = 12.
    rcParams['xtick.minor.size'] =  5.
    rcParams['ytick.major.size'] = 12.
    rcParams['ytick.minor.size'] =  5.
    rcParams['axes.labelsize']   = 22.
    rcParams['axes.labelsize']   = 22.
    rcParams['xtick.labelsize']  = 'xx-large'
    rcParams['ytick.labelsize']  = 'xx-large'
    rcParams['lines.markersize'] = 12.   
    rcParams['lines.markeredgewidth']  = 0.7   
    rcParams['font.size'] =                  12.0
    rcParams['font.weight'] =           'bold'
    rcParams['lines.linewidth']   = 6.
    rcParams['legend.fontsize']   = '26'
    #rcParams['axes.color_cycle'] = markers_C


def newFigure(style="Thesis", ratio="1.41:1", basesize=7):

  rcdefaults()

  if   style == "Normal" : pass
  elif style == "Large"  : setPlotOutputLarge()
  elif style == "Thesis" : setPlotOutputThesis()
  elif style == "Publication" : setPlotOutputPublication()
  else  : 
    print "No such style"
    1./0.

  #if aspect == 'equal' : subplot(111, aspect='equal')
  #else : TypeError("Aspect only 'equal' allowed")


  if   ratio == ""       : figsize = (1.00*basesize,basesize)
  elif ratio == "2.33:1" : figsize = (2.33*basesize,basesize) 
  elif ratio == "2.11:1" : figsize = (2.11*basesize,basesize) 
  elif ratio == "1.85:1" : figsize = (1.85*basesize,basesize) 
  elif ratio == "1.50:1" : figsize = (1.50*basesize,basesize) 
  elif ratio == "1.41:1" : figsize = (1.41*basesize,basesize) 
  elif ratio == "1.33:1" : figsize = (1.33*basesize,basesize) 
  elif ratio == "5:3"    : figsize = (1.66*basesize,basesize) 
  elif ratio == "1.66:1"    : figsize = (1.66*basesize,basesize) 
  elif ratio == "1.00:1" : figsize = (1.00*basesize,basesize) 
  elif ratio == "1:1.33": figsize = (1.00*basesize, 1.33*basesize) 
  elif ratio == "1:1.55": figsize = (1.00*basesize, 1.55*basesize) 
  elif ratio == "1:1" : figsize = (1.00*basesize,basesize) 
  elif ratio == "1:3" : figsize = (basesize,3*basesize) 

  fig = figure(figsize=figsize)

  return fig




def plotContourWithColorbar(X,Y,data, **kwargs):
   """
        Plots time evolution of mode power.


        Optional keyword arguments:

        Keyword           Description
        ===============   ==============================================
         *dir*             Direction 'X' (for radial) or 'Y' for poloidal
         *modes*           List of modes (default plot modes). 
                           e.g. modes = [1,4,5]         - to plot all modes
                                modes = range(Nky)[::2] - to plot every second mode
         *field*           'phi' electric potential
                           'A' parallel magnetic vector potential
                           'B' parallel magnetic field
         *dir*             Direction 'X' (for radial) or 'Y' for poloidal
         *doCFL*           clear previous figure
         *label*           'ky' or 'm'
         *offset*          Offset due to zeroset to 2 .

   """
  
   norm          = kwargs.pop('norm', True)
   orientation   = kwargs.pop('orientation', 'horizontal')
   interpolation = kwargs.pop('interpolation', 'bicubic')
   cmap          = kwargs.pop('cmap' , cm.jet)
   clearPlot     = kwargs.pop('clearPlot' , True)
   CBexp         = kwargs.pop('CBExp' , True)
   lineLev       = kwargs.pop('lineLev' , 0)
 

   if CBexp == True:
        v_absmax =  max(abs(data.min()), abs(data.max()))
        exponent = int(log10(v_absmax)-1)
        data     = data / (10**exponent)

   if norm == True:
            v_minmax = abs(data).max()
            norm = mpl.colors.Normalize(vmin = -v_minmax, vmax = v_minmax)
   else : 
            norm = mpl.colors.Normalize(vmin = data.min(), vmax = data.max())
   if clearPlot == True  : clf()

   
   if (lineLev > 0) : 
     pylab.contour(X,Y,data, lineLev, colors='k', linewidths=1.)
     pylab.contourf(X,Y,data, lineLev, cmap=cmap, norm=norm, interpolation=interpolation, vmin=-v_minmax, vmax=v_minmax)
   else : 
     pylab.contourf(X,Y,data, 100, cmap=cmap, norm=norm, interpolation=interpolation, vmin=-v_minmax, vmax=v_minmax)
   # Set Colorbar so that in goes from -9.9 - +9.9 x exp
   #sf = ScalarFormatter(useOffset=True)
   #sf._set_data_interval(-10., 10.)
   #sf = FormatStrFormatter(useOffset=True, fmt="%1.1f")
   #v_norm = 10**int(log10(v_minmax))

   class MyStrFormatter(FormatStrFormatter):
         """
           Use a format string (fmt) to format the tick label, other options
           are the same as for ScalarFormatter
         """
         def format_data(self, x):
            #v_exp = int(log(v))
            return super(MyStrFormatter, self).format_data(1.) #x/10.**v_exp)
         def format_data_short(self, x):
            #v_exp = int(log(v))
            return super(MyStrFormatter, self).format_data_short(1.) #x/10.**v_exp)
   
   sf = MyStrFormatter(fmt="%1.1f")
   
   cb = colorbar(orientation=orientation, format=sf, norm=norm)
   
   if   CBexp and exponent == 1 : cb.ax.set_ylabel("$10\\times$")
   elif CBexp and exponent != 0 : cb.ax.set_ylabel("$10^{%i}\\times$" % exponent)

   xlim((X.min(), X.max()))
   ylim((Y.min(), Y.max()))

    

def plotZeroLine(x_start, x_end, direction='horizontal', color="#666666", lw=1.5):

        T = linspace(x_start, x_end, 301)
        if   direction == 'horizontal' : plot(T, zeros(len(T)), '-', linewidth=lw, color=color)
        elif direction == 'vertical'   : plot(zeros(len(T)), T, '-', linewidth=lw, color=color)
        else                           : TypeError("No such direction")

