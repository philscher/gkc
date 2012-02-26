import scipy
import pylab 
import numpy

import numpy

def smoothFilter(x,window='hanning',  window_len=11):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    #y=numpy.convolve(w/w.sum(),s,mode='valid')
    y=numpy.convolve(w/w.sum(),s,mode='valid')
    return y






#http://www.swharden.com/blog/2008-11-17-linear-data-smoothing-in-python/
def smoothList(list,strippedXs=False,degree=10):  
     if strippedXs==True: return Xs[0:-(len(list)-(len(list)-degree+1))]  
     smoothed=[0]*(len(list)-degree+1)  
     for i in range(len(smoothed)):  
         smoothed[i]=sum(list[i:i+degree])/float(degree)  
     return smoothed  

def smoothListTriangle(list,strippedXs=False,degree=5):  
     weight=[]  
     window=degree*2-1  
     smoothed=[0.0]*(len(list)-window)  
     for x in range(1,2*degree):weight.append(degree-abs(degree-x))  
     w=numpy.array(weight)  
     for i in range(len(smoothed)):  
         smoothed[i]=sum(numpy.array(list[i:i+window])*w)/float(sum(w))  
     return smoothed  

def smoothListGaussian(list,strippedXs=False,degree=5):  
     window=degree*2-1  
     weight=numpy.array([1.0]*window)  
     weightGauss=[]  
     for i in range(window):  
         i=i-degree+1  
         frac=i/float(window)  
         gauss=1/(numpy.exp((4*(frac))**2))  
         weightGauss.append(gauss)  
     weight=numpy.array(weightGauss)*weight  
     smoothed=[0.0]*(len(list)-window)  
     for i in range(len(smoothed)):  
         smoothed[i]=sum(numpy.array(list[i:i+window])*weight)/sum(weight)  
     return smoothed  
"""
def smoothFilter(data,filterType,  degree):
    if   filterType == 'None'      : return array(data)
    elif filterType == 'Triangle'  : return smoothTriangle(data, degree=degree)
    elif filterType == 'List'      : return smoothList(data, degree=degree)
    elif filterType == 'Gaussian'  : return smoothListGaussian(data, degree=degree)
    else : print "Error no such smooth Filter"
"""
