import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from math import sin
import random
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

#############################################
# function for generating marr wavelet
def marr(points,a):
    
    A = 1 / (np.sqrt(2 * np.pi) * (a**3)) # this is the only difference compared to the ricker wavelet
    wsq = a**2
    vec = np.arange(0, points) - (points - 1.0) / 2
    tsq = vec**2
    mod = (1 - tsq / wsq)
    gauss = np.exp(-tsq / (2 * wsq))
    total = A * mod * gauss
    return total
# end of function
#############################################

#############################################
# function for generating marr wavelet
def rickerTime(points,a,tspacing):
    
    
    A = 2 / (np.sqrt(3 * a) * (np.pi**0.25))
    wsq = a**2
    
    # The normal signal.ricker function assumes a certain
    # ranges of scales, so it just generates evenly spaced
    # points at increments of 1...this works fine if you're
    # using a > 1, but for smaller values of a this does
    # not sample the wavelet peak features nearly enough...
    #
    # I need to determine if some sort of relationship exists
    # between the value of a and the increment spacing
    
    #tmin = -5.0*a
    #tmax = 5.0*a
    #vec = np.linspace(tmin,tmax,points)  
    
    #vec = np.arange(0,points,tspacing) - tspacing*(points - 1.0) / 2
    vec = np.arange(0, points) - (points - 1.0) / 2
    vec *= tspacing
    
    #vec = np.arange(0, points) - (points - 1.0) / 2
    tsq = vec**2
    mod = (1 - tsq / wsq)
    gauss = np.exp(-tsq / (2 * wsq))
    total = A * mod * gauss
    return total
# end of function
#############################################

x = np.linspace(20,100,100)
npoints = len(x)
xspace = x[1] - x[0]
y = np.zeros([npoints])
for i,dx in enumerate(x):
    if x[i] > 51 and x[i] < 53:
        y[i] = sin(x[i]) + 1#+ (random.random()-0.5)*0.2
    elif x[i] > 35 and x[i] < 36:
        y[i] = 0.2*sin(2.0*x[i]) + 0.3
    else:
        y[i] = random.random()*0.2
        
# don't bother adjusting the x scaling, the signal.cwt function
# doesn't use any x data
#x -= x[0]
#x /= xspace
#print x
#y /= max(y)
#y *= 0.8
#y *= 2.0 / (2.5 * max(y))

# generates consecutive integers starting from one and not exceeding 10
#
# remember that the wavelet is generated using 10 * a points, so if a = 0.1
# then you're only constructing the wavelet from a single data point...
max_width = 8.0
min_width = 0.5 # should not need to even go lower than a=1.0 since applying normalization
wspacing = 0.1
widths = np.arange(min_width,max_width,wspacing) # remember this does not include max_width
nwidths = len(widths)

# reference to a Ricker (i.e. Mexican Hat) wavelet function
wavelet = rickerTime
#wavelet = marr

# calculate CWT for 10 widths on 20 random points ( 20x10 matrix )
# This function implicitly assumes that the data is regularly spaced
#cwtmatr = signal.cwt(y, wavelet, widths)
wavelet_points = 301
cwtmatr = np.zeros([len(widths), len(y)])
for ind, width in enumerate(widths):
    #wavelet_data = wavelet(min(10 * width, len(y)), width)
    wavelet_data = wavelet(wavelet_points, width, xspace)
    cwtmatr[ind, :] = signal.convolve(y, wavelet_data,
                                              mode='same')
#print cwtmatr.shape
#print len(cwtmatr)
#print cwtmatr
#print cwtmatr[9][0]

# Find max value in cwtmatr
max_val = 0.00
for i in range(nwidths):
    for j in range(npoints):
        if cwtmatr[i][j] > max_val:
            max_val = cwtmatr[i][j]
            max_index_i = i
            max_index_j = j
              
print "max_index_j (position, x-axis): ", max_index_j, " position: ",x[max_index_j]         
print "max_index_i (scale, y-axis): ", max_index_i, " scale: ",widths[max_index_i]
print "max_val= ",max_val


# Plot the signal
# Also plot the wavelet with the overall best match with the signal
a = widths[max_index_i]
#a = 1.0
#tmin = -5.0*a
#tmax = 5.0*a
#twavelet = np.linspace(tmin,tmax,wavelet_points)
twavelet = np.arange(0, wavelet_points) - (wavelet_points - 1.0) / 2
twavelet *= xspace
twavelet += x[max_index_j]
#a = 2
#tt_half = x[max_index_j]/2.0
#ttmax = 1.6*x[max_index_j]
#ttmin = 0.4*x[max_index_j]
#tstep = (ttmax-ttmin)/(200.0)
#twavelet = np.arange(ttmin,ttmax,tstep)
#twavelet = np.linspace(ttmin,ttmax,100)
Ft = rickerTime(wavelet_points,a,xspace)
print max(Ft)

f1 = plt.figure()
ax1 = f1.add_subplot(111)
ax1.plot(x,y)
ax1.plot(twavelet,Ft,'ro--')
ax1.set_xlabel("m/z")
ax1.set_ylabel("Intensity")
f1.show()

# Plot C(a,b)

# setup scaling on colormap
#
# It's important to get the extent correct because that will
# determine whether the scaling is correct (ie, aligned with the center of 
# the appropriate bin) on both axes.
spacing = x[1] - x[0]
min_x = x[0] - spacing/2.0 # want to go to lower bound of first bin
max_x = x[npoints-1] - spacing + spacing/2.0 # subtract tspacing to get final time point, add tspacing/2 to get upper bound of final bin
min_y = widths[0] - wspacing/2.0 # want to go to lower bound of first bin
max_y = max_width - wspacing + wspacing/2.0 # see two lines above
extent = [min_x,max_x,min_y,max_y]

f2 = plt.figure()
axim = f2.add_subplot(111)
# interpolation can be bilinear (default), nearest, bicubic
# bilinear and bicubic do some fancy gaussian weighting or something or another
# linear is simple - each box is only one color.
im = axim.imshow(cwtmatr,extent=extent,vmin=0,vmax=2.2,interpolation='nearest',origin='lower',aspect='auto')
f2.colorbar(im)
axim.set_xlabel("m/z")
axim.set_ylabel("Scale")
f2.show()

########################################################
# setup data arrays for 3d plot
X = np.zeros([len(widths),len(x)])
Y = np.zeros([len(widths),len(x)])
for ind1,val1 in enumerate(widths):
    for ind2,val2 in enumerate(x):
        X[ind1][ind2] = x[ind2]
        Y[ind1][ind2] = widths[ind1]
   
#X, Y = np.meshgrid(X, Y)
f3 = plt.figure()
ax = f3.add_subplot(111,projection='3d')
surf = ax.plot_surface(X, Y, cwtmatr, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
#ax.set_xlim(5500,5700)
ax.set_zlim(0,2.2)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
f3.colorbar(surf, shrink=0.5, aspect=5)
#ax.set_xlabel("m/z")
#ax.set_ylabel("Scale")
f3.show()
###########################################################