import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

# generate some random data points, at regularly spaced time points
max_t = 11.0
min_t = 1.0
tspacing = 0.02
t = np.arange(min_t,max_t,tspacing) # remember this does not include max_t
npoints = len(t)

# linspace is probably more convenient than arange in many cases
# linspace is also guaranteed to include min_t and max_t, whereas
# arange is guaranteed to include min_t but will never contain max_t
#testingit = np.linspace(min_t,max_t,npoints)
#print t
#print testingit

sig = np.random.rand(npoints) #- 0.5

# generates consecutive integers starting from one and not exceeding 10
max_width = 11.0
min_width = 1.0
wspacing = 0.02
widths = np.arange(min_width,max_width,wspacing) # remember this does not include max_width
nwidths = len(widths)

# reference to a Ricker (i.e. Mexican Hat) wavelet function
wavelet = signal.ricker

# calculate CWT for 10 widths on 20 random points ( 20x10 matrix )
# This function implicitly assumes that the data is regularly spaced
cwtmatr = signal.cwt(sig, wavelet, widths)
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
              
print "max_index_j (position, x-axis): ", max_index_j, " position: ",t[max_index_j]         
print "max_index_i (scale, y-axis): ", max_index_i, " scale: ",widths[max_index_i]


# Plot the signal
# Also plot the wavelet with the overall best match with the signal
a = widths[max_index_i]
tt_half = t[max_index_j]/2.0
ttmax = 1.5*t[max_index_j]
ttmin = 0.5*t[max_index_j]
tstep = (ttmax-ttmin)/(20.0*a)
twavelet = np.arange(ttmin,ttmax,tstep)
Ft = signal.ricker(len(twavelet),a)

f1 = plt.figure()
ax1 = f1.add_subplot(111)
ax1.plot(t,sig,'ro-')
ax1.plot(twavelet,Ft,'bo--')
ax1.set_xlabel("Time")
ax1.set_ylabel("Intensity")
f1.show()

# Plot C(a,b)

# setup scaling on colormap
#
# It's important to get the extent correct because that will
# determine whether the scaling is correct (ie, aligned with the center of 
# the appropriate bin) on both axes.
min_x = t[0] - tspacing/2.0 # want to go to lower bound of first bin
max_x = max_t - tspacing + tspacing/2.0 # subtract tspacing to get final time point, add tspacing/2 to get upper bound of final bin
min_y = widths[0] - wspacing/2.0 # want to go to lower bound of first bin
max_y = max_width - wspacing + wspacing/2.0 # see two lines above
extent = [min_x,max_x,min_y,max_y]

f2 = plt.figure()
axim = f2.add_subplot(111)
# interpolation can be bilinear (default), nearest, bicubic
# bilinear and bicubic do some fancy gaussian weighting or something or another
# linear is simple - each box is only one color.
im = axim.imshow(cwtmatr,extent=extent,vmin=0,vmax=1.0,interpolation='nearest',origin='lower',aspect='auto')
f2.colorbar(im)
axim.set_xlabel("Position or Time")
axim.set_ylabel("Scale")
f2.show()