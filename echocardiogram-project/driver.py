import numpy as np #matplotlib/pylab imports numpy automatically
import fit #auxillary file that contains the user-defined functions
import sys #enables system functions to be called (eg, sys.exit())
import matplotlib.pyplot as plt

# import data
# python uses backslash for escape character
torig, Corig = np.loadtxt('C:\\Users\\frenchwr\\Documents\\Python-Scripts\\ECG Project\\roi1-nostdev.dat').T
#Ctmp = np.zeros(len(C))
#Can't just use Ctmp = C, this seems to treat
#Ctmp as an alias for C
#Ctmp[:] = C[:]

# first determine an appropriate baseline
# for now let's just take the minimum raw
# intensity value. This avoids negative intensities.
Cmin = Corig.min()
Cshifted = Corig - Cmin

# Filter the data
# How about a low-pass filter
# This filters out the high-frequency (> 2 Hz) noise
# See: http://www.scipy.org/Cookbook/FIRFilter
(t,C) = fit.low_pass_filter(torig,Corig)

plt.plot(torig,Cshifted,linewidth=1.0)
plt.plot(t,C,linewidth=3.0)
plt.legend(['unfiltered','filtered'])
plt.show()
#sys.exit()

#C = np.zeros(len(Ctmp)) #remember this is zero-indexed
# w-point sliding average
#w = 1
#fit.sAverage( w , C )

# Now define the fit region
(istart,iend) = fit.getFitRange( C , 0.2, 0.1 ) #call function getFitRange
Crange = C[istart:iend]

# Remember: the time scale is arbitrary in the original data.
#           we want to plot everything relative to the injection
#           time. so we optimize the fit while shifting the time
#           scale, and then pick the time axis where the mean-square-
#           error is minimized. think of this less as the injection
#           time but the point at which contrast agent first enters
#           the region of interest.
tshift = np.linspace(0,t[istart],300) 
mse = np.zeros(len(tshift)) #mean square error array
ierr = 1
for i in range(len(tshift)):
  t_tmp = t - tshift[i]
  trange = t_tmp[istart:iend]
  (mse[i],ierr) = fit.runFit(trange,Crange)
  if ierr==-1:
    print "fitting stopped"
    break
  else:  
    print "i= ",i," Mean Square Error= ",mse[i]," tshift[i]= ",tshift[i]

if i==0:
  sys.exit("No fits obtained...stopping fitting program")

# use i as the upper bound, accounts for the situation where
# fitting fails halfway through
min_mse = mse[0:i].min()
print "\nFitting Complete\n"
print "min_mse= ",min_mse
print "index at min= ",mse[0:i].argmin()
print "injection time= ",tshift[mse[0:i].argmin()]

plotInjectionOpt = 1
if (plotInjectionOpt):
  f3 = plt.figure()
  ax3 = f3.add_subplot(111)
  ax3.plot(tshift[0:i-1],mse[0:i-1])
  ax3.set_xlabel("t - t$_{inj.}$ (sec)")
  ax3.set_ylabel("Mean Square Error (Echo Intensity)")
  f3.show()

plotResults = 1 # either 1 or 0
# plot the results
if (plotResults):

  # use optimal fit parameters
  t_tmp = t - tshift[mse[0:i].argmin()]
  trange = t_tmp[istart:iend]
  # calculate the X matrix and Y vector
  (X,Y) = fit.getXY(trange,Crange)
  # fit the data
  print " \nOptimal Parameters:\n "
  (Yfit,Cfit,sq_error,ierr) = fit.getFit(X,Y,Crange)
  # now plot everything
  f1 = plt.figure()
  f2 = plt.figure()
  ax1 = f1.add_subplot(111)
  ax1.plot(trange,Y,trange,Yfit)
  ax1.legend(['Data','Fit'])
  ax1.set_xlabel("Time (sec)")
  ax1.set_ylabel("ln(C(X1)) + 0.5*ln(X1)")
  ax2 = f2.add_subplot(111)
  ax2.plot(t_tmp,C,'black')
  ax2.plot(trange,Cfit,'red',linewidth=4.0)
  ax2.legend(['Data','Fit'])
  ax2.set_xlabel("t - t$_{inj.}$ (sec)")
  ax2.set_ylabel("Echo Intensity")
  f1.show()
  f2.show()
