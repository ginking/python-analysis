import numpy as np # interestingly I need this here
from scipy.optimize import leastsq
from math import exp, pi, sqrt # also surprising that i need this
from scipy.signal import kaiserord, lfilter, firwin, freqz
#import matplotlib.pyplot as plt
# note, pyton is zero-indexed

########################################################
# function for calculating npoints-point sliding average
def sAverage( npoints, Ct ) :
  for i in range(len(Ct)) : 
      if i>npoints and i<(len(Ct)-npoints-1) :
        for j in range(npoints) :
            Ct[i] += Ct[i-j] + Ct[i+j]
        Ct[i] = Ct[i]/(2*npoints + 1)      
      #print "Ct[i]= ",Ct[i]
# end of sliding average function
########################################################

########################################################
# function for defining the fit range
def getFitRange( Ct , weight_lo, weight_hi ) :
  Ct_max = Ct.max()
  print "Cshifted_max= ",Ct_max

  Cstart = weight_lo * Ct_max # 0.2 works reasonably well for dataset1/roi1
  Cend = weight_hi * Ct_max # 0.2 works reasonably well for dataset1/roi1

  istart = 0
  iend = len(Ct) - 1
  trigger = 0
  # now find 0.1Cshifted_max (on curve rise) 
  for i in range(len(Ct)): 
      if Ct[i] > Cstart and istart == 0:
        istart = i
        print "istart= ",istart
        #print "Cshifted[i]= ",Cshifted[i]
        print "Cshifted[i]/Cshifted_max= ",Ct[i]/Ct_max
      if Ct[i] == Ct_max:
        trigger = 1
      if Ct[i] < Cend and trigger == 1:
        iend = i
        trigger = 0
        print "iend= ",iend
        #print "Cshifted[i]= ",Cshifted[i]
        print "Cshifted[i]/Cshifted_max= ",Ct[i]/Ct_max

  return (istart,iend)     
# end of function 
########################################################

########################################################
# This function defines the residual between data and
# a curve that's fit to it of the form:
# y = a - b*x1 - c*x2 
# where a, b, and c are the optimized parameters
#
def residuals(p, y, x1, x2):
  a, b, c = p
  err = y - (a - b*x1 - c*x2)
  return err
# end of function residuals
########################################################

########################################################
# Sets up the X and Y arrays for fitting
########################################################
def getXY(trange,Crange):
  X = np.ones((len(Crange),3)) #remember this is zero-indexed
  Y = np.zeros(len(Crange)) #remember this is zero-indexed
  for i in range(len(Crange)): #this loop starts at i=0, works like for (i=0;i<n;i++)
      X[i,1] = trange[i]
      X[i,2] = 1/trange[i]
      Y[i] = np.log(Crange[i]) + 0.5*np.log(X[i,1]) #log is ln by default
  return (X,Y)
# end function
########################################################

########################################################
# Function for calculating the model prediction
########################################################
def getC(ParameterGamma,ParameterLambda,ParameterMoverQ,t):
	result = (ParameterMoverQ/ParameterGamma) * exp(ParameterLambda) 
	result *= sqrt( ParameterLambda*ParameterGamma / (2*pi*t))  
	result *= exp( -1.0*(ParameterLambda/2) * (t/ParameterGamma + ParameterGamma/t) )
	return result
# end function
########################################################

########################################################
# Performs fitting and extracts parameter values
# Returns the fitted curves for plotting, further 
# analysis, etc
########################################################
def getFit(X,Y,Crange):
  # for discussion/examples of optimization techniques see:
  # http://docs.scipy.org/doc/scipy/reference/tutorial/optimize.html
  p0 = [5, 1, 1]
  # call scipy least squares fit function
  P = leastsq( residuals, p0, args=(Y,X[:,1],X[:,2]) )
  print "P= ",P

  Yfit = np.zeros(len(Crange)) #remember this is zero-indexed
  Cfit = np.zeros(len(Crange)) #remember this is zero-indexed

  if P[0][0] < 0 or P[0][1] < 0 or P[0][2] < 0:
    return (Yfit,Cfit,100000,-1)
  else:
    # real parameters
    ParameterGamma = sqrt(P[0][2]/P[0][1])
    ParameterLambda = 2*ParameterGamma*P[0][1]
    ParameterMoverQ = sqrt(pi/P[0][1]) * exp(P[0][0]-2*sqrt(P[0][1]*P[0][2]))
    print "ParameterGamma= ",ParameterGamma
    print "ParameterLambda= ",ParameterLambda
    print "ParameterMoverQ= ",ParameterMoverQ

    # now write out the fitted results
    sq_error = 0
    for i in range(len(Crange)):
		# Yfit = ln(C) + 0.5ln(X1)
		Yfit[i] = P[0][0] - P[0][1]*X[i,1] - P[0][2]*X[i,2]
		Cfit[i] = getC(ParameterGamma,ParameterLambda,ParameterMoverQ,X[i,1])
		sq_error += (Cfit[i] - Crange[i]) * (Cfit[i] - Crange[i])

    return (Yfit,Cfit,sq_error,1)
# end function
########################################################

########################################################
# Controls the fitting procedure by calling a few
# helper functions
########################################################
def runFit(trange,Crange):

  # calculate the X matrix and Y vector
  (X,Y) = getXY(trange,Crange)

  # fit the data
  (Yfit,Cfit,sq_error,ierr) = getFit(X,Y,Crange)

  # normalize for mean square error of fit
  sq_error = sq_error / len(Crange)
  return (sq_error,ierr) # return mean square error
# end function
########################################################

############################################################
# function to filter high frequency noise from a signal
############################################################
def low_pass_filter(torig,Corig):
  sample_rate = 50 # per second
  nyq_rate = sample_rate / 2.0 # The Nyquist rate of the signal
  width = 5.0/nyq_rate # the desired width of the transition from pass to stop,
                     # relative to the Nyquist rate
  ripple_db = 60.0 # the desired attenuation in the stop band, in dB
  N, beta = kaiserord(ripple_db, width) # compute the order and kaiser parameter
                                      # for the FIR filter
  cutoff_hz = 2 # the cutoff frequency of the filter
  taps = firwin( N, cutoff_hz/nyq_rate, window=('kaiser', beta)) # use firwin with a Kaiser
                                                               # window to create a
                                                               # lowpass FIR filter
  C = lfilter(taps, 1.0, Corig)
  delay = 0.5 * (N-1) / sample_rate
  t = torig - delay
  return (t,C)
#############################################################
# end function
#############################################################
