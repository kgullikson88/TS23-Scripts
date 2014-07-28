#!/usr/bin/python

from astropy.io import fits as pyfits
import pylab
import numpy as np
import scipy
import scipy.interpolate
import scipy.optimize
import sys
import os
import subprocess
import FitTellurics as utils
import MakeModel
import DataStructures
import FitsUtils
import Units
import RotBroad


#Get stellar model first, rotationally broadened
homedir = os.environ['HOME']
Bstarfile = homedir + "/School/Research/McDonaldData/BstarModels/BG19000g425v2.vis.7"
starmodel = RotBroad.Broaden(Bstarfile, 150*Units.cm/Units.km)
BSTAR_first = starmodel.x[0]
BSTAR_last = starmodel.x[-1]
BSTAR_fcn = scipy.interpolate.UnivariateSpline(starmodel.x, starmodel.y/starmodel.cont, s=0)
def BSTAR(x):
  value = BSTAR_fcn(x)
  value[x < BSTAR_first] = 1.0
  value[x > BSTAR_last] = 1.0
  return value

UsedLineList = "UsedLines.log"
outfile = open(UsedLineList, "w")
outfile.close()

    
class fitpoints:
  def __init__(self):
    self.x = []
    self.y = []    

#Main part of the code
class Improve:
  def __init__(self,filename):
    self.filename = filename
    hdulist = pyfits.open(filename)
    self.header = hdulist[0].header
    self.orders = FitsUtils.MakeXYpoints(self.header, hdulist[0].data)
    
    
  def ImportTelluric(self, filename):
    wave, trans = np.loadtxt(filename, usecols=(0,1), unpack=True)
    self.telluric = DataStructures.xypoint(wave.size)
    self.telluric.x = wave[::-1]*Units.nm/Units.um
    self.telluric.y = trans[::-1]
    self.telluric = FittingUtilities.ReduceResolution(self.telluric, 70000)
    
    
    
  #Here is the really important function.
  def Fit(self, plot=False):
    #main function. Will plot each order separately, and allow user to interact if plot=True
    self.clicks = []
    #interpolate the telluric (earth's atmospheric transmission) function
    Telluric = scipy.interpolate.UnivariateSpline(self.telluric.x, self.telluric.y, s=0)
    print "Plotting... press i to begin clicking points, and d when done"
    outfile = open("residuals.log", "w")
    outfile2 = open("UsedLines.log", "a")

    linelist = np.loadtxt(utils.LineListFile)
    
    #Loop over the spectral orders
    for i in range(37,51):
      print "Fitting order #" + str(i+1)
      self.orderNum = i
      self.fitpoints = fitpoints()
      wave = self.orders[i].x
      flux = self.orders[i].y/self.orders[i].cont
      tell = Telluric(wave)
      print "wave = ", wave
      
      #Do a cross-correlation first, to get the wavelength solution close
      ycorr = scipy.correlate(flux-1.0, tell-1.0, mode="full")
      xcorr = np.arange(ycorr.size)
      lags = xcorr - (flux.size-1)
      distancePerLag = (wave[-1] - wave[0])/float(wave.size)
      offsets = -lags*distancePerLag
      offsets = offsets[::-1]
      ycorr = ycorr[::-1]

      fit = np.poly1d(np.polyfit(offsets, ycorr, ycorr.size/100))
      ycorr = ycorr - fit(offsets)
      left = np.searchsorted(offsets, -1.0)
      right = np.searchsorted(offsets, +1.0)
      maxindex = ycorr[left:right].argmax() + left
      print "maximum offset: ", offsets[maxindex], " nm"
      pylab.plot(offsets, ycorr)
      pylab.show()
      userin = raw_input("Apply Cross-correlation correction? ")
      if "y" in userin:
        self.orders[i].x = self.orders[i].x + offsets[maxindex]
      
      #Fit using the (GridSearch) utility function
      data = DataStructures.xypoint(self.orders[i].x.size)
      data.x = np.copy(self.orders[i].x)
      data.y = np.copy(self.orders[i].y)
      data.cont = np.copy(self.orders[i].cont)
      fitfcn, offset = FitWavelength2(data, self.telluric, linelist)
      self.orders[i].x = fitfcn(self.orders[i].x - offset)
      
      #Let user fix, if plot is true
      if plot:
        #First, just plot all at once so user can examine fit
        left = np.searchsorted(self.telluric.x, self.orders[i].x[0])
        right = np.searchsorted(self.telluric.x, self.orders[i].x[-1])
        pylab.plot(self.orders[i].x, self.orders[i].y/self.orders[i].cont, label="data")
        pylab.plot(self.telluric.x[left:right], self.telluric.y[left:right]*BSTAR(self.telluric.x[left:right]), label="model")
        pylab.legend(loc='best')
        pylab.title("Order "+str(self.orderNum+1))
        pylab.show()
      
        #We only want to plot about 3 nm at a time
        spacing = 3.0
        data_left = 0
        data_right = np.searchsorted(self.orders[i].x, self.orders[i].x[0]+spacing)
        while (data_left < self.orders[i].x.size):
          #Bind mouseclick:
          self.fig = pylab.figure()
          self.clickid = self.fig.canvas.mpl_connect('button_press_event', self.onclick)

          left = np.searchsorted(self.telluric.x, self.orders[i].x[data_left])
          right = np.searchsorted(self.telluric.x, self.orders[i].x[min(self.orders[i].x.size-1, data_right)])
          pylab.plot(self.orders[i].x[data_left:data_right], self.orders[i].y[data_left:data_right]/self.orders[i].cont[data_left:data_right], label="data")
          pylab.plot(self.telluric.x[left:right], self.telluric.y[left:right]*BSTAR(self.telluric.x[left:right]), label="model")
          pylab.legend(loc='best')
          pylab.title("Order "+str(self.orderNum+1))
          pylab.show()
          data_left = data_right
          data_right = np.searchsorted(self.orders[i].x, self.orders[i].x[min(self.orders[i].x.size-1, data_left)]+spacing)


        #Once you close the window, you will get past the pylab.show() command
        #Fit the points to a cubic
        #This is done in a loop, to iteratively remove outliers
        done = False
        while not done:
          #self.fitpoints is filled when you are clicking in the window
          if (len(self.fitpoints.x) > 3):
            pars = np.polyfit(self.fitpoints.x, self.fitpoints.y, 3)
          else:
            pars = [0,1,0] #y=x... meaning don't try to improve on this order
          func = np.poly1d(pars)
          ignorelist = []
          x = np.array(self.fitpoints.x)
          y = np.array(self.fitpoints.y)
          resid = y - func(x) #residuals from the fit
          mean = resid.mean()
          std_dev = resid.std()
        
          #Find outliers (points with residuals over 0.01 or more than 2.5
          #   standard deviations from the mean
          for j in range(len(self.fitpoints.x)):
            residual = self.fitpoints.y[j] - func(self.fitpoints.x[j])
            if np.abs(residual) > 0.01 or np.abs(residual) > std_dev*2.5:
              ignorelist.append(j)
          if len(ignorelist) == 0:
            done = True
          else:
            for index in ignorelist[::-1]:
              print "removing point ", index, " of ", len(self.fitpoints.x)
              self.fitpoints.x.pop(index)
              self.fitpoints.y.pop(index)
            
        #Done removing outliers. Apply fit to the wavelengths
        print "y = ",pars[0],"x^2 + ",pars[1],"x + ",pars[2]
        self.orders[i].x = func(self.orders[i].x)
      
        #Output the residuals, and plot them. Make sure they look alright
        for j in range(len(self.fitpoints.x)):
          outfile.write(str(self.fitpoints.x[j])+"\t"+str(self.fitpoints.y[j])+"\t"+str(self.fitpoints.y[j]-func(self.fitpoints.x[j]))+"\n")
        outfile.write("\n\n\n\n")
        pylab.plot(self.fitpoints.x, self.fitpoints.y-func(self.fitpoints.x), 'ro')
        pylab.show()

        #Finally, add the lines to UsedLineList.log
        for line in self.fitpoints.y:
          outfile2.write("%.10g\n" %line)
          
        #Output after every order, in case program crashes
        FitsUtils.OutputFitsFile(self.filename, self.orders, func_order=5)
    outfile.close()
    
    #Output calibrated spectrum to file
    return FitsUtils.OutputFitsFile(self.filename, self.orders, func_order=5)
    
    
    
  #This function gets called when hit a key
  def keypress(self, event):
    if (event.key == "i"):
      #allow user to click on the canvas
      print "Mouse press active!"
      self.clickid = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
      return
    elif (event.key == "d"):
      #Done clicking. allow user to zoom to another set of lines
      print "Mouse press deactivated"
      self.fig.canvas.mpl_disconnect(self.clickid)
    elif (event.key == "r"):
      #User made a mistake in the previous click. Remove
      print 'Removing last click, which was at x = '
      if len(self.fitpoints.x) == len(self.fitpoints.y):
        print self.fitpoints.y.pop()
      else:
        print self.fitpoints.x.pop()
  
  #This function gets called when you click on the canvas, if the binding is active
  def onclick(self, event):
    tol = 0.04 #This is how close you have to be to the lowest point in the peak, in nm
    
    if len(self.fitpoints.x) == len(self.fitpoints.y):
      #Find lowest point in observed array:
      left = np.searchsorted(self.orders[self.orderNum].x, event.xdata - tol)
      right = np.searchsorted(self.orders[self.orderNum].x, event.xdata + tol)
      i = self.orders[self.orderNum].y[left:right].argmin()
      left = left + i
      left = np.searchsorted(self.orders[self.orderNum].x, self.orders[self.orderNum].x[left] - tol)
      right = np.searchsorted(self.orders[self.orderNum].x, self.orders[self.orderNum].x[left] + 2*tol)      
      centroid = self.Centroid(self.orders[self.orderNum].x[left:right], self.orders[self.orderNum].y[left:right])
      
      #Fit to gaussian:
      cont = 1.0
      depth = cont - self.orders[self.orderNum].y[(left+right)/2]/self.orders[self.orderNum].cont[i]
      mu = self.orders[self.orderNum].x[(left+right)/2]
      sig = 0.025
      params = [cont, depth, mu, sig]
      params,success = scipy.optimize.leastsq(ErrFunction, params, args=(self.orders[self.orderNum].x[left:right], self.orders[self.orderNum].y[left:right]/self.orders[self.orderNum].cont[left:right]))
      
      print "mean: ", params[2]
      self.fitpoints.x.append(params[2])
    else:
      #Find lowest point in observed array:
      left = np.searchsorted(self.telluric.x, event.xdata - tol)
      right = np.searchsorted(self.telluric.x, event.xdata + tol)
      i = self.telluric.y[left:right].argmin()
      left = left + i
      left = np.searchsorted(self.telluric.x, self.telluric.x[left] - tol)
      right = np.searchsorted(self.telluric.x, self.telluric.x[left] + 2*tol)
      centroid = self.Centroid(self.telluric.x[left:right], self.telluric.y[left:right])
      #Fit to gaussian:
      cont = 1.0
      depth = cont - self.telluric.y[(left+right)/2]
      mu = self.telluric.x[(left+right)/2]
      sig = 0.025
      params = [cont, depth, mu, sig]
      params,success = scipy.optimize.leastsq(ErrFunction, params, args=(self.telluric.x[left:right], self.telluric.y[left:right]))
      
      self.fitpoints.y.append(params[2])
      print "mean = ", params[2]
    return
  
  #Find centroid of the absorption line. I don't think I use this
  def Centroid(self, x, y):
    if x.size != y.size:
      print "Error! x and y not same size!"
      sys.exit()
    centroid = 0
    norm = sum(1/y)
    for i in range(x.size):
        centroid = centroid + x[i]/y[i]
    return centroid/norm
  

#Gaussian absorption line
def FitFunction(x,params):
  cont = params[0]
  depth = params[1]
  mu = params[2]
  sig = params[3]
  return cont - depth*np.exp(-(x-mu)**2/(2*sig**2))

#Returns the residuals between the fit from above and the actual values
def ErrFunction(params, x, y):
  return FitFunction(x,params) - y


#Second wavelength-fitting function that just shifts lines, instead of fitting them to gaussians
def WavelengthErrorFunction(shift, data, model):
  modelfcn = scipy.interpolate.UnivariateSpline(model.x, model.y, s=0)
  weight = 1.0/np.sqrt(data.y)
  weight[weight < 0.01] = 0.0
  newmodel = modelfcn(model.x + float(shift))
  if shift < 0:
    newmodel[model.x - float(shift) < model.x[0]] = 0
  else:
    newmodel[model.x - float(shift) > model.x[-1]] = 0
  returnvec = (data.y - newmodel)*weight
  return returnvec


#Gaussian absorption line
def GaussianFitFunction(x,params):
  cont = params[0]
  depth = params[1]
  mu = params[2]
  sig = params[3]
  return cont - depth*np.exp(-(x-mu)**2/(2*sig**2))

#Returns the residuals between the fit from above and the actual values
def GaussianErrorFunction(params, x, y):
  return GaussianFitFunction(x,params) - y



#Function to Fit the wavelength solution, using a bunch of telluric lines
#This assumes that we are already quite close to the correct solution
#Note: it comes from FitTellurics and is just slightly modified. 
#      Therefore, it takes FitTellurics structures
def FitWavelength2(order, telluric, linelist, tol=0.05, oversampling = 4, fit_order=3, max_change=2.0, debug=False):
  old = []
  new = []

  #Interpolate to finer spacing
  DATA_FCN = scipy.interpolate.UnivariateSpline(order.x, order.y, s=0)
  CONT_FCN = scipy.interpolate.UnivariateSpline(order.x, order.cont, s=0)
  MODEL_FCN = scipy.interpolate.UnivariateSpline(telluric.x, telluric.y, s=0)
  data = DataStructures.xypoint(order.x.size*oversampling)
  data.x = np.linspace(order.x[0], order.x[-1], order.x.size*oversampling)
  data.y = DATA_FCN(data.x)
  data.cont = CONT_FCN(data.x)
  model = DataStructures.xypoint(data.x.size)
  model.x = np.copy(data.x)
  model.y = MODEL_FCN(model.x)*BSTAR(model.x)
  
  #Begin loop over the lines
  for line in linelist:
    if line-tol > data.x[0] and line+tol < data.x[-1]:
      #Find line in the model
      left = np.searchsorted(model.x, line - tol)
      right = np.searchsorted(model.x, line + tol)
      minindex = model.y[left:right].argmin() + left

      mean = model.x[minindex]
      left2 = np.searchsorted(model.x, mean - tol*2)
      right2 = np.searchsorted(model.x, mean + tol*2)

      argmodel = DataStructures.xypoint(right2 - left2)
      argmodel.x = np.copy(model.x[left2:right2])
      argmodel.y = np.copy(model.y[left2:right2])

      #Do the same for the data
      left = np.searchsorted(data.x, line - tol)
      right = np.searchsorted(data.x, line + tol)
      minindex = data.y[left:right].argmin() + left

      mean = data.x[minindex]

      argdata = DataStructures.xypoint(right2 - left2)
      argdata.x = np.copy(data.x[left2:right2])
      argdata.y = np.copy(data.y[left2:right2]/data.cont[left2:right2])

      #Do a cross-correlation first, to get the wavelength solution close
      ycorr = scipy.correlate(argdata.y-1.0, argmodel.y-1.0, mode="full")
      xcorr = np.arange(ycorr.size)
      maxindex = ycorr.argmax()
      lags = xcorr - (argdata.x.size-1)
      distancePerLag = (argdata.x[-1] - argdata.x[0])/float(argdata.x.size)
      offsets = -lags*distancePerLag
      shift = offsets[maxindex]
      shift, success = scipy.optimize.leastsq(WavelengthErrorFunction, shift, args=(argdata, argmodel))
      if (debug):
        print argdata.x[0], argdata.x[-1], argdata.x.size
        print "wave: ", mean, "\tshift: ", shift, "\tsuccess = ", success
        pylab.plot(model.x[left:right]-shift, model.y[left:right], 'g-')
        pylab.plot(argmodel.x, argmodel.y, 'r-')
        pylab.plot(argdata.x, argdata.y, 'k-')
      if (success < 5):
        old.append(mean)
        new.append(mean + float(shift))
  if debug:
    pylab.show()
    pylab.plot(old, new, 'ro')
    pylab.show()
  #fit = UnivariateSpline(old, new, k=1, s=0)
  #Iteratively fit to a cubic with sigma-clipping
  fit = np.poly1d((1,0))
  mean = 0.0
  done = False
  while not done and len(old) > fit_order:
    done = True
    mean = np.mean(old)
    fit = np.poly1d(np.polyfit(old - mean, new, fit_order))
    residuals = fit(old - mean) - new
    std = np.std(residuals)
    #if debug:
    #  pylab.plot(old, residuals, 'ro')
    #  pylab.plot(old, std*np.ones(len(old)))
    #  pylab.show()
    badindices = np.where(np.logical_or(residuals > 2*std, residuals < -2*std))[0]
    for badindex in badindices[::-1]:
      print "Deleting index ", badindex+1, "of ", len(old)
      del old[badindex]
      del new[badindex]
      done = False

  #Check if the function changed things by too much
  difference = np.abs(order.x - fit(order.x - mean))
  if np.any(difference > max_change):
    fit = np.poly1d((1,0))
    mean = 0.0

  if debug:
    pylab.plot(old, fit(old - mean) - new, 'ro')
    pylab.show()
    pylab.plot(fit(order.x - mean), order.y, 'k-')
    pylab.plot(model.x, model.y, 'g-')
    print mean
    print fit
    pylab.show()
  
  return fit, mean



if __name__ == "__main__":
  #You will need to change these filename to to reflect where you store the models
  if sys.platform.startswith("linux"):
    outfilename = "/media/FreeAgent_Drive/TelluricLibrary/transmission-743.15-283.38-60.0-40.0-368.50-4.00-1.71-1.40"
  else:
    outfilename = "/Users/kgulliks//School/Research/lblrtm/run_examples/MyModel/OutputFiles/Generic.dat"

  if len(sys.argv)>1:
    for fname in sys.argv[1:]:
      improve = Improve(fname)
      improve.ImportTelluric(outfilename)
      improve.Fit(True)
  else:
    improve = Improve(raw_input("Enter file to calibrate: "))
    improve.ImportTelluric(outfilename)
    improve.Fit(True)
