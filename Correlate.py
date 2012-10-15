import subprocess
import sys
import os
import numpy
from collections import defaultdict
from scipy.interpolate import UnivariateSpline
from scipy.signal import fftconvolve
import DataStructures
import Units
#import MakeModel

class Resid:
  def __init__(self, size):
    wave = numpy.zeros(size)
    rect = numpy.zeros(size)
    opt = numpy.zeros(size)
    recterr = numpy.zeros(size)
    opterr = numpy.zeros(size)
    cont = numpy.zeros(size)

#This function rebins (x,y) data onto the grid given by the array xgrid
def RebinData(data,xgrid):
  Model = UnivariateSpline(data.x, data.y, s=0)
  newdata = DataStructures.xypoint(xgrid.size)
  newdata.x = numpy.copy(xgrid)
  newdata.y = Model(newdata.x)
  
  left = numpy.searchsorted(data.x, (3*xgrid[0]-xgrid[1])/2.0)
  search = numpy.searchsorted
  mean = numpy.mean
  for i in range(xgrid.size-1):
    right = search(data.x, (xgrid[i]+xgrid[i+1])/2.0)
    newdata.y[i] = mean(data.y[left:right])
    left = right
  right = search(data.x, (3*xgrid[-1]-xgrid[-2])/2.0)
  newdata.y[xgrid.size-1] = numpy.mean(data.y[left:right])
  
  return newdata

#This function reduces the resolution by convolving with a gaussian kernel
def ReduceResolution(data,resolution, cont_fcn=None, extend=True, nsigma=8):
  centralwavelength = (data.x[0] + data.x[-1])/2.0
  xspacing = data.x[1] - data.x[0]   #NOTE: this assumes constant x spacing!
  FWHM = centralwavelength/resolution;
  sigma = FWHM/(2.0*numpy.sqrt(2.0*numpy.log(2.0)))
  left = 0
  right = numpy.searchsorted(data.x, 10*sigma)
  x = numpy.arange(0,nsigma*sigma, xspacing)
  gaussian = numpy.exp(-(x-float(nsigma)/2.0*sigma)**2/(2*sigma**2))
  if extend:
    #Extend array to try to remove edge effects (do so circularly)
    before = data.y[-gaussian.size/2+1:]
    after = data.y[:gaussian.size/2]
    extended = numpy.append(numpy.append(before, data.y), after)

    first = data.x[0] - float(int(gaussian.size/2.0+0.5))*xspacing
    last = data.x[-1] + float(int(gaussian.size/2.0+0.5))*xspacing
    x2 = numpy.linspace(first, last, extended.size) 
    
    conv_mode = "valid"

  else:
    extended = data.y.copy()
    x2 = data.x.copy()
    conv_mode = "same"

  newdata = DataStructures.xypoint(data.x.size)
  newdata.x = numpy.copy(data.x)
  if cont_fcn != None:
    cont1 = cont_fcn(newdata.x)
    cont2 = cont_fcn(x2)
    cont1[cont1 < 0.01] = 1
  
    newdata.y = numpy.convolve(extended*cont2, gaussian/gaussian.sum(), mode=conv_mode)/cont1

  else:
    newdata.y = numpy.convolve(extended, gaussian/gaussian.sum(), mode=conv_mode)
    
  return newdata



#Ensure a directory exists. Create it if not
def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)


currentdir = os.getcwd() + "/"
homedir = os.environ["HOME"]
outfiledir = currentdir + "Cross_correlations/"
modeldir = homedir + "/School/Research/Models/Sorted/Stellar/"
gridspacing = "2e-4"
minvel = -1000  #Minimum velocity to output, in km/s
maxvel = 1000  


star_list = ["M2", "M1", "M0", "K9", "K8", "K7", "K6", "K5", "K4", "K3", "K2", "K1","K0", "G9", "G8", "G7", "G6", "G5", "G4", "G3", "G2", "G1", "G0", "F9", "F8", "F7", "F6", "F5", "F4", "F3", "F2", "F1"]
temp_list = [3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000, 5100, 5200, 5225, 5310, 5385, 5460, 5545, 5625, 5700, 5770, 5860, 5940, 6117, 6250, 6395, 6512, 6650, 6775, 6925, 7050, 7185]
model_list = [modeldir + "lte30-3.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte32-3.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte34-3.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte36-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte38-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte40-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte42-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte44-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte46-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte48-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte50-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte51-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte52-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte52-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte53-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte54-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte55-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte55-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte56-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte57-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte58-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte59-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte59-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte61-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte63-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte64-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte65-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte67-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte68-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte69-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte70-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte72-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted"]
    

#This will do the correlation within python/numpy
#The combine keyword decides whether to combine the chips into a master cross-correlation 
#The normalize keyword decides whether to output as correlation power, or as significance
#The sigmaclip keyword decides whether to perform sigma-clipping on each chip before cross-correlating
#The nsigma keyword tells the program how many sigma to clip. This is ignored if sigmaclip = False
#The clip_order keyword tells what order polynomial to fit the flux to during sigma clipping. Ignored if sigmaclip = False
#The models keyword is a list of models to cross-correlate against (either filenames of two-column ascii files, or
#    each entry should be a list with the first element the x points, and the second element the y points
#The segments keyword controls which orders of the data to use, and which parts of them. Can be used to ignore telluric
#    contamination. Can be a string (default) which will use all of the orders, a list of integers which will
#    use all of the orders given in the list, or a dictionary of lists which gives the segments of each order to use.
def PyCorr(filename, combine=True, normalize=False, sigmaclip=False, nsigma=3, clip_order=3, models=model_list, segments="all", save_output=True, outdir=outfiledir, outfilename=None):
  
  ensure_dir(outdir)
  #1: Read in the datafile, if necessary
  if type(filename) == str:
    orders = FitsUtils.MakeXYpoints(filename)
  elif type(filename) == list:
    orders = list(filename)
  else:
    sys.exit("Error! Not sure what to do with input to PyCorr!!")

  makefname = False
  if outfilename == None:
    makefname = True
  
  #2: Interpolate data to a single constant wavelength grid in logspace
  maxsize = 0
  for order in orders:
    if order.size() > maxsize:
      maxsize = order.size()
  data = DataStructures.xypoint(len(orders)*maxsize)
  data.x = numpy.linspace(numpy.log10(orders[-1].x[0]), numpy.log10(orders[0].x[-1]), data.x.size)
  data.y = numpy.ones(data.x.size)
  data.err = numpy.ones(data.x.size)
  data.cont = numpy.ones(data.cont.size)
  firstindex = 1e9
  for i in range(len(orders)):
    order = orders[i]
    order_sections = [[-1, 1e9],]
    #Use this order? Use all of it?
    if type(segments) != str:
      if type(segments) == list:
        for element in segments:
          if element == i+1:
            #Use all of this order
            break
      elif type(segments) == defaultdict or type(segments) == dict:
        try:
          order_sections = segments[i+1]
        except KeyError:
          order_sections = [[-1, -1],]
    
    #Sigma-clipping?
    if sigmaclip:
      done = False
      wave = order.x.copy()
      flux = order.y.copy()
      while not done:
        done = True
        fit = numpy.poly1d(numpy.polyfit(wave, flux, clip_order))
        residuals = flux - fit(wave)
        mean = numpy.mean(residuals)
        std = numpy.std(residuals)
        badindices = numpy.where(numpy.abs(residuals - mean) > nsigma*std)[0]
        flux[badindices] = fit(wave[badindices])
        if badindices.size > 0:
          done = False
      order.y = flux.copy()

    #Interpolate to constant wavelength grid (in log-space)
    FLUX = UnivariateSpline(numpy.log10(order.x), order.y, s=0)
    ERR = UnivariateSpline(numpy.log10(order.x), order.err, s=0)
    CONT = UnivariateSpline(numpy.log10(order.x), order.cont, s=0)
    for section in order_sections:
      left = numpy.searchsorted(order.x, section[0])
      right = numpy.searchsorted(order.x, section[1])
      if right == left:
        continue
      if right > 0:
        right -= 1
      
      left = numpy.searchsorted(data.x, numpy.log10(order.x[left]))
      right = numpy.searchsorted(data.x, numpy.log10(order.x[right]))
      if right > firstindex:
        #Take the average of the two overlapping orders
        data.y[firstindex:right] = (data.y[firstindex:right]/data.cont[firstindex:right] + FLUX(data.x[firstindex:right])/CONT(data.x[firstindex:right]))/2.0
        right = firstindex
      data.y[left:right] = FLUX(data.x[left:right])
      data.err[left:right] = ERR(data.x[left:right])
      data.cont[left:right] = CONT(data.x[left:right])
      firstindex = left

  #3: Begin loop over model spectra
  for i in range(len(models)):
    modelfile = models[i]

    #a: Read in file
    if isinstance(modelfile, str):
      print "******************************\nReading file ", modelfile
      x,y = numpy.loadtxt(modelfile, usecols=(0,1), unpack=True)
      x = x*Units.nm/Units.angstrom
      y = 10**y
    else:
      x = modelfile[0]
      y = modelfile[1]

    left = numpy.searchsorted(x, 2*10**data.x[0] - 10**data.x[-1])
    right = numpy.searchsorted(x, 2*10**data.x[-1] - 10**data.x[0])
    #left = numpy.searchsorted(x, 10**data.x[0])
    #right = numpy.searchsorted(x, 10**data.x[-1])
    model = DataStructures.xypoint(right - left + 1)
    x2 = x[left:right].copy()
    y2 = y[left:right].copy()
    MODEL = UnivariateSpline(x2,y2, s=0)

    #b: Make wavelength spacing constant
    model.x = numpy.linspace(x2[0], x2[-1], right - left + 1)
    model.y = MODEL(model.x)

    #c: Convolve to a resolution of 60000
    model = ReduceResolution(model.copy(), 60000, extend=False)
    MODEL = UnivariateSpline(model.x, model.y, s=0)
    model.x = numpy.linspace(numpy.log10(model.x[0]), numpy.log10(model.x[-1]), model.x.size)
    model.y = MODEL(10**model.x)

    #d: Rebin to the same spacing as the data (but not the same pixels)
    xgrid = numpy.arange(model.x[0], model.x[-1], data.x[1] - data.x[0])
    model = RebinData(model.copy(), xgrid)

    #e: Find continuum by fitting model to a quadratic. Weight high points more
    weights = model.y**2
    done = False
    x2 = 10**model.x
    y2 = model.y
    fitorder = 3
    while not done:
      done = True
      fit = numpy.poly1d(numpy.polyfit(x2 - x2.mean(), y2, fitorder))
      residuals = y2 - fit(x2 - x2.mean())
      mean = numpy.mean(residuals)
      std = numpy.std(residuals)
      badpoints = numpy.where((residuals - mean) < -std)[0]
      if badpoints.size > 0 and x2.size - badpoints.size > 5*fitorder:
        done = False
        x2 = numpy.delete(x2, badpoints)
        y2 = numpy.delete(y2, badpoints)
    model.cont = fit(10**model.x - x2.mean())

    #f: Cross-correlate
    data_rms = numpy.sqrt(numpy.sum((data.y/data.cont-1)**2))
    model_rms = numpy.sqrt(numpy.sum((model.y/model.cont-1)**2))
    left = numpy.searchsorted(model.x, data.x[0])
    right = model.x.size - numpy.searchsorted(model.x, data.x[-1])
    delta = left - right
    print "Cross-correlating..."
   
    #ycorr = numpy.correlate(data.y/data.cont-1.0, model.y/model.cont-1.0, mode="full")
    ycorr = fftconvolve((data.y/data.cont-1.0)[::-1], model.y/model.cont-1.0, mode="full")[::-1]
    xcorr = numpy.arange(ycorr.size)
    lags = xcorr - (model.x.size + data.x.size - delta -1.0)/2.0
    distancePerLag = model.x[1] - model.x[0]
    offsets = -lags*distancePerLag
    velocity = offsets*3e5*numpy.log(10.0)   
    corr = DataStructures.xypoint(velocity.size)
    corr.x = velocity[::-1]
    corr.y = ycorr[::-1]/(data_rms*model_rms)
    #My version at home has a bug in numpy.correlate, reversing ycorr
    #BUG FIXED IN THE PYTHON VERSION I HAVE FOR LINUX MINT 13
    #if "linux" in sys.platform:
    #  corr.y = corr.y[::-1]
        
    #Fit low-order polynomal to cross-correlation
    left = numpy.searchsorted(corr.x, minvel)
    right = numpy.searchsorted(corr.x, maxvel)
    vel = corr.x[left:right]
    corr = corr.y[left:right]
    fit = numpy.poly1d(numpy.polyfit(vel, corr, 2))
        
    #Adjust correlation by fit
    corr = corr - fit(vel)
    if normalize:
      mean = numpy.mean(corr)
      std = numpy.std(corr)
      corr = (corr - mean)/std

    #j: Finally, output or return
    if save_output:
      if makefname:
        outfilename = outdir + filename.split("/")[-1] + "." + star
      print "Outputting to ", outfilename, "\n"
      numpy.savetxt(outfilename, numpy.transpose((vel, corr)), fmt="%.8g")
    else:
      return vel, corr
    



if __name__ == "__main__":
  if len(sys.argv) > 1:
    for fname in sys.argv[1:]:
      PyCorr(fname) #, combine=False, sigmaclip=False)
  
