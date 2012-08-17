"""
   Function to fit tellurics in the optical regime, given a McDonald-style reduced echelle spectrum
   Molecules to fit are H2O and O2
   Also fit wavelength (should already have a good first guess from the reduction),
     continuum (same), and resolution (should be ~60000)
"""


import pylab
import pyfits
import numpy
import sys
import os
import subprocess
import scipy
from scipy.interpolate import interp1d, UnivariateSpline
from scipy.optimize import leastsq, brute, fmin
from scipy.linalg import svd, diagsvd
from scipy import mat
from collections import defaultdict
import MakeModel
import DataStructures
import FitsUtils
import Units
import RotBroad


homedir = os.environ['HOME']
TelluricModelDir = homedir + "/School/Research/lblrtm/run_examples/MyModel/"
LineListFile = homedir + "/School/Research/McDonaldData/MyLinelist.dat"
ContinuumFile = homedir + "/School/Research/Models/PlanetFinder/src/CRIRES/ContinuumRegions.dat"
Bstarfiles = {#2.25: homedir + "/School/Research/McDonaldData/BstarModels/BG19000g225v2.vis.7",
              #2.50: homedir + "/School/Research/McDonaldData/BstarModels/BG19000g250v2.vis.7",
              #2.75: homedir + "/School/Research/McDonaldData/BstarModels/BG19000g275v2.vis.7",
              #3.00: homedir + "/School/Research/McDonaldData/BstarModels/BG19000g300v2.vis.7",
              3.25: homedir + "/School/Research/McDonaldData/BstarModels/BG19000g325v2.vis.7",
              3.50: homedir + "/School/Research/McDonaldData/BstarModels/BG19000g350v2.vis.7",
              3.75: homedir + "/School/Research/McDonaldData/BstarModels/BG19000g375v2.vis.7",
              4.00: homedir + "/School/Research/McDonaldData/BstarModels/BG19000g400v2.vis.7",
              4.25: homedir + "/School/Research/McDonaldData/BstarModels/BG19000g425v2.vis.7",
              4.50: homedir + "/School/Research/McDonaldData/BstarModels/BG19000g450v2.vis.7",
              4.75: homedir + "/School/Research/McDonaldData/BstarModels/BG19000g475v2.vis.7"}

#Define some bounds. Bound functions are defined below
o2_bounds = [1.8e5, 2.5e5]
humidity_bounds = [0,100]
resolution_bounds = [35000,90000]
angle_bounds = [0.0, 80.0]
vsini_bounds = [1.0,500*Units.cm/Units.km]
primary_rv_bounds = [-0.1*Units.c, 0.1*Units.c]
searchgrid = ((1.8e5,2.5e5,1e4),(20,90, 5))


#Read in the parent star file as a global variable (don't rotationally broaden yet!), as a function of surface gravity
bstars = defaultdict(DataStructures.xypoint)
print "Reading in primary star models"
for key in Bstarfiles:
  filename = Bstarfiles[key]
  print "Reading model ", filename
  bstars[key] = RotBroad.ReadFile(filename)

class ContinuumSegments:
  def __init__(self,size):
    self.low = numpy.zeros(size)
    self.high = numpy.zeros(size)

def Main(filename, humidity=None, resolution=None, angle=None, ch4=None, co=None, o2=None):
  
  #Get some information from the fits header
  hdulist = pyfits.open(filename)
  header = hdulist[0].header
  orders = FitsUtils.MakeXYpoints(header, hdulist[0].data)
  
  #Resolution:
  if resolution == None:
    resolution = 60000
    
  #Angle:
  if angle == None:
    angle = float(header["zd"])
  angle_bounds = [angle-2, angle+2] #Do not allow zenith angle to change much
  
  #Humidity:
  if humidity == None:
    humidity = 40.0
  
  #Always use fits header for the temperature and pressure
  temperature = (64 - 32)*5./9. + 273.15
  pressure_start = 23.55*Units.hPa/Units.inch_Hg
  pressure_end = 23.58*Units.hPa/Units.inch_Hg
  pressure = (pressure_start + pressure_end)/2.0
  
  wave_start = 400
  wave_end = 1100
  wavenum_start = int(1.0e7/wave_end)
  wavenum_end = int(1.0e7/wave_start+1)
  
  xspacing = 1e-4    #Change this to change the interpolation spacing
  #Generate interpolated chip data
  for order in orders:
    order.cont = numpy.ones(order.x.size)
    order.err[order.y > 0.0] = numpy.sqrt(order.y[order.y > 0.0])
    order.err[order.y <= 0.0] = 1e9
    
    #make sure there are no zeros (fancy indexing)
    order.err[order.err <= 0] = 1e10
  hdulist.close()
  
  #Read in line list:
  linelist = numpy.loadtxt(LineListFile)
  
  #Read in continuum database:
  segments = ReadContinuumDatabase(ContinuumFile)
  
  #Set up parameter lists
  if ch4 == None:
    ch4 = 1.6   #CH4 abundance in ppm
  co2 = 368.5 
  if co == None:
    co = 0.14
  o3 = 4e-2
  o2 = 2.12e5
  vsini = 140*Units.cm/Units.km
  
  continuum_fit_order = [5,3,3,3]
  const_pars = [temperature, pressure, co2, o3, wavenum_start, wavenum_end, resolution, angle, ch4, co, humidity, o2, vsini]
  pars = [humidity, o2, angle]
  
  #Make outfilename from the input fits file
  outfilename = "Corrected_" + filename[3:].split("-")[0] + ".dat"
  print "outputting to ", outfilename
  outfile = open(outfilename, "w")
  debug = False
  if (debug):
    ErrorFunctionBrute = lambda pars, chip, const_pars, linelist, contlist: numpy.sum(ErrorFunction(pars, chip, const_pars, linelist, contlist)**2)
  for i in range(3,len(orders)-2):
    #ErrorFunction(pars, chips[i], const_pars, linelist, segments)
    order = orders[i]
    
    print "Fitting order ", i+1, "with size ", order.x.size
    
    if not debug:
      #const_pars[8] = continuum_fit_order[i]
      model = FitFunction(order, pars, const_pars)
      model = MakeModel.ReduceResolution(model.copy(), resolution)
      model = MakeModel.RebinData(model.copy(), order.x.copy())

      order = FitContinuum3(order, model, 4)

      order = CCImprove(order, model)
      fitout = leastsq(ErrorFunction, pars, args=(order, const_pars, linelist, segments), full_output=True, epsfcn = 0.0005, maxfev=1000)
      fitpars = fitout[0]
      pars = fitpars
    else:
      fitout = brute(ErrorFunctionBrute, searchgrid, args=(order, const_pars, linelist, segments))
      fitpars = pars
      
    outfile2 = open("chisq_summary.dat", 'a')
    outfile2.write("\n\n\n\n")
    outfile2.close()
    
    print fitpars
    print "Done fitting chip ", i, "\n\n\n"
    model = FitFunction(order, fitpars, const_pars)
    model_original = model.copy()
  
    #convolve to approximate resolution (will fit later)
    resolution = const_pars[6]
    Continuum = UnivariateSpline(order.x, order.cont, s=0)
    model = MakeModel.ReduceResolution(model, resolution, Continuum)
    model = MakeModel.RebinData(model, order.x)
    
    #Fit Continuum of the chip, using the model
    #chips[i] = FitContinuum2(chips[i+2],model,segments)
    order = FitContinuum3(order, model,4)
    #if i == 0:
    #  chips[i] = FitContinuum(chips[i], model, condition=0.99)
    #else:
    #  chips[i] = FitContinuum(chips[i], model, order=3, condition=0.99)
  
    #Fit model wavelength to the chip
    #model = FitWavelength(chips[i], model,linelist)
    modelfcn, mean = FitWavelength2(order, model, linelist, const_pars[12], debug=False)
    model_original.x = modelfcn(model_original.x - mean)    

    #Fit resolution:
    model, resolution = FitResolution(order, model_original, resolution)
    Model = UnivariateSpline(model.x, model.y,s=0)
    
    #Estimate errors:
    model2 = Model(order.x)
    residuals = model2 - order.y
    std = numpy.std(residuals)
    covariance = numpy.sqrt(fitout[1]*std)

    #For very deep lines, just set the line cores equal to the model (will destroy any info in the line, but will prevent huge residuals)
    indices = numpy.where(model2 < 0.05)[0]
    order.y = order.y.astype('float64')
    if indices.size > 0:
      order.y[indices] = model2[indices]*order.cont[indices]

    #Output
    if not debug:
      outfile.write("#Temperature: " + str(const_pars[0]) + "\n")
      outfile.write("#Pressure: " + str(const_pars[1]) + "\n")
      outfile.write("#Humidity: " + str(fitpars[0]) + " +/- " + str(covariance[0][0]) + "\n")
      outfile.write("#O2: " + str(fitpars[1]) + " +/- " + str(covariance[1][1]) + "\n")
      outfile.write("#CH4: " + str(const_pars[8]) + "\n")
      outfile.write("#CO: " + str(const_pars[9]) + "\n")
      outfile.write("#Angle: " + str(fitpars[2]) + " +/- " + str(covariance[2][2]) + "\n")
      outfile.write("#Resolution: " + str(resolution) + "\n")
      outfile.write("#Vsini: " + str(const_pars[12]) + "\n")
      outfile.write("#Convergence message: " + fitout[3] + "\n")
      outfile.write("#Convergence code: " + str(fitout[4]) + "\n")
      for j in range(order.x.size):
        outfile.write("%.15f" %order.x[j] + "\t1.0\t%.15f" %(order.y[j]/model2[j]) + "\t%.15f" %order.y[j] + "\t%.15f" %model2[j] + "\t1.0\t%.15f" %order.err[j] + "\t%.15f" %order.cont[j] + "\n")
      outfile.write("\n\n\n")

      orders[i].y /= model2
      FitsUtils.OutputFitsFile(filename, orders)
    
    #pylab.plot(chips[i].x, chips[i].y/chips[i].cont, label="data")
    #pylab.plot(chips[i].x, model2, label="model")
    
  outfile.close()
  #pylab.legend()
  #pylab.show()
  
#This function reads the continuum database, and returns the segments
def ReadContinuumDatabase(filename):
  low,high = numpy.loadtxt(filename,usecols=(0,1),unpack=True)
  segments = ContinuumSegments(low.size)
  segments.low = low
  segments.high = high
  return segments  
  
########################################################################
#This function will Generate a telluric model with the given parameters#
########################################################################
def FitFunction(order, pars, const_pars):
  temperature = const_pars[0]
  pressure = const_pars[1]
  co2 = const_pars[2]
  o3 = const_pars[3]
  wavenum_start = const_pars[4]
  wavenum_end = const_pars[5]
  angle = const_pars[7]
  resolution = const_pars[6]
  ch4 = const_pars[7]
  co = const_pars[8]
  humidity = pars[0]
  o2 = pars[1]
  angle = pars[2]
  vsini = const_pars[12]
  plotflg = False

  wave_start = order.x[0] - 10.
  wave_end = order.x[-1] + 10.
  wavenum_start = int(1.0e7/wave_end)
  wavenum_end = int(1.0e7/wave_start+1)
  
  #Make sure certain variables are positive
  if o2 < 0:
    co = 0
    print "\nWarning! O2 was set to be negative. Resetting to zero before generating model!\n\n"
    #plotflg = True
  if ch4 < 0:
    ch4 = 0
    print "\nWarning! CH4 was set to be negative. Resetting to zero before generating model!\n\n"
    #plotflg = True
  if humidity < 0:
    humidity = 0
    print "\nWarning! Humidity was set to be negative. Resetting to zero before generating model!\n\n"
    #plotflg = True
  if angle < 0:
    angle = -angle
    print "\nWarning! Angle was set to be negative. Resetting to a positive value before generating model!\n\n"
    #plotflg = True

  #Generate the model:
  model = MakeModel.Main(pressure, temperature, humidity, wavenum_start, wavenum_end, angle, co2, o3, ch4, co, o2, order.x, resolution)
  if "FullSpectrum.freq" in os.listdir(TelluricModelDir):
    cmd = "rm "+ TelluricModelDir + "FullSpectrum.freq"
    command = subprocess.check_call(cmd, shell=True)

  if plotflg:
    pylab.plot(order.x, order.y)
    pylab.plot(model.x, model.y)
    pylab.show()
  return model
  
  
def ErrorFunction(pars, order, const_pars, linelist, contlist):
  model = FitFunction(order.copy(), pars, const_pars)
   
  model_original = model.copy()
  plotflg = False
  #if (pars[1] <= 0):
  #  plotflg = True
  
  #Reduce to initial guess resolution
  Continuum = UnivariateSpline(order.x, order.cont, s=0)
  resolution = const_pars[6]
  if (resolution - 10 < resolution_bounds[0] or resolution+10 > resolution_bounds[1]):
    resolution = 60000
  model = MakeModel.ReduceResolution(model.copy(), resolution, Continuum)
  model = MakeModel.RebinData(model.copy(), order.x.copy())

  PRIMARY_STAR = FitPrimary(bstars, model, order.copy(), vsini=150*Units.cm/Units.km, z=10.0)

  #Fit Continuum of the chip, using the model
  #chip = FitContinuum2(chip,model,contlist)
  model2 = model.copy()
  model2.y *= PRIMARY_STAR(model2.x)
  order.cont = numpy.ones(order.x.size)
  #order = FitContinuum3(order, model2, 4)
  #fit_order = const_pars[8]
  #chip = FitContinuum(chip, model, condition=0.99, tol=3, order=fit_order)
  
  #Fit model wavelength to the chip
  #model = FitWavelength(chip,model,linelist)
  
  #order = CCImprove(order, model)
  modelfcn, mean = FitWavelength2(order.copy(), model.copy(), linelist, const_pars[12], primary_fcn = PRIMARY_STAR)
  model.x = modelfcn(model.x - mean)
  model_original.x = modelfcn(model_original.x - mean)
  model2 = model_original.copy()
  model2.y *= PRIMARY_STAR(model2.x)

  #Fit resolution
  model, resolution = FitResolution(order.copy(), model2, resolution, plotflg)
  const_pars[6] = resolution
  
  weights = 1.0/order.err
  weights = weights/weights.sum()
  return_array = ((order.y  - order.cont*model.y)**2*weights + bound(humidity_bounds,pars[0]) + 
                                          bound(o2_bounds,pars[1]) + 
  					  bound(angle_bounds, pars[2]))
  print "X^2 = ", numpy.sum(return_array)/float(weights.size)
  outfile = open("chisq_summary.dat", 'a')
  outfile.write(str(pars[0])+"\t"+str(pars[1])+"\t"+str(resolution)+"\t"+str(numpy.sum(return_array)/float(weights.size))+"\n")
  outfile.close()
  
  return return_array
  
  
#Define bounding functions:
# lower bound:            lbound(boundary_value, parameter)
# upper bound:            ubound(boundary_value, parameter)
# lower and upper bounds: bound([low, high], parameter)
# fixed parameter:        fixed(fixed_value, parameter)
lbound = lambda p, x: 1e4*numpy.sqrt(p-x) + 1e-3*(p-x) if (x<p) else 0
ubound = lambda p, x: 1e4*numpy.sqrt(x-p) + 1e-3*(x-p) if (x>p) else 0
bound  = lambda p, x: lbound(p[0],x) + ubound(p[1],x)
fixed  = lambda p, x: bound((p,p), x)


"""
Function to fit the radial velocity (z) and rotational velocity (vsini) of the primary star
"""
def FitPrimary(bstar_dict, model, order, vsini=150*Units.cm/Units.km, alpha=0.5, intervalsize=50.0, z=0.0):
  vsini = 100*Units.cm/Units.km
  z = 10*Units.cm/Units.km/Units.c
  pars = [vsini, z]
  order.y /= model.y
  const_pars = [bstar_dict, order, alpha, intervalsize]
  pars, success = leastsq(primary_errfunc, pars, args=const_pars, epsfcn=1)
  print "Best fit vsini: ", pars[0]
  print "Best fit radial velocity: ", pars[1]*Units.c*Units.km/Units.cm, "km/s"

  broadened, logg = FindBestGravity(const_pars[0], const_pars[1], pars[0], const_pars[2], const_pars[3], pars[1])
  fcn = UnivariateSpline(broadened.x, broadened.y/broadened.cont, s=0)

  return fcn


def FindBestGravity(bstar_dict, order, vsini, alpha, intervalsize, z):
  best_chisq = 1e18
  best_star = bstar_dict[bstar_dict.keys()[0]].copy()
  best_logg = 0
  search = numpy.searchsorted
  broaden = RotBroad.Broaden
  add = numpy.sum
  for logg in sorted(bstar_dict.keys()):
    left = search(bstar_dict[logg].x*(1+z), 2*order.x[0] - order.x[-1])
    right = search(bstar_dict[logg].x*(1+z), 2*order.x[-1] - order.x[0])
    bstar = DataStructures.xypoint(right - left + 1)
    bstar.x = bstar_dict[logg].x[left:right]*(1+z)
    bstar.y = bstar_dict[logg].y[left:right]
    bstar.cont = bstar_dict[logg].cont[left:right]
    broadened = broaden(bstar, vsini=vsini, alpha=alpha, intervalsize=intervalsize)
    fcn = UnivariateSpline(broadened.x, broadened.y/broadened.cont, s=0)

    chisq = add((order.y - fcn(order.x))**2/order.err**2)
    print "chisq = ", chisq
    if chisq < best_chisq:
      best_chisq = chisq
      best_star = broadened.copy()
      best_logg = logg

  return best_star, best_logg

  
    
"""
Error function for the leastsq call in the 'FitPrimary' function
"""
def primary_errfunc(pars, const_pars):
  order = const_pars[1]
  vsini = pars[0]
  RV = pars[1]
  print "vsini = ", vsini*Units.km/Units.cm, "km/s"
  print "RV = ", RV*Units.c*Units.km/Units.cm, "km/s"
  if vsini < vsini_bounds[0]:
    vsini = 100*Units.cm/Units.km
  if vsini > vsini_bounds[-1]:
    vsini = vsini_bounds[-1]
  broadened, logg = FindBestGravity(const_pars[0], order, vsini, const_pars[2], const_pars[3], RV)
  fcn = UnivariateSpline(broadened.x, broadened.y/broadened.cont, s=0)

  print "log(g) = ", logg
  print "X^2 = ", numpy.sum((order.y - fcn(order.x))**2/order.err**2), '\n'
  return (order.y - fcn(order.x))**2/order.err**2 + bound(vsini_bounds,pars[0]) + bound(primary_rv_bounds,pars[1])
  
  


#Function to fit the resolution
def FitResolution(data, model, resolution=75000, plotflg = False):
  ####resolution is the initial guess####
  #Interpolate to a constant wavelength grid
  print "Fitting resolution"
  if plotflg:
    #pylab.plot(data.x, data.y)
    pylab.plot(model.x, model.y)
    pylab.show()
  ModelFcn = UnivariateSpline(model.x, model.y, s=0)
  newmodel = DataStructures.xypoint(model.x.size*2)
  newmodel.x = numpy.linspace(model.x[0], model.x[-1], model.x.size*2)
  newmodel.y = ModelFcn(newmodel.x)

  Continuum = UnivariateSpline(data.x, data.cont, s=0)

  #errfunc = lambda R, data, model: (data.y - MakeModel.ReduceResolutionAndRebinData(model,R,data.x).y - bound(resolution_bounds, R))
  #Do a brute force grid search first, then refine with Levenberg-Marquardt
  searchgrid = (resolution_bounds[0], resolution_bounds[1], 5000)
  #resolution = brute(ResolutionFitErrorBrute,(searchgrid,), args=(data,newmodel,Continuum)) #, finish=fmin)
  resolution = brute(ResolutionFitErrorBrute,(searchgrid,), args=(data,newmodel))
  resolution, success = leastsq(ResolutionFitError, resolution, args=(data, newmodel), epsfcn=10)
  #resolution, success = leastsq(ResolutionFitError, resolution, args=(data, newmodel, Continuum), epsfcn=10)
  #resolution, success = leastsq(ResolutionFitError, resolution, args=(data,newmodel))
  print "Optimal resolution found at R = ", float(resolution)
  #newmodel = MakeModel.ReduceResolution(newmodel, float(resolution), Continuum)
  newmodel = MakeModel.ReduceResolution(newmodel, float(resolution))
  return MakeModel.RebinData(newmodel, data.x), float(resolution)
  
  
def ResolutionFitError(resolution, data, model, cont_fcn=None):
  newmodel = MakeModel.ReduceResolution(model, resolution, cont_fcn)
  newmodel = MakeModel.RebinData(newmodel, data.x)
  weights = 1.0/data.err
  weights = weights/weights.sum()
  returnvec = (data.y - data.cont*newmodel.y)**2*weights + bound(resolution_bounds, resolution)
  print "Resolution-fitting X^2 = ", numpy.sum(returnvec)/float(weights.size), "at R = ", resolution
  if numpy.isnan(numpy.sum(returnvec**2)):
    #sys.exit("NaN found in ResolutionFitError! Exiting...")
    print "Error! NaN found in ResolutionFitError!"
    outfile=open("ResolutionFitError.log", "a")
    for i in range(data.y.size):
      outfile.write("%.10g\t" %data.x[i] + "%.10g\t" %data.y[i] + "%.10g\t" %data.cont[i] + "%.10g\t" %newmodel.x[i] + "%.10g\n" %newmodel.y[i])
    outfile.write("\n\n\n\n")
    outfile.close()
  return returnvec
 
def ResolutionFitErrorBrute(resolution, data, model, cont_fcn=None):
  return numpy.sum(ResolutionFitError(resolution, data, model, cont_fcn))


#Fits the broadening profile using singular value decomposition
#oversampling is the oversampling factor to use before doing the SVD
#m is the size of the broadening function, in oversampled units
#dimension is the number of eigenvalues to keep in the broadening function. (Keeping too many starts fitting noise)
def Broaden(data, model, oversampling = 5, m = 101, dimension = 15):
  n = data.x.size*oversampling
  
  #resample data
  Spectrum = UnivariateSpline(data.x, data.y/data.cont, s=0)
  Model = UnivariateSpline(model.x, model.y, s=0)
  xnew = numpy.linspace(data.x[0], data.x[-1], n)
  ynew = Spectrum(xnew)
  #model_new = Model(xnew)
  model_new = MakeModel.RebinData(model, xnew).y

  #Make 'design matrix'
  design = numpy.zeros((n-m,m))
  for j in range(m):
    for i in range(m/2,n-m/2-1):
      design[i-m/2,j] = model_new[i-j+m/2]
  design = mat(design)
  #Do Singular Value Decomposition
  U,W,V_t = svd(design, full_matrices=False)
  
  #Invert matrices:
  #   U, V are orthonormal, so inversion is just their transposes
  #   W is a diagonal matrix, so its inverse is 1/W
  W1 = 1.0/W
  U_t = numpy.transpose(U)
  V = numpy.transpose(V_t)
  
  #Remove the smaller values of W
  W1[dimension:] = 0
  W2 = diagsvd(W1,m,m)
  #Solve for the broadening function
  spec = numpy.transpose(mat(ynew[m/2:n-m/2-1]))
  temp = numpy.dot(U_t, spec)
  temp = numpy.dot(W2,temp)
  Broadening = numpy.dot(V,temp)
  #Convolve model with this function
  spacing = xnew[2] - xnew[1]
  xnew = numpy.arange(model.x[0], model.x[-1], spacing)
  model_new = Model(xnew)
  Broadening = numpy.array(Broadening)[...,0]
  model.x = xnew
  Broadened = UnivariateSpline(xnew, numpy.convolve(model_new,Broadening, mode="same"),s=0)
  model.y = Broadened(model.x)

  return model


"""
Do a cross-correlation first, to get the wavelength solution close
"""
def CCImprove(order, model):
  ycorr = scipy.correlate(order.y-1.0, model.y-1.0, mode="full")
  xcorr = numpy.arange(ycorr.size)
  maxindex = ycorr.argmax()
  lags = xcorr - (order.y.size-1)
  distancePerLag = (order.x[-1] - order.x[0])/float(order.x.size)
  offsets = -lags*distancePerLag
  print "maximum offset: ", offsets[maxindex], " nm"

  if numpy.abs(offsets[maxindex]) < 0.2:
    #Apply offset
    order.x = order.x + offsets[maxindex]
  return order


"""
Function to Fit the wavelength solution, using a bunch of telluric lines
This assumes that we are already quite close to the correct solution
"""
def FitWavelength(data, model, linelist, tol=2e-2, FWHM=0.05, debug=False):
  print "Fitting wavelength"
  sigma = FWHM/2.65	#approximate relation, but good enough for initial guess
  
  #Fitting function: gaussian on a linear background:
  fitfunc = lambda p,x: p[0] + p[1]*x + p[2]*numpy.exp(-(x-p[3])**2/(2*p[4]**2))
  errfunc = lambda p,x,y: y - fitfunc(p,x)
  fitfunc2 = lambda p,x: p[0] + p[1]*numpy.exp(-(x-p[2])**2/(2*p[3]**2))
  errfunc2 = lambda p,x,y,err: (y - fitfunc2(p,x))/(err)

  data_points = []
  model_points = []
  data_rms = []
  model_rms = []
  
  if debug:
    print "Plotting in debug mode"
    pylab.plot(data.x, data.y)
    pylab.plot(model.x, model.y)
  
  #Begin loop over the lines
  for line in linelist:
    if line-tol-FWHM/2.0 > data.x[0] and line+tol+FWHM/2.0 < data.x[-1]:
      #This line falls on the chip. Use it
      #First, find the lowest point within tol of the given line position
      left = numpy.searchsorted(data.x, line - tol)
      right = numpy.searchsorted(data.x, line + tol)
      minindex = data.y[left:right].argmin() + left
      
      #Set up guess parameters for line fit:
      mean = data.x[minindex]
      left2 = numpy.searchsorted(data.x, mean - FWHM/2.0)
      right2 = numpy.searchsorted(data.x, mean + FWHM/2.0)
      slope = (data.y[right2] - data.y[left2])/(data.x[right2] - data.x[left2])
      a = data.y[right2] - slope*data.x[right2]
      b = slope
      amplitude = (a + b*mean) - data.y[minindex]
      #pars = [a,b,-amplitude,mean,sigma]
      pars=[1,-amplitude, mean, sigma]
      print "Number of points in fit: ", right2 - left2 + 1

      #Do the fit
      #pars, success = leastsq(errfunc, pars, args=(data.x[left2:right2], data.y[left2:right2]))
      pars, success = leastsq(errfunc2, pars, args=(data.x[left2:right2], data.y[left2:right2]/data.cont[left2:right2], data.err[left2:right2]))      
      data_fit_failed = True
      if success < 5:
        #Determine rms error
        resid = data.y[left2:right2] - fitfunc2(pars, data.x[left2:right2])
        data_rms.append(numpy.sqrt(numpy.sum(resid**2)))
        
        #Save the mean value:
        mean = pars[2]
        data_points.append(mean)
        data_fit_failed = False
        if debug:
          print "Plotting in debug mode"
          pylab.plot(data.x[left2:right2], fitfunc2(pars, data.x[left2:right2]))
          pylab.plot(mean, fitfunc2(pars, mean), 'ro')
      elif debug:
        print "Data fit bad"
      ################################################################
      #  Now, do the same thing for the model
      ################################################################
      minindex = model.y[left:right].argmin() + left
      
      #Set up guess parameters for line fit:
      mean = model.x[minindex]
      left2 = numpy.searchsorted(model.x, mean - FWHM/2.0)
      right2 = numpy.searchsorted(model.x, mean + FWHM/2.0)
      slope = (model.y[right2] - model.y[left2])/(model.x[right2] - model.x[left2])
      a = model.y[right2] - slope*model.x[right2]
      b = slope
      amplitude = (a + b*mean) - model.y[minindex]
      #pars = [a,b,-amplitude,mean,sigma]
      pars = [1,-amplitude, mean, sigma]
      
      #Do the fit
      #pars, success = leastsq(errfunc, pars, args=(model.x[left2:right2], model.y[left2:right2]))
      pars, success = leastsq(errfunc2, pars, args=(model.x[left2:right2], model.y[left2:right2]/model.y[left2:right2], numpy.ones(model.y[left2:right2].size)))
      
      if success < 5 and not data_fit_failed:
        #Determine rms error
        resid = model.y[left2:right2] - fitfunc2(pars, model.x[left2:right2])
        model_rms.append(numpy.sqrt(numpy.sum(resid**2)))

        #Save the mean value:
        mean=pars[2]
        model_points.append(mean)
        if debug:
          print "Plotting in debug mode"
          pylab.plot(model.x[left2:right2], fitfunc2(pars, model.x[left2:right2]))
          pylab.plot(mean, fitfunc2(pars,mean), 'bo')
      elif success >=5 and not data_fit_failed:
        data_points.pop()
      elif debug:
        print "Model fit bad"
        
  #Output figure:
  figs = os.listdir("FitPictures")
  j = 0
  for fname in figs:
   if "lines-" in fname:
     i = fname.split("lines-")[-1].split(".")[0]
     if int(i) > j:
       j = int(i)
  linefig = "FitPictures/lines-" + str(j+1) + ".png"
  fitfig = "FitPictures/fit-" + str(j+1) + ".png"
  if (debug):
    pylab.savefig(linefig, dpi=600)
    pylab.show()
    pylab.cla()
  #Remove points with too large an rms
  done = False
  order = 3
  while not done:
    done = True
    pars = numpy.polyfit(model_points, data_points, order)
    fit = numpy.poly1d(pars)
    residuals = list(data_points - fit(model_points))
    mean = numpy.mean(residuals)
    std = numpy.std(residuals)
    for i in range(len(data_points)-1,-1,-1):
      if numpy.abs(residuals[i]) > mean + 2.5*std:
        if (debug):
          print "Removing point near ", data_points[i]
          pylab.plot(model_points[i], data_points[i], 'bx')
        data_points.pop(i)
        model_points.pop(i)
        data_rms.pop(i)
        model_rms.pop(i)
        residuals.pop(i)
        done = False
  if debug:  
    for i in range(len(data_points)):
      pylab.plot(model_points[i], data_points[i], 'ro')
  
  #Fit datapoints and modelpoints to a quadratic
  #pars = numpy.polyfit(model_points, data_points, order)
  #fit = numpy.poly1d(pars)
  fit = UnivariateSpline(model_points, data_points,k=2,s=0)
  if (debug):
    pylab.plot(model.x, fit(model.x))
    pylab.savefig(fitfig, dpi=600)
    pylab.show()
    residuals = data_points - fit(model_points)
    pylab.plot(data_points, residuals, 'ro')
    pylab.show()
    pylab.cla()
  
  model.x = fit(model.x)
  return model


#Second wavelength-fitting function that just shifts lines, instead of fitting them to gaussians
def WavelengthErrorFunction(shift, data, model):
  modelfcn = UnivariateSpline(model.x, model.y, s=0)
  weight = 1.0/numpy.sqrt(data.y)
  weight[weight < 0.01] = 0.0
  newmodel = modelfcn(model.x + float(shift))
  if shift < 0:
    newmodel[model.x - float(shift) < model.x[0]] = 0
  else:
    newmodel[model.x - float(shift) > model.x[-1]] = 0
  returnvec = (data.y - newmodel)**2*weight
  return returnvec

def FitWavelength2(order, telluric, linelist, parent_vsini, primary_fcn=None, tol=0.05, oversampling = 4, debug=False, max_change=2.0):
  print "Fitting wavelength"
  old = []
  new = []

  if primary_fcn == None:
    primary_fcn = numpy.poly1d([1.0,])

  #Interpolate to finer spacing
  DATA_FCN = UnivariateSpline(order.x, order.y, s=0)
  CONT_FCN = UnivariateSpline(order.x, order.cont, s=0)
  MODEL_FCN = UnivariateSpline(telluric.x, telluric.y, s=0)
  data = DataStructures.xypoint(order.x.size*oversampling)
  data.x = numpy.linspace(order.x[0], order.x[-1], order.x.size*oversampling)
  data.y = DATA_FCN(data.x)
  data.cont = CONT_FCN(data.x)
  model = DataStructures.xypoint(data.x.size)
  model.x = numpy.copy(data.x)
  model.y = MODEL_FCN(model.x)*primary_fcn(model.x)
  
  
  #Begin loop over the lines
  for line in linelist:
    if line-tol > data.x[0] and line+tol < data.x[-1]:
      #Find line in the model
      left = numpy.searchsorted(model.x, line - tol)
      right = numpy.searchsorted(model.x, line + tol)
      minindex = model.y[left:right].argmin() + left

      mean = model.x[minindex]
      left2 = numpy.searchsorted(model.x, mean - tol*2)
      right2 = numpy.searchsorted(model.x, mean + tol*2)

      argmodel = DataStructures.xypoint(right2 - left2)
      argmodel.x = numpy.copy(model.x[left2:right2])
      argmodel.y = numpy.copy(model.y[left2:right2])

      #Do the same for the data
      left = numpy.searchsorted(data.x, line - tol)
      right = numpy.searchsorted(data.x, line + tol)
      minindex = data.y[left:right].argmin() + left

      mean = data.x[minindex]

      argdata = DataStructures.xypoint(right2 - left2)
      argdata.x = numpy.copy(data.x[left2:right2])
      argdata.y = numpy.copy(data.y[left2:right2]/data.cont[left2:right2])

      #Fit argdata to gaussian:
      cont = 1.0
      depth = cont - argdata.y[argdata.y.size/2]
      mu = argdata.x[argdata.x.size/2]
      sig = 0.025
      params = [cont, depth, mu, sig]
      params,success = leastsq(GaussianErrorFunction, params, args=(argdata.x, argdata.y))
      
      mean = params[2]
      #Do a cross-correlation first, to get the wavelength solution close
      ycorr = scipy.correlate(argdata.y-1.0, argmodel.y-1.0, mode="full")
      xcorr = numpy.arange(ycorr.size)
      maxindex = ycorr.argmax()
      lags = xcorr - (argdata.x.size-1)
      distancePerLag = (argdata.x[-1] - argdata.x[0])/float(argdata.x.size)
      offsets = -lags*distancePerLag
      shift = offsets[maxindex]
      shift, success = leastsq(WavelengthErrorFunction, shift, args=(argdata, argmodel))
      if (debug):
        print argdata.x[0], argdata.x[-1], argdata.x.size
        print "wave: ", mean, "\tshift: ", shift, "\tsuccess = ", success
        pylab.plot(model.x[left:right]-shift, model.y[left:right], 'g-')
        pylab.plot(argmodel.x, argmodel.y, 'r-')
        pylab.plot(argdata.x, argdata.y, 'k-')
      if (success < 5):
        old.append(mean)
        new.append(mean - float(shift))
  if debug:
    pylab.show()
    pylab.plot(old, new, 'ro')
    pylab.show()
  #fit = UnivariateSpline(old, new, k=1, s=0)
  #Iteratively fit to a cubic with sigma-clipping
  order = 3
  done = False
  print "Number of wavelength points: ", len(old), "\t", len(new)
  if len(old) < order:
    fit = numpy.poly1d((1,0))
    return fit, 0.0 

  while not done:
    done = True
    mean = numpy.mean(old)
    fit = numpy.poly1d(numpy.polyfit(old - mean, new, order))
    residuals = fit(old - mean) - new
    std = numpy.std(residuals)
    #if debug:
    #  pylab.plot(old, residuals, 'ro')
    #  pylab.plot(old, std*numpy.ones(len(old)))
    #  pylab.show()
    badindices = numpy.where(numpy.logical_or(residuals > 2*std, residuals < -2*std))[0]
    for badindex in badindices[::-1]:
      del old[badindex]
      del new[badindex]
      done = False

  #Check if the function changed things by too much
  difference = numpy.abs(telluric.x - fit(telluric.x - mean))
  if numpy.any(difference > max_change):
    fit = numpy.poly1d((1,0))
    mean = 0.0

  if debug:
    pylab.plot(old, fit(old - mean) - new, 'ro')
    pylab.show()

  return fit, mean
  

#Gaussian absorption line
def GaussianFitFunction(x,params):
  cont = params[0]
  depth = params[1]
  mu = params[2]
  sig = params[3]
  return cont - depth*numpy.exp(-(x-mu)**2/(2*sig**2))

#Returns the residuals between the fit from above and the actual values
def GaussianErrorFunction(params, x, y):
  return GaussianFitFunction(x,params) - y


#This function will fit the continuum in the regions given
def FitContinuum(order,model,contlist):
  print "Fitting continuum"
  wave = []
  cont = []
  weight = []
  minimum = order.x[0]
  maximum = order.x[-1]
    
  #loop over the segments
  fig1 = pylab.figure(1)
  pylab.plot(order.x, order.y, 'r-')
  fig2 = pylab.figure(2)
  pylab.plot(model.x, model.y, 'r-')
  for i in range(contlist.low.size):
    if contlist.high[i] > minimum and contlist.low[i] < maximum:
      left = numpy.searchsorted(order.x, contlist.low[i])
      right = numpy.searchsorted(order.x, contlist.high[i])
      data = (numpy.mean(order.y[left:right]))
      pylab.figure(1)
      pylab.plot(order.x[left:right], order.y[left:right], 'b-')
        
      wave.append(order.x[(left+right)/2])
        
      left = numpy.searchsorted(model.x, contlist.low[i])
      right = numpy.searchsorted(model.x, contlist.high[i])
      telluric_model = (numpy.mean(model.y[left:right]))
      pylab.figure(2)
      pylab.plot(model.x[left:right], model.y[left:right], 'b-')
      
      cont.append(data/telluric_model)
      weight.append(float(right-left))
  
  pylab.show()
  
  #Interpolate/Extrapolate the continuum:
  wave = numpy.array(wave)
  cont = numpy.array(cont)
  weight = numpy.array(weight)
  weight = weight/numpy.max(weight)
  chisq = numpy.ones(5)
  #fit = extrap1d(interp1d(wave,cont,kind='linear'))
  """
  for i in range(chisq.size):
    fit = UnivariateSpline(wave,cont,w=weight,k=i+1)
    chisq[i] = numpy.sum((model.y - order.y/fit(order.x))**2)
    #chisq[i] = fit.get_residual()
    print "k = ",i+1,": ",chisq[i]
 
  #Use the best-fitting spline order
  fit = UnivariateSpline(wave,cont,w=weight,k=chisq.argmin()+1)  
  print "Using spline of order k=", chisq.argmin()+1, " to fit continuum"
  """
  pars = numpy.polyfit(wave,cont,3)
  pars, success = leastsq(ContFitFcn, pars, args=(wave,cont,weight))
  fit = numpy.poly1d(pars)
  
  order.cont = fit(order.x)
  
  pylab.plot(wave, cont, 'ro')
  pylab.plot(order.x, order.cont)
  pylab.plot(order.x, order.y)
  pylab.show()
  pylab.plot(order.x, order.y/order.cont)
  pylab.plot(model.x, model.y)
  pylab.show()
  
  #pylab.plot(order.x, order.y, label="data")
  #pylab.plot(wave, cont/fit(wave), 'ro')
  #pylab.plot(model.x, model.y, label="model")
  #pylab.legend(loc=3)
  #pylab.show()

  return order
  
#This function will fit the continuum in the regions given
def FitContinuum2(order,model,contlist):
  print "Fitting continuum"
  wave = []
  cont = []
  data = []
  weight = []
  minimum = order.x[0]
  maximum = order.x[-1]
  model_high = 0.0  
  data_high = 0.0
  
  Model = UnivariateSpline(model.x, model.y, s=0)

  #loop over the segments
  for i in range(contlist.low.size):
    if contlist.high[i] > minimum and contlist.low[i] < maximum:
      left = numpy.searchsorted(order.x, contlist.low[i])
      right = numpy.searchsorted(order.x, contlist.high[i])
      for i in range(left,right):
        #data.append(float(order.y[i]/Model(order.x[i])))
        data.append(order.y[i]/model.y[i])
        wave.append(order.x[i])
        weight.append(numpy.sqrt(model.y[i]))
 
  #Interpolate/Extrapolate the continuum:
  wave = numpy.array(wave)
  #cont = numpy.array(cont)/model_high
  data = numpy.array(data)
  weight = numpy.array(weight)
  weight = weight/numpy.sum(weight)
  #weight = numpy.ones(data.size)
  
  pars = numpy.polyfit(wave,data,4)
  pars, success = leastsq(ContFitFcn, pars, args=(wave,data,weight))
  fit = numpy.poly1d(pars)
  order.cont = fit(order.x)
  return order
 

def ContFitFcn(pars, x, y, w):
  #retval = pars[0]
  retval = 0.0
  for i in range(len(pars)):
    retval = retval + pars[i]*x**float(len(pars) - 1 - i)
  return (retval - y)*w
  
def FitContinuum3(order,model,fitorder=2):
  print "Fitting continuum"
  wave = order.x.copy()
  flux = order.y.copy()
  model2 = model.y.copy()
  weight = (model2/model2.max())**10
  done = False
  while not done:
    done = True
    wave_mean = numpy.mean(wave)
    pars = numpy.polyfit(wave - wave_mean, flux/model2, fitorder)
    pars, success = leastsq(ContFitFcn, pars, args=(wave - wave_mean, flux/model2, weight))
    fit = numpy.poly1d(pars)
    residuals = flux/(model2*fit(wave - wave_mean))
    std = numpy.std(residuals)
    mean = numpy.mean(residuals)
    badindices = numpy.where(numpy.abs(residuals-mean) > std*3)[0]
    wave = numpy.delete(wave, badindices)
    flux = numpy.delete(flux, badindices)
    model2 = numpy.delete(model2, badindices)
    weight = numpy.delete(weight, badindices)
    if badindices.size > 0:
      done = False
  order.cont = fit(order.x - wave_mean)
  return order

 
def FitContinuum(order, model, condition=0.95, tol=3, fitorder=5):
  print "Fitting continuum"
  #Fit continuum, using all points with model transmission > condition
  done = False
  while not done:
    done = True
    wave = model.x[model.y > condition]
    cont = order.y[model.y > condition]/model.y[model.y > condition]
  
    #make sure there are some points on the edges of the order
    if wave[0] > order.x[0] + 1.0 or wave[-1] < order.x[-1] - 1.0:
      condition = condition - 0.005
      done = False
  
  done = False
  while not done:
    done = True
    fit = numpy.poly1d(numpy.polyfit(wave - wave.mean(), cont, fitorder))
    resid = cont - fit(wave - wave.mean())
    mean = numpy.mean(resid)
    std = numpy.std(resid)
    badvals = numpy.abs(resid - mean) > std*tol
    if numpy.sum(badvals) > 0:
      done = False
    deleteindices = []
    for i in range(badvals.size):
      if badvals[i]:
        deleteindices.append(i)
        
    cont = numpy.delete(cont, deleteindices)
    wave = numpy.delete(wave, deleteindices)
      
  fit = numpy.poly1d(numpy.polyfit(wave - wave.mean(), cont, order))
  order.cont = fit(order.x - wave.mean())
  
  return order
  
  
#allows interp1d to extrapolate
def extrap1d(interpolator):
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)

    def ufunclike(xs):
        return numpy.array(map(pointwise, numpy.array(xs)))

    return ufunclike


for arg in sys.argv:
  if "debug" in arg:
    debug = True

if __name__=="__main__":
  #Parse command-line arguments for filename as well as optional abundances
  try:
    filename = sys.argv[1]
  except IndexError:
    filename = raw_input("Enter filename to fit: ") 
  outfile=open("chisq_summary.dat", "w")
  outfile.write("#CH4\tH20\tCO\tR\tX^2\n")
  outfile.close()
  j=2
  water = None
  ch4 = None
  resolution=None
  angle=None
  co=None
  o2 = None
  for i in range(j,len(sys.argv)):
    if j >=len(sys.argv):
      break
    elif "h2o" in sys.argv[j]:
      water = float(sys.argv[j+1])
      j = j+1
      print "water = ", water
    elif "ch4" in sys.argv[j]:
      ch4 = float(sys.argv[j+1])
      j = j+1
      print "CH4 = ", ch4
    elif "resolution" in sys.argv[j]:
      resolution = float(sys.argv[j+1])
      j=j+1
      print "Resolution = ", resolution
    elif "angle" in sys.argv[j]:
      angle = float(sys.argv[j+1])
      j = j+1
      print "angle = ", angle
    elif "co" in sys.argv[j]:
      co = float(sys.argv[j+1])
      j = j+1
      print "CO = ", co
    elif "o2" in sys.argv[j]:
      o2 = float(sys.argv[j+1])
      j = j+1
      print "O2 = ", o2
    j = j+1 
  Main(filename, humidity=water, resolution=resolution, angle=angle, ch4=ch4, co=co, o2=o2)
  
  
