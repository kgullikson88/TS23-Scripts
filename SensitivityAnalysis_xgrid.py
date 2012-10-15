print "Very beginning of script!"
import numpy
import os
import sys
from scipy.interpolate import UnivariateSpline
import SpectralTypeRelations
from collections import defaultdict
import Units
import DataStructures
import Correlate
#import MakeModel
from subprocess import check_call
import FitsUtils

"""
   This function performs a sensitivity analysis on reduced spectra
   This version is specifically tailored to use on the x-grid
   You MUST give the name of the input .fits file as the FIRST command line argument
   After that, you can give optional command line arguments that specify ranges of spectral types
      for both the primary and secondary stars.
   Example: To do a sensitivity analysis on file 'foo.fits', for primary spectral types from B4-A5, and
            secondary spectral types from G0-K5, you type the following:
            python SensitivityAnalysis.py foo.fits -primary=B4-A5 -secondary=G0-K5


def Planck(x,T):
  h = 6.626e-27
  c = 3e10
  k = 1.38e-16
  pi = numpy.pi
  return 2*pi*h*c**2/x**5*1.0/(numpy.exp(h*c/(x*k*T)) - 1.0)


print "Very beginning of script!"
#homedir = os.environ['HOME'] + "/"
outfiledir = os.getcwd() + "/Sensitivity/"

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)
        
ensure_dir(outfiledir)
print "Sensitivity directory created"

#Define the sections of each order to use (those without telluric contamination)
good_sections = {1: [[-1,-1],],
                 2: [[-1,-1],],
		 3: [[-1,-1],],
		 4: [[-1,-1],],
		 5: [[-1,-1],],
		 6: [[-1,-1],],
		 7: [[-1,-1],],
		 8: [[-1,-1],],
		 9: [[-1,-1],],
		 10: [[-1,-1],],
		 11: [[-1,-1],],
		 12: [[-1,-1],],
		 13: [[-1,-1],],
		 14: [[-1,-1],],
		 15: [[-1,-1],],
		 16: [[-1,-1],],
		 17: [[-1,-1],],
		 18: [[-1,-1],],
		 19: [[-1,686.7],],
		 20: [[-1, 1e9],],
		 21: [[661.4, 667.3],],
		 22: [[-1,-1],],
		 23: [[-1, 1e9]],
		 24: [[-1, 627.6],],
		 25: [[-1, 1e9],],
		 26: [[-1, 1e9],],
		 27: [[-1,-1],],
		 28: [[-1, -1],],
		 29: [[-1, 1e9],],
		 30: [[-1, 567.6],],
		 31: [[-1, 1e9],],
		 32: [[-1, 1e9],],
		 33: [[-1, 1e9],],
		 34: [[-1, 1e9],],
		 35: [[-1, 1e9],],
		 36: [[-1, 1e9],],
		 37: [[-1, 1e9],],
		 38: [[-1, 501.1],],
		 39: [[-1, 491.6], [492.8, 1e9]],
		 40: [[-1, 1e9],],
		 41: [[-1, 1e9],],
		 42: [[-1, -1],],
		 43: [[-1, -1],],
		 44: [[-1, -1],],
		 45: [[-1, 1e9],],
		 46: [[-1, 1e9],],
		 47: [[-1, 1e9],],
		 48: [[-1, 1e9],],
		 49: [[-1, 1e9],],
		 50: [[-1, 1e9],],
		 51: [[-1, 1e9],],
                 52: [[-1,-1],],
		 43: [[-1,-1],] }

"""

"""
def Add(data, model, prim_spt, sec_spt, age="MS", vel=0.0):
    
    This function will add a model to the data. The flux ratio is determined
    from the primary spectral type (prim_spt), the secondary spectral type
    (sec_spt), and the age of the system. The age keyword argument can be either
    a string, in which case the main sequence spectral-type - radius relations are used,
    or a number in which case the radii of the two stars are determined from 
    evolutionary tracks. The 'vel' keyword gives a radial velocity at which
    the model should be added (MUST be given in cm/s)
    
  if type(age) == str:
    #Main sequence stars!
    MS = SpectralTypeRelations.MainSequence()
    prim_radius = MS.Interpolate(MS.Radius, prim_spt)
    prim_temp = MS.Interpolate(MS.Temperature, prim_spt)
    sec_radius = MS.Interpolate(MS.Radius, sec_spt)
    sec_temp = MS.Interpolate(MS.Temperature, sec_spt)
  else:
    sys.exit("Sorry! Can only handle Main Sequence stars right now!")

  
  #Interpolate Model, and note the endpoints
  first = model.x[0]*(1.+vel/Units.c)
  last = model.x[-1]*(1.+vel/Units.c)
  model_fcn = UnivariateSpline(model.x, model.y, s=0)

  data2 = []
  for order in data:
    data2.append(order.copy())


  ##############################################################
  #Begin main loop over the orders
  for i in range(len(data2)):
    order = data2[i]
    prim_flux = Planck(order.x*Units.cm/Units.nm, prim_temp)*prim_radius**2
    sec_flux = Planck(order.x*Units.cm/Units.nm, sec_temp)*sec_radius**2
    fluxratio = sec_flux/prim_flux
    print "order %i flux ratio = %.3g" %(i+1, numpy.mean(fluxratio))

    model2 = DataStructures.xypoint(x=order.x, y=model_fcn(order.x*(1.+vel/Units.c)))
    
    #Reduce resolution
    model2 = Correlate.ReduceResolution(model2.copy(), 60000)
    
    #Get model continuum in this section
    weights = model2.y**2
    done = False
    x2 = model2.x.copy()
    y2 = model2.y.copy()
    fit_order = 2
    while not done:
      done = True
      fit = numpy.poly1d(numpy.polyfit(x2 - x2.mean(), y2, fit_order))
      residuals = y2 - fit(x2 - x2.mean())
      mean = numpy.mean(residuals)
      std = numpy.std(residuals)
      badpoints = numpy.where((residuals - mean) < -std)[0]
      if badpoints.size > 0 and x2.size - badpoints.size > 5*fit_order:
        done = False
        x2 = numpy.delete(x2, badpoints)
        y2 = numpy.delete(y2, badpoints)
    model2.cont = fit(model2.x - x2.mean())

    #Scale the model by the above scale factor and normalize
    scaled_model = (model2.y/model2.cont)*fluxratio

    #Add to the data
    data2[i].y = (scaled_model)*order.cont + order.y


  return data2

"""

if __name__ == "__main__":
  print os.getcwd()
  """
  try:
    datafile = sys.argv[1]
    print "Using ", datafile, "as template"
  except IndexError:
    print "Error! Must give .fits file!"
    sys.exit()
  
  tolerance = 5    #allow highest cross-correlation peak to be anywhere within 5 km/s of the correct velocity
  MS = SpectralTypeRelations.MainSequence()  #Class for interpolating stellar info from the spectral type

  logfile = open(outfiledir+"summary.dat", "w")
  logfile.write("Parent SpT\tSecondary SpT\tParent Mass\tSecondary Mass\tMass Ratio\tPercent Detected\tAverage Significance\n")

  #Make some lists that we will loop through in the analysis
  parent_spts = ["B0", "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9",
                 "A0", "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9"]
  sec_spts = ["F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9",
              "G0", "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9",
              "K0", "K1", "K2", "K3", "K4", "K5", "K6", "K7", "K8", "K9",
              "M0", "M1", "M2", "M3", "M4", "M5"]
  velocitylist = [-400,-440,-360,-300,-250,-210,-140,-90-30,0,50,110,140,200,260,310,350,390]

  #Get the models. Must unpack from archive first
  print "Unpacking model archive"
  check_call("tar -xzf models.tar.gz", shell=True)
  modeldir = "./models/"
  files = os.listdir(modeldir)
  modelfiles = defaultdict(list)
  for fname in files:
    temperature = float(fname.split("lte")[-1][:2])*100
    modelfiles[temperature].append(modeldir+fname)

  #Read in data
  orders_original = tuple(FitsUtils.MakeXYpoints(datafile))

  #Check for command line arguments specifying spectral type endpoints
  p_left = 0
  p_right = len(parent_spts)
  s_left = 0
  s_right = len(sec_spts)
  print len(sys.argv)
  if len(sys.argv) > 2:
    for arg in sys.argv[2:]:
      if "primary" in arg:
	first = arg.split("=")[-1].split("-")[0]
	last = arg.split("=")[-1].split("-")[1]
	for i in range(len(parent_spts)):
	  spt = parent_spts[i]
	  if spt == first:
	    p_left = i
	  if spt == last:
	    p_right = i+1
	#if p_right < len(p_spt)-1:
	#  p_right += 1
      elif "secondary" in arg:
        first = arg.split("=")[-1].split("-")[0]
	last = arg.split("=")[-1].split("-")[1]
	for i in range(len(sec_spts)):
	  spt = sec_spts[i]
	  if spt == first:
	    s_left = i
	  if spt == last:
	    s_right = i+1
	#if s_right < len(s_spt)-1:
	#  s_right += 1

  if p_left > p_right:
    temp = p_left
    p_left = p_right-1
    p_right = temp+1
  if s_left > s_right:
    temp = s_left
    s_left = s_right-1
    s_right = temp+1

  ############################################
  #Start looping!
  ############################################
  print "Beginning main loop now!"
  for p_spt in parent_spts[p_left:p_right]:
    for s_spt in sec_spts[s_left:s_right]:
      found = 0.0
      sig = []
      for velocity in velocitylist:
        #Figure out what the best model is from the secondary spectral type
        s_mass = MS.Interpolate(MS.Mass, s_spt)
        p_mass = MS.Interpolate(MS.Mass, p_spt)
        radius = MS.Interpolate(MS.Radius, s_spt)
        temperature = MS.Interpolate(MS.Temperature, s_spt)
        logg = numpy.log10(Units.G*s_mass*Units.Msun/(radius*Units.Rsun)**2)
        best_key = modelfiles.keys()[0]
        for key in modelfiles.keys():
          if numpy.abs(temperature - key) < numpy.abs(temperature - best_key):
            best_key = key
        best_logg = float(modelfiles[best_key][0].split("lte")[-1][3:7])
        modelfile = modelfiles[best_key][0]
        for fname in modelfiles[best_key]:
          model_logg = float(fname.split("lte")[-1][3:7])
          if numpy.abs(logg - model_logg) < numpy.abs(logg - best_logg):
            best_logg = logg
            modelfile = fname
        
        x,y = numpy.loadtxt(modelfile, usecols=(0,1), unpack=True)
        x = x*Units.nm/Units.angstrom
        y = 10**y
        model = DataStructures.xypoint(x=x, y=y)

        orders = Add(list(orders_original), model, p_spt, s_spt, vel=velocity*Units.cm/Units.km)
        outfilebase = outfiledir+p_spt+"_"+s_spt+"_v%i" %velocity
        FitsUtils.OutputFitsFile(datafile, orders, outfilename=outfilebase+".fits")

        #Cross-correlate with original model
        Correlate.PyCorr(orders, models=[[x,y],], segments=good_sections, outfilename=outfilebase+"_CC.dat")

        vel, corr = numpy.loadtxt(outfilebase+"_CC.dat", unpack=True)
        maxindex = corr.argmax()
        std = corr.std()
        mean = corr.mean()
        corr = (corr - mean)/std
        maxvel = vel[maxindex]
        significance = corr[maxindex]
        if maxvel-tolerance <= velocity and maxvel+tolerance >= velocity:
          sig.append(significance)
          found += 1.

      print "Found %i signals with a mean significance of %.3g" %(found, numpy.mean(sig))
      outfilestring = p_spt+"\t\t"+s_spt+"\t\t%.3f\t\t%.3f\t\t%.3f\t\t%.5f\t\t%.5f\n" %(p_mass, s_mass, s_mass/p_mass, found*100./float(len(velocitylist)), numpy.mean(sig))
      print outfilestring
      logfile.write(outfilestring)
  

  logfile.close()
  """
