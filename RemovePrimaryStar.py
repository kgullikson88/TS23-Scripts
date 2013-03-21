#!/opt/local/bin/python
import numpy
import os
import sys
import pylab
from scipy.interpolate import UnivariateSpline
import SpectralTypeRelations
from collections import defaultdict
from PlotBlackbodies import Planck
import Units
import DataStructures
import Correlate
import MakeModel
import RotBroad
import FindContinuum
import time
import FitTellurics_McDonald as Outputter

"""
   This function removes a primary model from data  
"""



def GetModel(data, model, vel=0.0, vsini=15*Units.cm/Units.km):
  #Broaden
  model = RotBroad.Broaden(model, vsini=vsini, findcont=True)
  
  #Interpolate Model, and note the endpoints
  first = model.x[0]*(1.+vel/Units.c)
  last = model.x[-1]*(1.+vel/Units.c)
  a = []
  for i in range(1, model.x.size):
    a.append(model.x[i] - model.x[i-1])
  model_fcn = UnivariateSpline(model.x, model.y, s=0)

  data2 = []
  for order in data:
    data2.append(order.copy())

  ##############################################################
  #Begin main loop over the orders
  ##############################################################
  output_orders = []
  for i in range(len(data2)):
    order = data2[i]
    order.cont = FindContinuum.Continuum(order.x, order.y, lowreject=2, highreject=5)

    model2 = DataStructures.xypoint(x=order.x, y=model_fcn(order.x*(1.+vel/Units.c)))
    
    #Get model continuum in this section
    model2.cont = FindContinuum.Continuum(model2.x, model2.y, lowreject=1, fitorder=3)
    
    #Reduce resolution
    model2 = MakeModel.ReduceResolution(model2.copy(), 60000)

    output_orders.append(model2.copy())
    
  return output_orders




if __name__ == "__main__":
  import FitsUtils
  import os
  import sys
  home = os.environ["HOME"]
  try:
    datafile = sys.argv[1]
    print "Removing primary star from  ", datafile
  except IndexError:
    print "Error! Must give .fits file!"
    sys.exit()
  
  modeldir = os.environ["HOME"] + "/School/Research/Models/Sorted/Stellar/Vband/"
  files = os.listdir(modeldir)
  modelfiles = defaultdict(list)
  for fname in files:
    try:
      temperature = float(fname.split("lte")[-1].split("-")[0])*100
    except ValueError:
      print "Skipping file %s" %fname
    modelfiles[temperature].append(modeldir+fname)

  #Read in data
  orders_original = tuple(FitsUtils.MakeXYpoints(datafile, extensions=True, x="wavelength", y="flux", errors="error"))
  orders_original = orders_original[::-1]

  #Check for command line arguments
  p_temp = 5800
  vsini = 15.0
  vel = 0.0
  if len(sys.argv) > 2:
    for arg in sys.argv[2:]:
      if "primary" in arg:
        p_temp = float(arg.split("=")[-1])
      elif "vsini" in arg:
        vsini = float(arg.split("=")[-1])
      elif "vel" in arg:
        vel = float(arg.split("=")[-1])
  

  #Get the best logg for a main sequence star with the given temperature
  MS = SpectralTypeRelations.MainSequence()
  p_spt = MS.GetSpectralType(MS.Temperature, p_temp)
  p_mass = MS.Interpolate(MS.Mass, p_spt)
  radius = MS.Interpolate(MS.Radius, p_spt)
  logg = numpy.log10(Units.G*p_mass*Units.Msun/(radius*Units.Rsun)**2)
  best_key = modelfiles.keys()[0]
  
  for key in modelfiles.keys():
    if numpy.abs(p_temp - key) < numpy.abs(p_temp - best_key):
      best_key = key
      logg_values = [float(fname.split("lte")[-1].split("-")[1]) for fname in modelfiles[key]]
 
  logg_index = numpy.argmin(numpy.array(logg_values - logg))
  modelfile = modelfiles[best_key][logg_index]

  #Read in model
  print "Model file: %s" %modelfile
  x,y = numpy.loadtxt(modelfile, usecols=(0,1), unpack=True)
  x = x*Units.nm/Units.angstrom
  y = 10**y
  model = DataStructures.xypoint(x=x, y=y)

  orders = GetModel(list(orders_original), model, vel=vel*Units.cm/Units.km, vsini=vsini*Units.cm/Units.km)
  
  outfilename = "%s-0.fits" %(datafile.split(".fits")[0])
  print "Outputting to %s" %outfilename

  for i, model in enumerate(orders[2:-1]):
    original = orders_original[i+2]
    original.cont = FindContinuum.Continuum(original.x, original.y, lowreject=2, highreject=2)
    pylab.figure(1)
    pylab.plot(original.x, original.y/original.cont, 'k-')
    pylab.plot(model.x, model.y/model.cont, 'r-')
    pylab.figure(2)
    pylab.plot(original.x, original.y/(original.cont * model.y/model.cont), 'k-')
    original.y /= model.y/model.cont
    
    columns = {"wavelength": original.x,
               "flux": original.y,
               "continuum": original.cont,
               "error": original.err}
    
    #Output
    if i == 0:
      Outputter.OutputFitsFile(columns, datafile, outfilename, mode="new")
    else:
      Outputter.OutputFitsFile(columns, outfilename, outfilename)
  pylab.show()

