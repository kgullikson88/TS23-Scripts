import pyfits
import FitsUtils
import sys
from scipy.interpolate import InterpolatedUnivariateSpline as interp
import matplotlib.pyplot as plt
import DataStructures
import os
import FitTellurics_McDonald
import FindContinuum
import numpy


def ReadCorrectedFile(fname):
  orders = []
  headers = []
  hdulist = pyfits.open(fname)
  numorders = len(hdulist)
  for i in range(1, numorders):
    order = hdulist[i].data
    xypt = DataStructures.xypoint(x=order.field("wavelength"),
                                  y=order.field("model"),
                                  cont=order.field("continuum"),
                                  err=order.field("error"))

    orders.append(xypt)
    headers.append(hdulist[i].header)
  return orders, headers


def Correct(original, corrected, offset=None):
  #Read in the data and model
  original_orders = FitsUtils.MakeXYpoints(original, errors=2)
  corrected_orders, corrected_headers = ReadCorrectedFile(corrected)
  if offset == None:
    offset = len(original_orders) - len(corrected_orders)

  #Find and read in blaze function (MUST BE IN CURRENT DIRECTORY!)
  files = os.listdir("./")
  blazefile = [fname for fname in files if fname.startswith("BLAZE")][0]
  blaze_orders = FitsUtils.MakeXYpoints(blazefile, errors=2)
  blaze_functions = []
  for order in blaze_orders:
    blaze_functions.append( interp(order.x, order.y) )

  for i in range(len(original_orders)):
    original_orders[i].y /= blaze_functions[i](original_orders[i].x)

  for i in range(offset, len(original_orders)):
    data = original_orders[i]
    data.cont = FindContinuum.Continuum(data.x, data.y)
    try:
      model = corrected_orders[i-offset]
      header = corrected_headers[i-offset]
      print "Order = %i\nHumidity: %g\nO2 concentration: %g\n" %(i, header['h2oval'], header['o2val'])
    except IndexError:
      model = DataStructures.xypoint(x=data.x, y=numpy.ones(data.x.size))
      print "Warning!!! Telluric Model not found for order %i" %i
    data.y /= model.y
    original_orders[i] = data.copy()
  return original_orders

if __name__ == "__main__":
  for number in sys.argv[1:]:
    #Number is for for e.g. KG11186.fits
    original = "KG%s.fits" %number
    corrected = "Corrected_%s-0.fits" %number
    #original = "HIP_16147.fits"
    #corrected = "HIP_16147-0.fits"
    corrected = "Corrected_11319-0.fits"
    print original
    print corrected
    outfilename = "%s_telluric_corrected.fits" %(original.split(".fits")[0])
    print "Outputting to %s" %outfilename

    try:
      corrected_orders = Correct(original, corrected, offset=2)
    except IOError:
      continue
    for i, data in enumerate(corrected_orders):
      #Set up data structures for OutputFitsFile
      columns = {"wavelength": data.x,
                 "flux": data.y,
                 "continuum": data.cont,
                 "error": data.err}
      if i == 0:
        FitTellurics_McDonald.OutputFitsFile(columns, original, outfilename, mode="new")
      else:
        FitTellurics_McDonald.OutputFitsFile(columns, outfilename, outfilename)
