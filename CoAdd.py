import numpy
from astropy.io import fits as pyfits
import FitsUtils
import sys
import os
from scipy.interpolate import InterpolatedUnivariateSpline as interp
import FitTellurics_McDonald
import matplotlib.pyplot as plt
import FindContinuum
import DataStructures


def ReadCorrectedFile(fname):
  orders = []
  headers = []
  hdulist = pyfits.open(fname)
  numorders = len(hdulist)
  for i in range(1, numorders):
    order = hdulist[i].data
    xypt = DataStructures.xypoint(x=order.field("wavelength"),
                                  y=order.field("flux"),
                                  cont=order.field("continuum"),
                                  err = order.field("error"))

    orders.append(xypt)
    headers.append(hdulist[i].header)
  return orders, headers


def CoAdd(files):
  orders, headers = ReadCorrectedFile(files[0])
  orders = FitsUtils.MakeXYpoints(files[0], errors=2)
  for i, order in enumerate(orders):
    orders[i].err = order.err**2

  fnum = 1
  for fname in files[1:]:
    orders2, headers2 = ReadCorrectedFile(fname)
    orders2 = FitsUtils.MakeXYpoints(fname, errors=2)
    for i, order in enumerate(orders2):
      fcn = interp(order.x, order.y)
      errfcn = interp(order.x, order.err**2)
      orders[i].y += fcn(orders[i].x)
      orders[i].err += errfcn(orders[i].x)
    fnum += 1

  for i, order in enumerate(orders):
    orders[i].err = numpy.sqrt(order.err)
  return orders

if __name__ == "__main__":
  fileList = sys.argv[1:]
  fileList = []
  for arg in sys.argv[1:]:
    if "name" in arg:
      starname = arg.split("=")[-1]
      allfiles = [fname for fname in os.listdir("./") if fname.endswith(".fits") and fname.startswith("KG")]
      for fname in allfiles:
        header = pyfits.getheader(fname)
        object_name = header["OBJECT"]
        if object_name == starname:
          print "Adding file %s to list" %fname
          fileList.append(fname)
    else:
      fileList.append(arg)
  orders = CoAdd(fileList)
  print "There are %i orders" %len(orders)

  #Find and read in blaze function (MUST BE IN CURRENT DIRECTORY!)
  files = os.listdir("./")
  blazefile = [fname for fname in files if fname.startswith("BLAZE")][0]
  blaze_orders = FitsUtils.MakeXYpoints(blazefile, errors=2)
  blaze_functions = []
  for order in blaze_orders:
    blaze_functions.append( interp(order.x, order.y) )

  for i, order in enumerate(orders):
    orders[i].y /= blaze_functions[i](order.x)
    orders[i].cont = FindContinuum.Continuum(orders[i].x, orders[i].y, lowreject=2, highreject=5)
    plt.plot(orders[i].x, orders[i].y/orders[i].cont)
  plt.show()

  #Get the name of the star 
  header = pyfits.getheader(fileList[0])
  object_name = header["OBJECT"]
  outfilename = "%s.fits" %(object_name.replace(" ", "_") )
  print "Outputting co-added data to %s" %outfilename

  for i, data in enumerate(orders):
    print "Outputting order %i of %i" %(i, len(orders))
    #Set up data structures for OutputFitsFile
    columns = {"wavelength": data.x,
               "flux": data.y,
               "continuum": data.cont,
               "error": data.err}
    if i == 0:
      FitTellurics_McDonald.OutputFitsFile(columns, fileList[0], outfilename, mode='new')
    else:
      FitTellurics_McDonald.OutputFitsFile(columns, outfilename, outfilename)
