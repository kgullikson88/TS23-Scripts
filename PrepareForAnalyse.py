"""
  This code takes a fits file that is telluric-corrected but not smoothed,
and prepares it for input to AnalyseBstar. It does the following:
==========================================================================

1: Fits normalization for each order
2: Combines orders
3: Outputs the full spectrum as an ascii file in the appropriate directory

"""

import sys
import os
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy
import HelperFunctions
import FittingUtilities
import DataStructures

homedir = os.environ["HOME"]
rootdir = "%s/Applications/AnalyseBstar" %homedir
outdir = "%s/spectra" %rootdir


def CombineOrders(orders):
  if len(orders) == 1:
    return orders[0]

  # Combine orders
  xgrid = numpy.arange(orders[0].x[0], orders[-1].x[-1], 0.01)
  output = DataStructures.xypoint(x=xgrid)
  output.y = numpy.ones(output.size())
  norm = numpy.zeros(output.size())
  for order in orders:
    fcn = spline(order.x, order.y, k=1)
    left = numpy.searchsorted(xgrid, order.x[0])
    right = numpy.searchsorted(xgrid, order.x[-1])
    output.y[left:right] += fcn(output.x[left:right])
    norm[left:right] += 1.0
  output.y /= norm
  return output



if __name__ == "__main__":
  fileList = []
  contorder = 1
  
  for arg in sys.argv[1:]:
    if '-c' in arg:
      contorder = int(arg.split("=")[-1])
    else:
      fileList.append(arg)

  for fname in fileList:
    orders = HelperFunctions.ReadFits(fname, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")

    for order in orders:
      order = order[20:-20]

    # Shift overlapping orders up so that the edges line up well
    for i, order in enumerate(orders[:-1]):
      if order.x[-1] > orders[i+1].x[0]:
        left = numpy.searchsorted(order.x, orders[i+1].x[0])
        fcn = spline(orders[i+1].x, orders[i+1].y, k=1)
        factor = order.y[left:]/fcn(order.x[left:])
        poly = FittingUtilities.Continuum(order.x[left:], factor, fitorder=1, lowreject=2, highreject=2)
        order.y[:left] /= poly[0]
        for j in range(i):
          orders[j].y /= poly[0]
        order.y[left:] /= poly

    # Re-fit continuum
    plt.figure(1)
    for i, order in enumerate(orders):
      # Use neighboring orders if possible
      if i > 0 and orders[i-1].x[-1] > order.x[0]:
        if i < len(orders)-1 and orders[i+1].x[0] < order.x[-1]:
          combine=(orders[i-1], order, orders[i+1])
        else:
          combine = (orders[i-1], order)
      elif i < len(orders)-1 and orders[i+1].x[0] < order.x[-1]:
        combine = (order, orders[i+1])
      else:
        combine = (order,)
      combination = CombineOrders(combine)
      combination.cont = FittingUtilities.Continuum(combination.x, combination.y, fitorder=contorder, lowreject=2, highreject=5)
      contfcn = spline(combination.x, combination.cont, k=1)
      order.cont = contfcn(order.x)

      
          

    # Combine orders
    xgrid = numpy.arange(orders[0].x[0], orders[-1].x[-1], 0.01)
    output = DataStructures.xypoint(x=xgrid)
    output.y = numpy.ones(output.size())
    norm = numpy.zeros(output.size())
    for order in orders:
      fcn = spline(order.x, order.y, k=1)
      left = numpy.searchsorted(xgrid, order.x[0])
      right = numpy.searchsorted(xgrid, order.x[-1])
      output.y[left:right] += fcn(output.x[left:right])
      #cont = FittingUtilities.Continuum(output.x[max(left-1000, 0):min(right+1000
      cfcn = spline(order.x, order.cont, k=1)
      output.cont[left:right] += cfcn(output.x[left:right])
      norm[left:right] += 1.0
    output.y /= output.cont
    plt.plot(output.x, output.y)
    plt.show()

    # Get instrument name from the header
    header = fits.getheader(fname)
    observatory = header["OBSERVAT"]
    if "ctio" in observatory.lower():
      instrument = "CHIRON"
      star = header["OBJECT"].replace(" ", "")
    else:
      instrument = header["INSTRUME"]
      if "ts23" in instrument.lower():
        instrument = "TS23"
        star = header["OBJECT"].replace(" ", "")
      elif "hrs" in instrument.lower():
        instrument = "HRS"
        star = header["OBJECT"].split()[0].replace("_", "")
      else:
        raise ValueError ("Unknown instrument: %s" %instrument)

    outfilename = "%s/%s/%s/%s.txt" %(outdir, instrument, star, star)
    print outfilename
    HelperFunctions.ensure_dir(outfilename)
    numpy.savetxt(outfilename, numpy.transpose((output.x*10.0, output.y)))

    #for i, order in enumerate(orders):
    #  outfilename = "%s/%s/%s/order%i.txt" %(outdir, instrument, star, i+1)
    #  numpy.savetxt(outfilename, numpy.transpose((order.x*10.0, order.y/order.cont)))
    
    
      
