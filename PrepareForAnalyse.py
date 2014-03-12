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
import matplotlib.pyplot as plt
import numpy
import HelperFunctions
import FittingUtilities
import DataStructures

homedir = os.environ["HOME"]
rootdir = "%s/Applications/AnalyseBstar" %homedir
outdir = "%s/spectra" %rootdir


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

    # Re-fit continuum
    plt.figure(1)
    for i, order in enumerate(orders):
      order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=contorder, lowreject=1, highreject=10)

      plt.figure(1)
      plt.plot(order.x, order.y/order.cont, 'k-')
    plt.show()

    # Find overlap
    for i, order in enumerate(orders):
      if i < len(orders)-1:
        first = numpy.searchsorted(order.x, orders[i+1].x[0])
        if first < order.size() - 2:
          temp = FittingUtilities.RebinData(orders[i+1], order.x[first:])
          plt.plot(temp.x, (order.y[first:]/order.cont[first:]) / (temp.y/temp.cont))
      #plt.plot(order.x, order.cont, 'r--')
    plt.show()

    # Combine orders
    xgrid = numpy.arange(orders[0].x[0], orders[-1].x[-1], 0.01)
    output = DataStructures.xypoint(x=xgrid)
    #for order in orders:
      
