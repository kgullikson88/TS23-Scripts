import FitsUtils
import FindContinuum
import sys
import os
import pylab
import pyfits
import numpy

goodregions = [[670.13, 670.73],
               [671.04, 671.7],
               [672.9, 673.3],
               [673.55, 674.3],
               [674.6, 674.9],
               [675.9, 676.7],
               [677.5, 679.5],
               [679.8, 680.4]]

if __name__ == "__main__":
  fileList = []
  for arg in sys.argv[1:]:
    fileList.append(arg)

  blazefile = [f for f in os.listdir("./") if "BLAZE" in f][0]
  blaze = FitsUtils.MakeXYpoints(blazefile, errors=2)

  emfiles = [f for f in os.listdir("./") if f.startswith("em") and f.endswith(".tbl")]

  iteration = 0
  outfile = open("SNR.dat", "w")
  for fname in fileList:
    number = int(fname[2:7])
    hdulist = pyfits.open(fname)
    try:
      if len(hdulist) == 1:
        hdulist.close()
        orders = FitsUtils.MakeXYpoints(fname, errors=2)
      else:
        hdulist.close()
        orders = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux")
    except TypeError:
      continue

    """
    print len(orders)
    for i, order in enumerate(orders):
      order.y /= blaze[i].y
      order.cont = FindContinuum.Continuum(order.x, order.y, lowreject=2, highreject=2)
      orders[i] = order.copy()
    """
      
    order = orders[7]
    order.y /= blaze[7].y
    order.cont = FindContinuum.Continuum(order.x, order.y, lowreject=2, highreject=2)

    for j, region in enumerate(goodregions):
      left = numpy.searchsorted(order.x, region[0])
      right = numpy.searchsorted(order.x, region[1])
      if j == 0:
        goodindices = numpy.arange(left,right)
      else:
        goodindices = numpy.r_[goodindices, numpy.arange(left, right)]0.5

    #Get S/N
    snr = 1.0 / numpy.std( (order.y/order.cont)[goodindices] )

    #Get Emeter counts
    emfile = "em%i.tbl" %number
    if emfile in emfiles:
      counts = numpy.loadtxt(emfile, usecols=(1,))
      emcounts = numpy.sum(counts)
    else:
      emcounts = 0
    
    print fname, emcounts, snr, '\n\n\n'
    outfile.write("%s\t%i\t%g\n" %(fname, emcounts, snr))
    
    pylab.plot(order.x, order.y/order.cont + iteration*0.2, 'k-')
    iteration += 1

  pylab.show()
  outfile.close()
