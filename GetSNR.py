import sys
import os
import pylab
from astropy.io import fits as pyfits
import numpy as np
import FittingUtilities
import HelperFunctions

goodregions = [[670.13, 670.73],
               [671.04, 671.7],
               [672.9, 673.3],
               [673.55, 674.3],
               [674.6, 674.9],
               [675.9, 676.7],
               [677.5, 679.5],
               [679.8, 680.4]]

def main1():
  fileList = []
  for arg in sys.argv[1:]:
    fileList.append(arg)

  blazefile = [f for f in os.listdir("./") if "BLAZE" in f][0]
  blaze = HelperFunctions.ReadFits(blazefile, errors=2)

  emfiles = [f for f in os.listdir("./") if f.startswith("em") and f.endswith(".tbl")]

  iteration = 0
  outfile = open("SNR.dat", "w")
  for fname in fileList:
    number = int(fname[2:7])
    hdulist = pyfits.open(fname)
    try:
      if len(hdulist) == 1:
        hdulist.close()
        orders = HelperFunctions.ReadFits(fname, errors=2)
      else:
        hdulist.close()
        orders = HelperFunctions.ReadFits(fname, extensions=True, x="wavelength", y="flux")
    except TypeError:
      continue

    """
    print len(orders)
    for i, order in enumerate(orders):
      order.y /= blaze[i].y
      order.cont = FittingUtilities.Continuum(order.x, order.y, lowreject=2, highreject=2)
      orders[i] = order.copy()
    """
      
    order = orders[7]
    order.y /= blaze[7].y
    order.cont = FittingUtilities.Continuum(order.x, order.y, lowreject=2, highreject=2)

    for j, region in enumerate(goodregions):
      left = np.searchsorted(order.x, region[0])
      right = np.searchsorted(order.x, region[1])
      if j == 0:
        goodindices = np.arange(left,right)
      else:
        goodindices = np.r_[goodindices, np.arange(left, right)]

    #Get S/N
    snr = 1.0 / np.std( (order.y/order.cont)[goodindices] )

    #Get Emeter counts
    emfile = "em%i.tbl" %number
    if emfile in emfiles:
      counts = np.loadtxt(emfile, usecols=(1,))
      emcounts = np.sum(counts)
    else:
      emcounts = 0
    
    print fname, emcounts, snr, '\n\n\n'
    outfile.write("%s\t%i\t%g\n" %(fname, emcounts, snr))
    
    pylab.plot(order.x, order.y/order.cont + iteration*0.2, 'k-')
    iteration += 1

  pylab.show()
  outfile.close()



def main2():
  fileList = []
  for arg in sys.argv[1:]:
    fileList.append(arg)

  for fname in fileList:
    orders = HelperFunctions.ReadFits(fname, extensions=True, x="wavelength", y="flux", cont="continuum")
    bestorder = 0
    bestsnr = -1
    for i, order in enumerate(orders):
      order.cont = FittingUtilities.Continuum(order.x, order.y, lowreject=2, highreject=3)
      snr = 1.0 / np.std( (order.y/order.cont) )
      #print "S/N in order %i = %g" %(i+1, snr)
      if snr > bestsnr:
        bestorder = i
        bestsnr = snr
    print "\n******  %s  ******" %fname
    print "******  Best order: %i ******" %(bestorder+1)
    print "******  with S/N = %g  ******" %bestsnr
    order = orders[bestorder]
    #pylab.plot(order.x, order.y/order.cont)
    #pylab.show()



if __name__ == "__main__":
  main2()
