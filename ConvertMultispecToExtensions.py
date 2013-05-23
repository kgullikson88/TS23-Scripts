import FitsUtils
import FindContinuum
import pyfits
import sys
import os
import numpy
import pylab


if __name__ == "__main__":
  fileList = []
  blazecorrect = False
  for arg in sys.argv[1:]:
    if "blaze" in arg:
      blazecorrect = True
    else:
      fileList.append(arg)
  for fname in fileList:
    outfilename = "%s-0.fits" %(fname.split(".fits")[0])
    header = pyfits.getheader(fname)

    try:
      orders = FitsUtils.MakeXYpoints(fname)
    except ValueError:
      orders = FitsUtils.MakeXYpoints(fname, errors=2)
    orders = orders[::-1]    #Reverse order so the bluest order is first
    if blazecorrect:
      header = pyfits.getheader(fname)
      try:
        blazefile = "%s.fits" %header['BLAZE']
      except KeyError:
        allfiles = os.listdir("./")
        blazefile = [f for f in allfiles if "BLAZE" in f][0]
      try:
        blaze = FitsUtils.MakeXYpoints(blazefile)
        blaze = blaze[::-1]
      except ValueError:
        blaze = FitsUtils.MakeXYpoints(blazefile, errors=2)
        blaze = blaze[::-1]
      except IOError:
        print "Error! blaze file %s does not exist!" %blazefile
        print "Not converting file %s" %fname
        continue
    for i, order in enumerate(orders):
      #Don't use the stuff near the picket fence
      if i >= 14 and i <=17:
        continue
      

      #Blaze correction
      if blazecorrect:
        order.y /= blaze[i].y
        order.err /= blaze[i].y

      order.cont = FindContinuum.Continuum(order.x, order.y, fitorder=3, lowreject=2, highreject=4)
      columns = columns = {"wavelength": order.x,
                           "flux": order.y,
                           "continuum": order.cont,
                           "error": order.err}

      
      if i == 0:
        FitsUtils.OutputFitsFileExtensions(columns, fname, outfilename, mode="new")
      else:
        FitsUtils.OutputFitsFileExtensions(columns, outfilename, outfilename, mode="append")
      
      
