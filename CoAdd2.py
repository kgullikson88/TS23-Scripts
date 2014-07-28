import FitsUtils
import DataStructures
import FindContinuum
import sys
from scipy.interpolate import InterpolatedUnivariateSpline as interp
import numpy as np
import pylab
from astropy.io import fits as pyfits



if __name__ == "__main__":
  fileList = []
  for arg in sys.argv[1:]:
    fileList.append(arg)

  all_data = []
  numorders = []
  header = pyfits.getheader(fileList[0])
  name = header["OBJECT"]
  for fname in fileList:
    header = pyfits.getheader(fname)
    if header["OBJECT"] != name:
      sys.exit("Error! Trying to CoAdd observations of different stars!")
    observation = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", cont="continuum", errors="error")
    all_data.append(observation)
    numorders.append(len(observation))

  if any(n != numorders[0] for n in numorders):
    print "Error! Some of the files had different numbers of orders!"
    for i in range(len(fileList)):
      print fileList[i], numorders[i]
    sys.exit()

  #If we get this far, all is well. Add each order indidually
  numorders = numorders[0]
  outfilename = "%s.fits" %(name.replace(" ", "_"))
  print "outputting to %s" %outfilename
  column_list = []
  for i in range(numorders):
    total = all_data[0][i].copy()
    total.err = total.err**2
    for observation in all_data[1:]:
      flux = interp(observation[i].x, observation[i].y)
      error = interp(observation[i].x, observation[i].err**2, k=1)
      total.y += flux(total.x)
      total.err += error(total.x)
    total.err = np.sqrt(total.err)
    total.cont = FindContinuum.Continuum(total.x, total.y, fitorder=3, lowreject=2, highreject=5)

     #Set up data structures for OutputFitsFile
    columns = {"wavelength": total.x,
               "flux": total.y,
               "continuum": total.cont,
               "error": total.err}
    column_list.append(columns)
    pylab.plot(total.x, total.y)
    pylab.plot(total.x, total.cont)

  pylab.show()
  FitsUtils.OutputFitsFileExtensions(column_list, fileList[0], outfilename, mode="new")
    
