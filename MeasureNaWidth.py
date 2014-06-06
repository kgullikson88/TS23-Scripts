"""
This is a script to measure the equivalent width of interstellar sodium.
The main function is Measure(filename), which returns the equivalent
widths of both sodium D lines, as well as the reddening E(B-V).
The relationship from EW to E(B-V) is taken from 
Poznanski et al. 2012 MNRAS, 426, 1465
"""

import matplotlib.pyplot as plt
import numpy
from scipy.optimize import leastsq
import SpectralTypeRelations
import HelperFunctions
import pyspeckit
import sys
import os
from astropy.io import fits

ordernum = 30  #The order number with the Sodium lines in it
D1_wave = 589.592
D2_wave = 588.995

#Make matplotlib interactive
plt.ion()


def GetReddening(D1_ew, D2_ew):
  #Use the sum of the lines, since it has the lowest error
  EB_V = 10**(1.17*(D1_ew + D2_ew) - 1.85)
  return EB_V


def Measure(filename, figure, title):
  """
  Given a filename, this will plot the appropriate order
  and help the user measure the equivalent width of
  both sodium D lines
  """
  orders = HelperFunctions.ReadExtensionFits(filename)
  order = orders[ordernum]
  left = numpy.searchsorted(order.x, 588)
  right = numpy.searchsorted(order.x, 590.5)
  order = order[left:right]

  x = pyspeckit.units.SpectroscopicAxis(order.x, units='nm')
  spec = pyspeckit.Spectrum(data=order.y/order.cont, xarr=x, error=order.err/order.cont, units='nm')
  pl = spec.plotter(autorefresh=True, figure=figure)
  plt.title(title)

  # Fit the baseline in this region (continuum)
  spec.baseline(subtract=False, order=3, interactive=True)
  plt.show()
  done = raw_input("hit enter ")

  # Fit voigt profiles to the lines
  fitguesses=[-0.5, D1_wave, 0.08, 0.0,
              -0.5, D2_wave, 0.08, 0.0]
  spec.specfit(fittype='voigt', multifit=True, guesses=fitguesses)
  plt.draw()

  # Grab the fitted parameters and their errors
  fitpars = spec.specfit.modelpars
  fiterrs = spec.specfit.modelerrs

  # Sometimes, the lines 'switch' places, so check that
  if fitpars[1] > fitpars[5]:
    fitpars = fitpars[4:] + fitpars[:4]
    fiterrs = fiterrs[4:] + fiterrs[:4]

  #Determine the start and end of each line
  start1 = fitpars[1] - 7*numpy.sqrt(fitpars[2]**2 + fitpars[3]**2)
  end1 = fitpars[1] + 7*numpy.sqrt(fitpars[2]**2 + fitpars[3]**2)
  start2 = fitpars[5] - 7*numpy.sqrt(fitpars[6]**2 + fitpars[7]**2)
  end2 = fitpars[5] + 7*numpy.sqrt(fitpars[6]**2 + fitpars[7]**2)

  #Find the equivalent width for both lines
  ews = []
  for line in [[start1, end1], [start2, end2]]:
    xmin=numpy.searchsorted(order.x, line[0])
    xmax=numpy.searchsorted(order.x, line[1])
    ew = spec.specfit.EQW(xmin=xmin, xmax=xmax, fitted=True, continuum=1.0, plot=True)
    ews.append(ew*10)  #Save the equivalent width in angstroms
  done = raw_input("hit enter ")
  print "\n\n"

  return ews, GetReddening(*ews)


if __name__ == "__main__":
  #Read in the lines in the current log file
  logfilename = "%s/Dropbox/School/Research/AstarStuff/TargetLists/Na_Params.csv" %os.environ['HOME']
  logfile = open(logfilename, "r")
  lines = logfile.readlines()
  logfile.close()

  #Parse command line arguments
  fileList = []
  for arg in sys.argv[1:]:
    fileList.append(arg)

  #For each file, determine the EW and E(B-V)
  fig = plt.figure()
  for fname in fileList:
    #First, check to make sure that this target was not already measured
    measure = True
    overwrite = 9999999
    starname = fits.getheader(fname)["object"]
    print starname
    #Check if this star already has an entry in the logfile
    for linenum, line in enumerate(lines):
      star = line.split("|")[0].strip()
      if star.lower() == starname.lower():
        measure = False
        print "%s found in the master logfile." %starname
        inp = raw_input("Do you want to overwrite (y/n)? ")
        if inp.lower() == "y":
          measure = True
          overwrite = linenum

    if measure:
      ews, reddening = Measure(fname, fig, starname)
      print ews, reddening

      outline = "%s | %.4f | %.4f | %.4f\n" %(starname, ews[0], ews[1], reddening)
      if overwrite > len(lines):
        lines.append(outline)
      else:
        lines[overwrite] = outline

  #Finally, output the new logfile
  logfile = open(logfilename, "w")
  logfile.writelines(lines)
  logfile.close()