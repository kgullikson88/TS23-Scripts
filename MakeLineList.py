import numpy as np
import matplotlib.pyplot as plt
import FitsUtils
from astropy.io import fits as pyfits
from astropy import units
import os
import scipy.signal
import DataStructures

bclength = 1000  #Boxcar smoothing length

def main1():
  fitsfile = "../17Vul/20120714/Star_17_Vul.fits"
  hdulist = pyfits.open(fitsfile)
  orders = FitsUtils.MakeXYpoints(hdulist[0].header, hdulist[0].data)[::-1]
  hdulist.close()

  boxcar = np.ones(bclength)/float(bclength)


  lines = []
  for order, index in zip(orders, range(len(orders))):
    smoothed = np.convolve(order.y, boxcar, mode='same')
    residuals = order.y - smoothed
    std = np.std(residuals)
    linepoints = np.where(np.logical_and(residuals[bclength:-bclength] - residuals.mean() < std, order.y[bclength:-bclength] > 0.9*np.max(order.y[bclength:-bclength])))[0] + bclength

    #Find all sets of consecutive points
    points = []
    for line in linepoints:
      if len(points) == 0 or int(line) - 1 == points[-1]:
        points.append(int(line))
      else:
        lines.append(order.x[int(np.median(points) + 0.5)])
        points = [int(line),]
        print lines[-1]
        yval = order.y[int(np.median(points) + 0.5)]
        plt.plot((lines[-1], lines[-1]), (yval-0.1, yval-0.2), 'r-')
    
    plt.figure()
    plt.title("Order number %i" %(len(orders) - index))
    plt.plot(order.x, order.y, 'k-')
    plt.plot(order.x, smoothed, 'r-')

    #plt.show()

  np.savetxt("Linelist.dat", lines, fmt="%.8f")


def main2():
  filename = os.environ["HOME"] + "/School/Research/Useful_Datafiles/Telluric_Visible.dat"
  print "Reading telluric model"
  x,trans = np.loadtxt(filename, unpack=True)
  x = x[::-1]
  trans = trans[::-1]

  boxcar = np.ones(bclength)/float(bclength)
  smoothed = np.convolve(trans, boxcar, mode='same')
  residuals = trans - smoothed
  std = np.std(residuals[bclength:-bclength])

  #plt.plot(x, residuals - residuals[bclength:-bclength].mean())
  #plt.plot(x, smoothed)
  #plt.plot(x, (-std)*np.ones(x.size))
  #plt.show()

  #linepoints = np.where(np.logical_and(residuals[bclength:-bclength] - residuals.mean() < std, trans[bclength:-bclength] > 0.9*np.max(trans[bclength:-bclength])))[0] + bclength
  linepoints = np.where(residuals[bclength:-bclength] - residuals[bclength:-bclength].mean() < -std)[0] + bclength
  linepoints = np.where(trans < 0.98)[0]

  print "Finding lines"
  points = []
  lines = []
  for line in linepoints:
    #print len(points)
    if len(points) == 0 or int(line) - 1 == points[-1]:
      points.append(int(line))
    else:
      index = int(np.median(points) + 0.5)
      if len(points) > 1:
        minindex = trans[points[0]:points[-1]].argmin() + points[0]
      else:
        minindex = points[0]
      if trans[minindex] < 0.95 and trans[minindex] > 0.1:
        lines.append(x[minindex])
        yval = trans[minindex]
      points = [int(line),]

  """
  #Make sure there are no points too close to each other
  tol = 0.05
  lines = sorted(lines)
  for i in range(len(lines) - 2, 0, -1):
    if np.abs(lines[i] - lines[i-1]) < tol:
      del lines[i]
      del lines[i-1]
    elif np.abs(lines[i] - lines[i+1]) < tol:
      del lines[i+1]
      del lines[i]
    else:
      index = np.searchsorted(x,lines[i]) - 1
      yval = trans[index]
      plt.plot((lines[i], lines[i]), (yval-0.05, yval-0.1), 'r-')
  """
  plt.plot(x, trans, 'k-')
  for line in lines:
    idx = np.searchsorted(x, line)
    plt.plot([x[idx], x[idx]], [trans[idx]-0.05, trans[idx]-0.1], 'r-')
  plt.show()
  np.savetxt("Linelist3.dat", lines, fmt="%.8f")


"""
  Function to find the spectral lines, given a model spectrum
  spectrum:        An xypoint instance with the model
  tol:             The line strength needed to count the line (0 is a strong line, 1 is weak)
  linespacing:     The minimum spacing between two consecutive lines
"""
def FindLines(spectrum, tol=0.99, linespacing = 0.01, debug=False):
  distance = 0.01
  xspacing = float(max(spectrum.x) - min(spectrum.x))/float(spectrum.size())
  N = int( linespacing / xspacing + 0.5)
  lines = list(scipy.signal.argrelmin(spectrum.y, order=N)[0])
  if debug:
    plt.plot(spectrum.x, spectrum.y)
  for i in range(len(lines)-1, -1, -1):
    idx = lines[i]
    xval = spectrum.x[idx]
    yval = spectrum.y[idx]
    if yval < tol:
      plt.plot([xval, xval], [yval-0.01, yval-0.03], 'r-')
    else:
      lines.pop(i)
  plt.show()
  return np.array(lines)
#np.savetxt("Linelist4.dat", lines, fmt="%.8f")



if __name__ == "__main__":
  filename = os.environ["HOME"] + "/School/Research/Useful_Datafiles/Telluric_Visible.dat"
  filename = "tell.dat"
  print "Reading telluric model"
  x,trans = np.loadtxt(filename, unpack=True)
  model = DataStructures.xypoint(x=x[::-1], y=trans[::-1])
  FindLines(model, debug=True)
