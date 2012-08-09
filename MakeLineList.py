import numpy
import matplotlib.pyplot as plt
import FitsUtils
import pyfits
import Units


bclength = 1000  #Boxcar smoothing length

def main1():
  fitsfile = "../17Vul/20120714/Star_17_Vul.fits"
  hdulist = pyfits.open(fitsfile)
  orders = FitsUtils.MakeXYpoints(hdulist[0].header, hdulist[0].data)[::-1]
  hdulist.close()

  boxcar = numpy.ones(bclength)/float(bclength)


  lines = []
  for order, index in zip(orders, range(len(orders))):
    smoothed = numpy.convolve(order.y, boxcar, mode='same')
    residuals = order.y - smoothed
    std = numpy.std(residuals)
    linepoints = numpy.where(numpy.logical_and(residuals[bclength:-bclength] - residuals.mean() < std, order.y[bclength:-bclength] > 0.9*numpy.max(order.y[bclength:-bclength])))[0] + bclength

    #Find all sets of consecutive points
    points = []
    for line in linepoints:
      if len(points) == 0 or int(line) - 1 == points[-1]:
        points.append(int(line))
      else:
        lines.append(order.x[int(numpy.median(points) + 0.5)])
        points = [int(line),]
        print lines[-1]
        yval = order.y[int(numpy.median(points) + 0.5)]
        plt.plot((lines[-1], lines[-1]), (yval-0.1, yval-0.2), 'r-')
    
    plt.figure()
    plt.title("Order number %i" %(len(orders) - index))
    plt.plot(order.x, order.y, 'k-')
    plt.plot(order.x, smoothed, 'r-')

    #plt.show()

  numpy.savetxt("Linelist.dat", lines, fmt="%.8f")


def main2():
  filename = "/Users/kgulliks/School/Research/lblrtm/run_examples/MyModel/OutputFiles/transmission-743.15-283.38-50.0-35.0-368.50-4.00-1.71-1.40"
  x,trans = numpy.loadtxt(filename, unpack=True)
  x = x[::-1]*Units.nm/Units.micron

  boxcar = numpy.ones(bclength)/float(bclength)
  smoothed = numpy.convolve(trans, boxcar, mode='same')
  residuals = trans - smoothed
  std = numpy.std(residuals[bclength:-bclength])

  #plt.plot(x, residuals - residuals[bclength:-bclength].mean())
  #plt.plot(x, smoothed)
  #plt.plot(x, (-std)*numpy.ones(x.size))
  #plt.show()

  #linepoints = numpy.where(numpy.logical_and(residuals[bclength:-bclength] - residuals.mean() < std, trans[bclength:-bclength] > 0.9*numpy.max(trans[bclength:-bclength])))[0] + bclength
  linepoints = numpy.where(residuals[bclength:-bclength] - residuals[bclength:-bclength].mean() < -std)[0] + bclength

  points = []
  lines = []
  for line in linepoints:
    #print len(points)
    if len(points) == 0 or int(line) - 1 == points[-1]:
      points.append(int(line))
    else:
      index = int(numpy.median(points) + 0.5)
      if len(points) > 1:
        minindex = trans[points[0]:points[-1]].argmin() + points[0]
      else:
        minindex = points[0]
      if trans[minindex] < 0.95 and trans[minindex] > 0.1:
        lines.append(x[minindex])
        yval = trans[minindex]
      points = [int(line),]

  #Make sure there are no points too close to each other
  tol = 0.05
  lines = sorted(lines)
  for i in range(len(lines) - 2, 0, -1):
    if numpy.abs(lines[i] - lines[i-1]) < tol:
      del lines[i]
      del lines[i-1]
    elif numpy.abs(lines[i] - lines[i+1]) < tol:
      del lines[i+1]
      del lines[i]
    else:
      index = numpy.searchsorted(x,lines[i]) - 1
      yval = trans[index]
      plt.plot((lines[i], lines[i]), (yval-0.05, yval-0.1), 'r-')
  
  plt.plot(x, trans, 'k-')
  plt.show()
  numpy.savetxt("Linelist2.dat", lines, fmt="%.8f")

if __name__ == "__main__":
  main2()
