import numpy
import matplotlib.pyplot as plt
import FitsUtils
import pyfits

bclength = 10  #Boxcar smoothing length

if __name__ == "__main__":
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
    
    plt.figure()
    plt.title("Order number %i" %(len(orders) - index))
    plt.plot(order.x, order.y, 'b-')
    plt.plot(order.x, smoothed, 'r-')

    #plt.show()

  numpy.savetxt("Linelist.dat", lines, fmt="%.8f")
