import numpy
import scipy
import pylab
import FitTellurics
import MakeModel
import FitsUtils
import pyfits
import Units
import sys
import os


homedir = os.environ['HOME']
LineListFile = homedir + "/School/Research/McDonaldData/LineList.dat"


#Do a cross-correlation first, to get the wavelength solution close
def CCImprove(order, model):
  ycorr = scipy.correlate(order.y-1.0, model.y-1.0, mode="full")
  xcorr = numpy.arange(ycorr.size)
  maxindex = ycorr.argmax()
  lags = xcorr - (order.y.size-1)
  distancePerLag = (order.x[-1] - order.x[0])/float(order.x.size)
  offsets = -lags*distancePerLag
  print "maximum offset: ", offsets[maxindex], " nm"

  if numpy.abs(offsets[maxindex]) < 1.0:
    #Apply offset
    order.x = order.x + offsets[maxindex]
  return order


if __name__ == "__main__":
  fitsfile = sys.argv[1]
  hdulist = pyfits.open(fitsfile)
  orders = FitsUtils.MakeXYpoints(hdulist[0].header, hdulist[0].data)[::-1][3:-2]

  resolution = 60000.0
  angle = float(hdulist[0].header["zd"])
  humidity = 40.0
  temperature = (64 - 32)*5./9. + 273.15
  pressure_start = 23.55*Units.hPa/Units.inch_Hg
  pressure_end = 23.58*Units.hPa/Units.inch_Hg
  pressure = (pressure_start + pressure_end)/2.0
  wave_start = 400
  wave_end = 1100
  wavenum_start = int(1.0e7/wave_end)
  wavenum_end = int(1.0e7/wave_start+1)
  ch4=1.6
  co2 = 368.5
  co=0.14
  o3=4e-2
  o2=2.12e5
  const_pars = [temperature, pressure, co2, o3, wavenum_start, wavenum_end, resolution, angle, ch4, co, 5]
  pars = [humidity, o2]

  for order in orders:
    order.x *= Units.nm/Units.angstrom
    order.cont = numpy.ones(order.x.size)
    order.err = numpy.sqrt(order.y)
  hdulist.close()

  #Read in line list:
  linelist = numpy.loadtxt(LineListFile)

  
  for order in orders[-2:]:
    wave_start = order.x[0] - 10.
    wave_end = order.x[-1] + 10.
    wavenum_start = int(1.0e7/wave_end)
    wavenum_end = int(1.0e7/wave_start+1)
    const_pars[4] = wavenum_start
    const_pars[5] = wavenum_end
    model = FitTellurics.FitFunction(order.copy(), pars, const_pars)

    model = MakeModel.ReduceResolution(model.copy(), resolution)
    model = MakeModel.RebinData(model.copy(), order.x.copy())

    order = FitTellurics.FitContinuum3(order, model, 4)

    order = CCImprove(order, model)
    pylab.plot(order.x, order.y)
    pylab.plot(model.x, model.y)
    pylab.show()

