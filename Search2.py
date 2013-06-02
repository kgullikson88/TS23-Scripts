import Correlate
import FitsUtils
import numpy
from scipy.interpolate import InterpolatedUnivariateSpline as interp
import os
import sys
import DataStructures
import FindContinuum
import matplotlib.pyplot as plt
#import Units
from astropy import units, constants


homedir = os.environ["HOME"]
modeldir = homedir + "/School/Research/Models/Sorted/Stellar/Vband/"

#Define regions contaminated by telluric residuals or other defects. We will not use those regions in the cross-correlation
badregions = [[588.8, 589.9],
              [627.1, 635.4]]
badregions = [[0, 466],
              [567.5, 575.5],
              [587.5, 593],
              [627, 634.5],
              [686, 706],
              [716, 742],
              [759, 775]]

#Set up model list
"""
model_list = [ modeldir + "lte30-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte32-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte34-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte35-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte36-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte37-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte38-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte39-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte40-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte42-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte44-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte46-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte48-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte50-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte51-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte52-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte53-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte54-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte55-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte56-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte57-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte58-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte59-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte60-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte61-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte62-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte63-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
"""
model_list = [ modeldir + "lte64-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte65-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte66-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte67-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte68-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte69-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte69-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted"]
"""               modeldir + "lte70-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte70-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte72-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte74-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte74-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte76-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte78-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte30-4.0-0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte30-4.0+0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte31-4.0-0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte31-4.0+0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte32-4.0-0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte32-4.0+0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte33-4.0-0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte33-4.0+0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte34-4.0-0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte34-4.0+0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte35-4.0-0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte35-4.0+0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte36-4.0-0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte36-4.0+0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte37-4.0-0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte37-4.0+0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte38-4.0-0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte38-4.0+0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte39-4.0-0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte39-4.0+0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte40-4.0-0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte40-4.0+0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte41-4.0-0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte41-4.0+0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte42-4.0-0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte42-4.0+0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte43-4.0-0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte43-4.0+0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte44-4.0-0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte44-4.0+0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte45-4.0-0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte45-4.0+0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte46-4.0-0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte46-4.0+0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte47-4.0-0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte47-4.0+0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte48-4.0-0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte48-4.0+0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte49-4.0-0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte49-4.0+0.5.Cond.PHOENIX2004.tab.7.sorted",
               modeldir + "lte50-3.5-0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte50-4.0+0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte51-4.0-0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte51-4.0+0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte52-4.0-0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte52-4.0+0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte53-4.0+0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte54-4.0-0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte54-4.0+0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte55-4.0-0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte55-4.0+0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte56-3.5-0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte56-4.0+0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte57-4.0-0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte57-4.0+0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte58-4.0-0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte58-4.0+0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte59-4.0+0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte60-4.0+0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte61-4.0-0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte61-4.0+0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte62-4.0+0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte63-4.0-0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte63-4.0+0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte64-4.0-0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte64-4.0+0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte65-4.0+0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte66-4.0-0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte66-4.0+0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte67-4.0+0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte68-4.0-0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte68-4.0+0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte69-4.0-0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte69-4.0+0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte70-3.5-0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte70-4.0+0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte72-3.5-0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte72-4.0+0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte74-4.0+0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte76-4.0+0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte78-3.5-0.5.Cond.PHOENIX2004.direct.7.sorted",
               modeldir + "lte78-4.0+0.5.Cond.PHOENIX2004.direct.7.sorted"]
"""               
               
star_list = []
temp_list = []
gravity_list = []
metal_list = []
model_data = []
for fname in model_list:
  if "PHOENIX2004" in fname:
    temp = int(fname.split("lte")[-1][:2])*100
    gravity = float(fname.split("lte")[-1][3:6])
    metallicity = float(fname.split("lte")[-1][6:10])
  elif "PHOENIX-ACES" in fname:
    temp = int(fname.split("lte")[-1][:2])*100
    gravity = float(fname.split("lte")[-1][3:7])
    metallicity = float(fname.split("lte")[-1][7:11])
  print "Reading in file %s" %fname
  x,y = numpy.loadtxt(fname, usecols=(0,1), unpack=True)
  model_data.append( DataStructures.xypoint(x=x*units.angstrom.to(units.nm), y=10**y) )
  star_list.append(str(temp))
  temp_list.append(temp)
  gravity_list.append(gravity)
  metal_list.append(metallicity)
  

if __name__ == "__main__":
  #Parse command line arguments:
  fileList = []
  extensions=True
  tellurics=False
  trimsize = 100
  for arg in sys.argv[1:]:
    if "-e" in arg:
      extensions=False
    if "-t" in arg:
      tellurics=True  #telluric lines modeled but not removed
    else:
      fileList.append(arg)

  for fname in fileList:
    if extensions:
      orders = FitsUtils.MakeXYpoints(fname, extensions=extensions, x="wavelength", y="flux", errors="error")
      if tellurics:
        model_orders = FitsUtils.MakeXYpoints(fname, extensions=extensions, x="wavelength", y="model")
        for i, order in enumerate(orders):
          orders[i].cont = FindContinuum.Continuum(order.x, order.y, lowreject=2, highreject=2)
          orders[i].y /= model_orders[i].y
          
    else:
      orders = FitsUtils.MakeXYpoints(fname, errors=2)
    numorders = len(orders)
    for i, order in enumerate(orders[::-1]):
      DATA = interp(order.x, order.y)
      CONT = interp(order.x, order.cont)
      ERROR = interp(order.x, order.err)
      order.x = numpy.linspace(order.x[trimsize], order.x[-trimsize], order.size() - 2*trimsize)
      order.y = DATA(order.x)
      order.cont = CONT(order.x)
      order.err = ERROR(order.x)
      
      #Remove bad regions from the data
      for region in badregions:
        left = numpy.searchsorted(order.x, region[0])
        right = numpy.searchsorted(order.x, region[1])
        order.x = numpy.delete(order.x, numpy.arange(left, right))
        order.y = numpy.delete(order.y, numpy.arange(left, right))
        order.cont = numpy.delete(order.cont, numpy.arange(left, right))
        order.err = numpy.delete(order.err, numpy.arange(left, right))

      #Remove whole order if it is too small
      if order.x.size > 10:
        order.cont = FindContinuum.Continuum(order.x, order.y, lowreject=3, highreject=3)
        orders[numorders -1 -i] = order.copy()
      else:
        orders.pop(numorders - 1 - i)

    """
    for i, order in enumerate(orders):
      plt.plot(order.x, order.y/order.cont+i)
    plt.show()
    sys.exit()
    """
        
    #Do the cross-correlation
    for vsini in [10, 20, 30, 40]:
      Correlate.PyCorr2(orders, resolution=60000, outdir="Cross_correlations/%s" %(fname.split(".fits")[0]), models=model_data, stars=star_list, temps=temp_list, gravities=gravity_list, metallicities=metal_list, vsini=vsini*units.km.to(units.cm), debug=True, outfilebase=fname)


