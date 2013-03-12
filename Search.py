import Correlate
import FitsUtils
import numpy
from scipy.interpolate import InterpolatedUnivariateSpline as interp
import os
import sys
import DataStructures
import FindContinuum
import matplotlib.pyplot as plt

homedir = os.environ["HOME"]
modeldir = homedir + "/School/Research/Models/Sorted/Stellar/Vband/"

#Define regions contaminated by telluric residuals or the picket fence. We will not use those regions in the cross-correlation
badregions = [[0,389],
              [454,483],
              [626,632],
              [685,696],
              [715,732]]

#Set up model list
model_list = [modeldir + "lte30-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
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
               modeldir + "lte64-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte65-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte66-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte67-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte68-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte69-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte69-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte70-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte70-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte72-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte74-4.00-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte74-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte76-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
               modeldir + "lte78-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted"]
star_list = []
for fname in model_list:
  temp = int(fname.split("lte")[-1][:2])*100
  star_list.append(str(temp))

if __name__ == "__main__":
  #Parse command line arguments:
  fileList = []
  extensions=True
  for arg in sys.argv[1:]:
    if "-e" in arg:
      extensions=False
    else:
      fileList.append(arg)

  for fname in fileList:
    if extensions:
      orders = FitsUtils.MakeXYpoints(fname, extensions=extensions, x="wavelength", y="flux", errors="error")
    else:
      orders = FitsUtils.MakeXYpoints(fname, errors=2)

    xspacing = 1e9
    for i, order in enumerate(orders):
      orders[i].cont = FindContinuum.Continuum(order.x, order.y, lowreject=2, highreject=2)
      spacing = (order.x[-1] - order.x[0]) / float(order.size() )
      if spacing < xspacing:
        xspacing = spacing
      
    data = DataStructures.CombineXYpoints(orders, xspacing=xspacing)

    #Remove bad regions from the data
    for region in badregions:
      left = numpy.searchsorted(data.x, region[0])
      right = numpy.searchsorted(data.x, region[1])
      data.y[left:right] = data.cont[left:right]
    plt.plot(data.x, data.y/data.cont)
    plt.show()

    #Do the cross-correlation
    Correlate.PyCorr(data, combine=False, resolution=60000, outdir="Cross_correlations/%s" %(fname.split(".fits")[0]), models=model_list, stars=star_list)


#This will do the correlation within python/numpy
#The combine keyword decides whether to combine the chips into a master cross-correlation 
#The normalize keyword decides whether to output as correlation power, or as significance
#The sigmaclip keyword decides whether to perform sigma-clipping on each chip before cross-correlating
#The nsigma keyword tells the program how many sigma to clip. This is ignored if sigmaclip = False
#The clip_order keyword tells what order polynomial to fit the flux to during sigma clipping. Ignored if sigmaclip = False
#The models keyword is a list of models to cross-correlate against (either filenames of two-column ascii files, or
#    each entry should be a list with the first element the x points, and the second element the y points
#The segments keyword controls which orders of the data to use, and which parts of them. Can be used to ignore telluric
#    contamination. Can be a string (default) which will use all of the orders, a list of integers which will
#    use all of the orders given in the list, or a dictionary of lists which gives the segments of each order to use.
#save_output determines whether the cross-correlation is saved to a file, or just returned
#def PyCorr(filename, combine=True, normalize=False, sigmaclip=False, nsigma=3, clip_order=3, stars=star_list, temps=temp_list, models=model_list, corr_mode='valid', process_model=True, vsini=15*Units.cm/Units.km, resolution=100000, segments="all", save_output=True, outdir=outfiledir, outfilename=None)
