import Bootstrap
import HelperFunctions
import FittingUtilities
import matplotlib.pyplot as plt
import os
import sys
import numpy
from scipy.stats.mstats import mquantiles


homedir = os.environ["HOME"]
modeldir = homedir + "/School/Research/Models/Sorted/Stellar/Vband/"


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


temp_list = []
gravity_list = []
metal_list = []
for fname in model_list:
  if "PHOENIX2004" in fname:
    temp = int(fname.split("lte")[-1][:2])*100
    gravity = float(fname.split("lte")[-1][3:6])
    metallicity = float(fname.split("lte")[-1][6:10])
  elif "PHOENIX-ACES" in fname:
    temp = int(fname.split("lte")[-1][:2])*100
    gravity = float(fname.split("lte")[-1][3:7])
    metallicity = float(fname.split("lte")[-1][7:11])
  temp_list.append(temp)
  gravity_list.append(gravity)
  metal_list.append(metallicity)
temp_list = numpy.array(temp_list)
gravity_list = numpy.array(gravity_list)
metal_list = numpy.array(metal_list)



if __name__ == "__main__":
  # Parse command line arguments
  fileList = []
  T = 3500
  logg = 4.0
  metal = 0.0
  Nboot = 10
  for arg in sys.argv[1:]:
    if "-T" in arg:
      T = float(arg.split("=")[1])
    elif "-logg" in arg:
      logg = float(arg.split("=")[1])
    elif "-metal" in arg:
      metal = float(arg.split("=")[1])
    elif "-N" in arg:
      Nboot = int(arg.split("=")[1])
    else:
      fileList.append(arg)
  
  
  #Find the model with the closest T, logg, and metallicity
  idx = numpy.argmin(5*(T-temp_list)**2 + 
                    2*(metal-metal_list)**2 + 
                    1*(logg-gravity_list)**2)
  T = temp_list[idx]
  logg = gravity_list[idx]
  metal = metal_list[idx]
  
  
  #Now, loop over infiles
  for fname in fileList:
    orders = HelperFunctions.ReadFits(fname, 
                                      extensions=True, 
                                      x="wavelength",
                                      y="flux",
                                      errors="error",
                                      cont="continuum")[:50]
    for order in orders:
      order.cont = FittingUtilities.Continuum(order.x, order.y, lowreject=3, highreject=3)
    
    
    #Get bootstrap samples
    y = Bootstrap.GetSamples(orders, model_list[idx].split("/")[-1], Nboot, vsini=10, resolution=60000.0)
    
    numpy.savetxt("Bootstrap_samples.dat", numpy.transpose((y,)) )
    quantiles = mquantiles(y)
    print "The mean CCF height is %g" %numpy.mean(y)
    print "The median CCF height is %g" %quantiles[1]
    print "The standard deviation of the peak height is %g" %numpy.std(y)
    print "The 25% quantile is %g" %quantiles[0]
    print "The 75% quantile is %g" %quantiles[2]
    plt.hist(y)
    plt.show()
  
  
  
  
  
  
  
  
  
  
  
  
