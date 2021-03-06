from scipy.interpolate import InterpolatedUnivariateSpline as interp
import scipy.signal
import os
import sys
from collections import defaultdict

import numpy as np
import DataStructures
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
from astropy import units, constants
import FittingUtilities

import StarData
import SpectralTypeRelations
from PlotBlackbodies import Planck
import RotBroad_Fast as RotBroad
import Smooth
import HelperFunctions
import Correlate
import Sensitivity



# Ensure a directory exists. Create it if not
def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)


homedir = os.environ["HOME"]
modeldir = homedir + "/School/Research/Models/Sorted/Stellar/Vband/"

#Define regions contaminated by telluric residuals or other defects. We will not use those regions in the cross-correlation
badregions = [[567.5, 575.5],
              [588.5, 598.5],
              [627, 632],
              [647, 655],
              [686, 706],
              [716, 734],
              [759, 9e9]]

vsini_values = [10, ]
vsini_values = [v * units.km.to(units.cm) for v in vsini_values]

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
              modeldir + "lte60-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted"]
""",
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
               modeldir + "lte78-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted"]"""

model_list = model_list[5:7]

modeldict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))
fullmodels = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(DataStructures.xypoint))))

model_data = []
for fname in model_list:
    if "PHOENIX2004" in fname:
        temp = int(fname.split("lte")[-1][:2]) * 100
        gravity = float(fname.split("lte")[-1][3:6])
        metallicity = float(fname.split("lte")[-1][6:10])
    elif "PHOENIX-ACES" in fname:
        temp = int(fname.split("lte")[-1][:2]) * 100
        gravity = float(fname.split("lte")[-1][3:7])
        metallicity = float(fname.split("lte")[-1][7:11])
    print "Reading in file %s" % fname
    x, y = np.loadtxt(fname, usecols=(0, 1), unpack=True)
    model = DataStructures.xypoint(x=x * units.angstrom.to(units.nm) / 1.00026, y=10 ** y)
    model2 = FittingUtilities.RebinData(model, np.linspace(model.x[0], model.x[-1], model.size()))
    model_data.append(model2)
    for vsini in vsini_values:
        #model_orders = Correlate.Process(model, orders_original, vsini, 60000.0)
        #modeldict[temp][gravity][metallicity][vsini] = model_orders
        fullmodels[temp][gravity][metallicity][vsini] = model2

if __name__ == "__main__":
    #Parse command line arguments:
    fileList = []
    extensions = True
    tellurics = False
    trimsize = 1
    windowsize = 101
    MS = SpectralTypeRelations.MainSequence()
    PMS = SpectralTypeRelations.PreMainSequence()
    vel_list = range(-400, 400, 50)
    outdir = "Sensitivity/"
    for arg in sys.argv[1:]:
        if "-e" in arg:
            extensions = False
        if "-t" in arg:
            tellurics = True  #telluric lines modeled but not removed
        else:
            fileList.append(arg)

    HelperFunctions.ensure_dir(outdir)
    outfile = open(outdir + "logfile.dat", "w")
    outfile.write("Sensitivity Analysis:\n*****************************\n\n")
    outfile.write(
        "Filename\t\t\tPrimary Temperature\tSecondary Temperature\tMass (Msun)\tMass Ratio\tVelocity\tPeak Correct?\tSignificance\n")

    for fname in fileList:
        if extensions:
            orders_original = HelperFunctions.ReadFits(fname, extensions=extensions, x="wavelength", y="flux",
                                                       errors="error")
            if tellurics:
                model_orders = HelperFunctions.ReadFits(fname, extensions=extensions, x="wavelength", y="model")
                for i, order in enumerate(orders_original):
                    orders_original[i].cont = FindContinuum.Continuum(order.x, order.y, lowreject=2, highreject=2)
                    orders_original[i].y /= model_orders[i].y

        else:
            orders_original = HelperFunctions.ReadFits(fname, errors=2)

        #Loop over orders, removing bad parts
        numorders = len(orders_original)
        for i, order in enumerate(orders_original[::-1]):
            #Linearize
            DATA = interp(order.x, order.y)
            CONT = interp(order.x, order.cont)
            ERROR = interp(order.x, order.err)
            left, right = trimsize, order.size() - 1 - trimsize
            order.x = np.linspace(order.x[left], order.x[right], order.size())
            order.y = DATA(order.x)
            order.cont = CONT(order.x)
            order.err = ERROR(order.x)


            #Remove bad regions from the data
            for region in badregions:
                left = np.searchsorted(order.x, region[0])
                right = np.searchsorted(order.x, region[1])
                if left > 0 and right < order.size():
                    print "Warning! Bad region covers the middle of order %i" % i
                    print "Removing full order!"
                    left = 0
                    right = order.size()
                order.x = np.delete(order.x, np.arange(left, right))
                order.y = np.delete(order.y, np.arange(left, right))
                order.cont = np.delete(order.cont, np.arange(left, right))
                order.err = np.delete(order.err, np.arange(left, right))

            #Remove whole order if it is too small
            remove = False
            if order.x.size <= windowsize:
                remove = True
            else:
                velrange = 3e5 * (np.median(order.x) - order.x[0]) / np.median(order.x)
                if velrange <= 1000.0:
                    remove = True
            if remove:
                print "Removing order %i" % (numorders - 1 - i)
                orders_original.pop(numorders - 1 - i)
            else:
                order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=1, lowreject=1.5, highreject=5)
                orders_original[numorders - 1 - i] = order.copy()



        #Read in the name of the star from the fits header
        header = pyfits.getheader(fname)
        starname = header["OBJECT"]
        print starname

        #Get spectral type of the primary from the name and simbad
        stardata = StarData.GetData(starname)
        primary_temp = [MS.Interpolate(MS.Temperature, stardata.spectype[:2]), ]
        #age = 'MS'   #Play with this later!
        primary_mass = MS.Interpolate(MS.Mass, stardata.spectype[:2])
        age = PMS.GetMainSequenceAge(primary_mass)


        #Check for close companions
        companions = HelperFunctions.CheckMultiplicityWDS(starname)
        if companions:
            for configuration in companions:
                component = companions[configuration]
                if component["Separation"] < 3.0 and component["Secondary SpT"] != "Unknown":
                    print "Known %s companion with a separation of %g arcseconds!" % (
                    component["Secondary SpT"], component["Separation"])
                    primary_temp.append(MS.Interpolate(MS.Temperature, component["Secondary SpT"]))
        companion = HelperFunctions.CheckMultiplicitySB9(starname)
        if companion:
            if companion['K1'] != 'Unknown' and companion['K2'] != 'Unknown':
                mass = primary_mass * companion['K1'] / companion['K2']
                spt = MS.GetSpectralType(MS.Mass, mass, interpolate=True)
                primary_temp.append(MS.Interpolate(MS.Temperature, spt))

        # Process all of the model spectra up front.
        # As long as you do more than one star at once, this is faster
        for temp in sorted(fullmodels.keys()):
            for gravity in sorted(fullmodels[temp].keys()):
                for metallicity in sorted(fullmodels[temp][gravity].keys()):
                    for vsini in vsini_values:
                        if modeldict[temp][gravity][metallicity][vsini] == []:
                            model_orders = Correlate.Process(model, orders_original, vsini, 60000, debug=True)
                            modeldict[temp][gravity][metallicity][vsini] = model_orders
                        else:
                            model_orders = modeldict[temp][gravity][metallicity][vsini]
                        model = fullmodels[temp][gravity][metallicity][vsini]

                        #Get info about the secondary star for this model temperature
                        secondary_spt = MS.GetSpectralType(MS.Temperature, temp)
                        secondary_radius = PMS.Interpolate(secondary_spt, age, key='Radius')
                        secondary_mass = PMS.Interpolate(secondary_spt, age, key='Mass')
                        massratio = secondary_mass / primary_mass


                        #Check sensitivity to this star
                        orders = [order.copy() for order in orders_original]  #Make a copy of orders
                        output_dir = "%s/Sensitivity/" % (os.path.dirname(fname))
                        if output_dir.startswith("/"):
                            output_dir = output_dir[1:]
                        found, significance = Sensitivity.Analyze(orders,
                                                                  model,
                                                                  prim_temp=primary_temp,
                                                                  sec_temp=temp,
                                                                  age=age,
                                                                  smoothing_windowsize=101,
                                                                  smoothing_order=5,
                                                                  outdir=output_dir,
                                                                  outfilebase=fname.split(".fits")[0],
                                                                  process_model=False,
                                                                  model_orders=model_orders,
                                                                  debug=False)

            #Write to logfile
            vels = range(-400, 450, 50)
            for v, f, s in zip(vels, found, significance):
                if f:
                    #Signal found!
                    outfile.write("%s\t%i\t\t\t%i\t\t\t\t%.2f\t\t%.4f\t\t%i\t\tyes\t\t%.2f\n" % (
                    fname, primary_temp[0], temp, secondary_mass, massratio, v, s))
                else:
                    outfile.write("%s\t%i\t\t\t%i\t\t\t\t%.2f\t\t%.4f\t\t%i\t\tno\t\tN/A\n" % (
                    fname, primary_temp[0], temp, secondary_mass, massratio, v))

    outfile.close()
      
    


