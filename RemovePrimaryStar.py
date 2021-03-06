#!/opt/local/bin/python
import os
import sys
from scipy.interpolate import UnivariateSpline
from collections import defaultdict
import time

import numpy as np
import pylab
import DataStructures
import MakeModel

import SpectralTypeRelations
from PlotBlackbodies import Planck
import Units
import Correlate
import RotBroad
import FindContinuum
import FitTellurics_McDonald as Outputter


"""
   This function removes a primary model from data  
"""


def GetModel(data, model, vel=0.0, vsini=15 * Units.cm / Units.km):
    # Broaden
    model = RotBroad.Broaden(model, vsini=vsini, findcont=True)

    #Interpolate Model, and note the endpoints
    first = model.x[0] * (1. + vel / Units.c)
    last = model.x[-1] * (1. + vel / Units.c)
    a = []
    for i in range(1, model.x.size):
        a.append(model.x[i] - model.x[i - 1])
    model_fcn = UnivariateSpline(model.x, model.y, s=0)

    data2 = []
    for order in data:
        data2.append(order.copy())

    ##############################################################
    #Begin main loop over the orders
    ##############################################################
    output_orders = []
    for i in range(len(data2)):
        order = data2[i]
        order.cont = FindContinuum.Continuum(order.x, order.y, lowreject=2, highreject=5)

        model2 = DataStructures.xypoint(x=order.x, y=model_fcn(order.x * (1. + vel / Units.c)))

        #Get model continuum in this section
        model2.cont = FindContinuum.Continuum(model2.x, model2.y, lowreject=1, fitorder=3)

        #Reduce resolution
        model2 = FittingUtilities.ReduceResolution(model2.copy(), 60000)

        #Fit velocity with a straight shift via cross-correlation
        ycorr = np.correlate(order.y / order.cont - 1.0, model2.y / model2.cont - 1.0, mode="full")
        xcorr = np.arange(ycorr.size)
        lags = xcorr - (order.x.size - 1)
        distancePerLag = (order.x[-1] - order.x[0]) / float(order.x.size)
        offsets = -lags * distancePerLag
        offsets = offsets[::-1]
        ycorr = ycorr[::-1]
        fit = np.poly1d(np.polyfit(offsets, ycorr, 3))
        ycorr = ycorr - fit(offsets)
        maxindex = ycorr.argmax()
        model2 = DataStructures.xypoint(x=order.x, y=model_fcn(order.x * (1. + vel / Units.c) + offsets[maxindex]))
        model2.cont = FindContinuum.Continuum(model2.x, model2.y, lowreject=1, fitorder=3)
        model2 = FittingUtilities.ReduceResolution(model2.copy(), 60000)

        #Scale using Beer's Law
        line_indices = np.where(model2.y / model2.cont < 0.96)[0]
        #w line_indices = np.where(order.y / order.cont < 0.90)[0]
        if len(line_indices) > 0:
            #scale_fcn = np.poly1d(np.polyfit(model2.x[line_indices], np.log(order.y[line_indices]/order.cont[line_indices]) / np.log(model2.y[line_indices] / model2.cont[line_indices]), 2))
            scale = np.median(np.log(order.y[line_indices] / order.cont[line_indices]) / np.log(
                model2.y[line_indices] / model2.cont[line_indices]))
            #pylab.plot(model2.x[line_indices], np.log(order.y[line_indices]/order.cont[line_indices]) / np.log(model2.y[line_indices] / model2.cont[line_indices]))
            #pylab.plot(model2.x[line_indices], scale_fcn(model2.x[line_indices]))
            model2.y = model2.y ** scale
            model2.cont = model2.cont ** scale
            #model2.y = model2.y**scale_fcn(model2.x)
            #model2.cont = model2.cont**scale_fcn(model2.x)

        print "Order %i: Scale = %g" % (i, scale)

        output_orders.append(model2.copy())

    pylab.show()
    return output_orders


if __name__ == "__main__":
    import FitsUtils
    import os
    import sys

    home = os.environ["HOME"]
    try:
        datafile = sys.argv[1]
        print "Removing primary star from  ", datafile
    except IndexError:
        print "Error! Must give .fits file!"
        sys.exit()

    modeldir = os.environ["HOME"] + "/School/Research/Models/Sorted/Stellar/Vband/"
    files = os.listdir(modeldir)
    modelfiles = defaultdict(list)
    for fname in files:
        try:
            temperature = float(fname.split("lte")[-1].split("-")[0]) * 100
        except ValueError:
            try:
                temperature = float(fname.split("lte")[-1].split("+")[0]) * 100
            except ValueError:
                print "Skipping file %s" % fname
        modelfiles[temperature].append(modeldir + fname)

    # Read in data
    orders_original = tuple(FitsUtils.MakeXYpoints(datafile, extensions=True, x="wavelength", y="flux", errors="error"))
    orders_original = orders_original[::-1]

    #Check for command line arguments
    p_spt = "A0"
    vsini = 15.0
    vel = 0.0
    metal = 0.0
    if len(sys.argv) > 2:
        for arg in sys.argv[2:]:
            if "primary" in arg:
                p_spt = arg.split("=")[-1]
            elif "vsini" in arg:
                vsini = float(arg.split("=")[-1])
            elif "vel" in arg:
                vel = float(arg.split("=")[-1])
            elif "metal" in arg:
                metal = float(arg.split("=")[-1])

    #Get the best logg and temperature for a main sequence star with the given spectral type
    MS = SpectralTypeRelations.MainSequence()
    p_temp = MS.Interpolate(MS.Temperature, p_spt)
    p_mass = MS.Interpolate(MS.Mass, p_spt)
    radius = MS.Interpolate(MS.Radius, p_spt)
    logg = np.log10(Units.G * p_mass * Units.Msun / (radius * Units.Rsun) ** 2)

    #Find the best fit temperature
    temperature = modelfiles.keys()[np.argmin((modelfiles.keys() - p_temp) ** 2)]
    print p_temp, temperature

    #Find the best logg
    best_logg = 9e9
    indices = []
    for fname in modelfiles[temperature]:
        logg_val = float(fname.split("lte")[-1].split("-")[1][:3])
        if np.abs(logg_val - logg) < np.abs(best_logg - logg):
            best_logg = logg_val
    for i, fname in enumerate(modelfiles[temperature]):
        logg_val = float(fname.split("lte")[-1].split("-")[1][:3])
        if logg_val == best_logg:
            indices.append(i)
    filenames = []
    for i in indices:
        filenames.append(modelfiles[temperature][i])

    chisquareds = []
    models = []
    for modelfile in filenames:
        #Read in model
        print "Model file: %s" % modelfile
        x, y = np.loadtxt(modelfile, usecols=(0, 1), unpack=True)
        x = x * Units.nm / Units.angstrom
        y = 10 ** y
        model = DataStructures.xypoint(x=x, y=y)

        orders = GetModel(list(orders_original), model, vel=vel * Units.cm / Units.km,
                          vsini=vsini * Units.cm / Units.km)
        chisq = 0.0
        i = 0
        for model, data in zip(orders, orders_original):
            if i != 15 and i != 16 and i != 44:
                data.cont = FindContinuum.Continuum(data.x, data.y, lowreject=2, highreject=2)
                chisq += np.sum((data.y - model.y * data.cont / model.cont) ** 2 / data.err ** 2) / float(data.size())
                i += 1
        chisquareds.append(chisq)
        models.append(orders)

    best_index = np.argmin(chisquareds)
    orders = models[best_index]
    print chisquareds

    outfilename = "%s-0.fits" % (datafile.split(".fits")[0])
    print "Outputting to %s" % outfilename

    orders_original = orders_original[::-1]
    orders = orders[::-1]

    for i, original in enumerate(orders_original[2:-1]):
        #original = orders_original[i+2]
        original.cont = FindContinuum.Continuum(original.x, original.y, lowreject=2, highreject=2)
        model = orders[i + 2]
        pylab.figure(1)
        pylab.plot(original.x, original.y / original.cont, 'k-')
        pylab.plot(model.x, model.y / model.cont, 'r-')
        pylab.figure(2)
        pylab.plot(original.x, original.y / (original.cont * model.y / model.cont), 'k-')
        original.y /= model.y / model.cont

        columns = {"wavelength": original.x,
                   "flux": original.y,
                   "continuum": original.cont,
                   "error": original.err}

        #Output
        if i == 0:
            Outputter.OutputFitsFile(columns, datafile, outfilename, mode="new")
        else:
            Outputter.OutputFitsFile(columns, outfilename, outfilename)

    pylab.show()

