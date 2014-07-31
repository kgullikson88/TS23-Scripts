#!/opt/local/bin/python
from scipy.interpolate import UnivariateSpline
from collections import defaultdict

import numpy
import pylab
import DataStructures

import SpectralTypeRelations
from PlotBlackbodies import Planck
import Units
import RotBroad
import FindContinuum


"""
   This function adds a given model spectrum to the data at the requested vsini and
     detector resolution. Useful for checking cross-correlation results to see if you
     should be able to visibly distinguish the secondary lines   
"""

homedir = os.environ['HOME'] + "/"
outfiledir = os.getcwd() + "/Sensitivity/"


def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)


ensure_dir(outfiledir)


# Define the sections of each order to use (those without telluric contamination)
good_sections = {1: [[-1, -1], ],
                 2: [[-1, -1], ],
                 3: [[-1, -1], ],
                 4: [[-1, -1], ],
                 5: [[-1, -1], ],
                 6: [[-1, -1], ],
                 7: [[-1, -1], ],
                 8: [[-1, -1], ],
                 9: [[-1, -1], ],
                 10: [[-1, -1], ],
                 11: [[-1, -1], ],
                 12: [[-1, -1], ],
                 13: [[-1, -1], ],
                 14: [[-1, -1], ],
                 15: [[-1, -1], ],
                 16: [[-1, -1], ],
                 17: [[-1, -1], ],
                 18: [[-1, -1], ],
                 19: [[-1, 686.7], ],
                 20: [[-1, 1e9], ],
                 21: [[661.4, 667.3], ],
                 22: [[-1, -1], ],
                 23: [[-1, 1e9]],
                 24: [[-1, 627.6], ],
                 25: [[-1, 1e9], ],
                 26: [[-1, 1e9], ],
                 27: [[-1, -1], ],
                 28: [[-1, -1], ],
                 29: [[-1, 1e9], ],
                 30: [[-1, 567.6], ],
                 31: [[-1, 1e9], ],
                 32: [[-1, 1e9], ],
                 33: [[-1, 1e9], ],
                 34: [[-1, 1e9], ],
                 35: [[-1, 1e9], ],
                 36: [[-1, 1e9], ],
                 37: [[-1, 1e9], ],
                 38: [[-1, 501.1], ],
                 39: [[-1, 491.6], [492.8, 1e9]],
                 40: [[-1, 1e9], ],
                 41: [[-1, 1e9], ],
                 42: [[-1, -1], ],
                 43: [[-1, -1], ],
                 44: [[-1, -1], ],
                 45: [[-1, 1e9], ],
                 46: [[-1, 1e9], ],
                 47: [[-1, 1e9], ],
                 48: [[-1, 1e9], ],
                 49: [[-1, 1e9], ],
                 50: [[-1, 1e9], ],
                 51: [[-1, 1e9], ],
                 52: [[-1, -1], ],
                 43: [[-1, -1], ]}

#Define regions contaminated by telluric residuals or the picket fence. We will not use those regions in the cross-correlation
badregions = [[0, 389],
              [454, 483],
              [626, 632],
              [685, 696],
              [715, 732]]


def Add(data, model, prim_spt, sec_spt, age="MS", vel=0.0, SNR=1e6, SN_order=19, sensitivity=lambda x: 1.0,
        vsini=15.0 * Units.cm / Units.km):
    """
      This function will add a model to the data. The flux ratio is determined
      from the primary spectral type (prim_spt), the secondary spectral type
      (sec_spt), and the age of the system.
      The age keyword argument can be either a string, in which case the main sequence
             spectral-type - radius relations are used, or a number in which case the
             radii of the two stars are determined from evolutionary tracks (NOT YET IMPLEMENTED)
      The 'vel' keyword gives a radial velocity at which the model should be added (MUST be given in cm/s)
      The SNR keyword gives the desired signal to noise of the spectrum, in order SN_order. Note
             that it is only possible to create a spectrum with signal to noise LOWER than that
             in the data. SN_order should be given in fortran-style numbering (starting at 1, not 0)
    """
    if type(age) == str:
        #Main sequence stars!
        MS = SpectralTypeRelations.MainSequence()
        prim_radius = MS.Interpolate(MS.Radius, prim_spt)
        prim_temp = MS.Interpolate(MS.Temperature, prim_spt)
        sec_radius = MS.Interpolate(MS.Radius, sec_spt)
        sec_temp = MS.Interpolate(MS.Temperature, sec_spt)
    else:
        sys.exit("Sorry! Can only handle Main Sequence stars right now!")


    #Interpolate Model, and note the endpoints
    first = model.x[0] * (1. + vel / Units.c)
    last = model.x[-1] * (1. + vel / Units.c)
    a = []
    for i in range(1, model.x.size):
        a.append(model.x[i] - model.x[i - 1])
    print min(a)
    model_fcn = UnivariateSpline(model.x, model.y, s=0)

    data2 = []
    for order in data:
        data2.append(order.copy())

    #Figure out how to add noise to the spectrum to get the desired S/N
    flux = Planck(data2[SN_order - 1].x * Units.cm / Units.nm, prim_temp).mean()
    SN_factor = SNR / numpy.sqrt(flux * sensitivity(data2[SN_order - 1].x.mean()))

    ##############################################################
    #Begin main loop over the orders
    output_orders = []
    for i in range(len(data2)):
        order = data2[i]
        order.cont = FindContinuum.Continuum(order.x, order.y, lowreject=2, highreject=5)
        prim_flux = Planck(order.x * Units.cm / Units.nm, prim_temp) * prim_radius ** 2
        sec_flux = Planck(order.x * Units.cm / Units.nm, sec_temp) * sec_radius ** 2
        fluxratio = sec_flux / prim_flux
        print "order %i flux ratio = %.3g" % (i + 1, numpy.mean(fluxratio))

        model2 = DataStructures.xypoint(x=order.x, y=model_fcn(order.x * (1. + vel / Units.c)))

        #Get model continuum in this section
        model2.cont = FindContinuum.Continuum(model2.x, model2.y)

        #Rotationally broaden
        #model2 = RotBroad.Broaden(model2, vsini)

        #Reduce resolution
        model2 = FittingUtilities.ReduceResolution(model2.copy(), 60000)

        #fluxratio = 0.5

        #Scale the model by the above scale factor and normalize
        scaled_model = (model2.y / model2.cont - 1.0) * fluxratio + 1.0
        model2.y = scaled_model
        model2.cont = numpy.ones(model2.size())
        output_orders.append(model2.copy())

        pylab.plot(order.x, order.y / order.cont, 'k-')
        pylab.plot(model2.x, scaled_model + 0.01, 'r-')


        #Add noise to the data
        noise = numpy.random.normal(loc=0, scale=1.0 / (
        SN_factor * numpy.sqrt(numpy.mean(prim_flux * sensitivity(order.x.mean()) / prim_radius ** 2))),
                                    size=data2[i].x.size)
        #data2[i].y += noise


        #Add model to the data
        data2[i].y = (scaled_model) * order.cont + order.y
        data2[i].cont = FindContinuum.Continuum(data2[i].x, data2[i].y, lowreject=2, highreject=5)

        #pylab.plot(data2[i].x, data2[i].y - 0.02, 'b-')

    #pylab.xlabel("Wavelength (nm)")
    #pylab.ylabel("Normalized Flux")
    #pylab.ylim((0.97, 1.01))
    #pylab.xlim((510, 570))
    pylab.show()
    #sys.exit()
    return data2


if __name__ == "__main__":
    import FitsUtils
    import os
    import sys

    home = os.environ["HOME"]
    try:
        datafile = sys.argv[1]
        print "Using ", datafile, "as template"
    except IndexError:
        print "Error! Must give .fits file!"
        sys.exit()
    tolerance = 5  #allow highest cross-correlation peak to be anywhere within 5 km/s of the correct velocity
    MS = SpectralTypeRelations.MainSequence()  #Class for interpolating stellar info from the spectral type

    #Make some lists that we will loop through in the analysis
    parent_spts = ["B0", "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9",
                   "A0", "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9"]
    modeldir = homedir + "School/Research/Models/Sorted/Stellar/Vband/"
    files = os.listdir(modeldir)
    modelfiles = defaultdict(list)
    for fname in files:
        if "PHOENIX2004" in fname:
            temperature = int(fname.split("lte")[-1][:2]) * 100
            metallicity = float(fname.split("lte")[-1][6:10])
        elif "PHOENIX-ACES" in fname:
            temperature = int(fname.split("lte")[-1][:2]) * 100
            metallicity = float(fname.split("lte")[-1][7:11])
        else:
            continue

        if numpy.abs(metallicity) < 0.1:
            modelfiles[temperature].append(modeldir + fname)

    #Read in data
    orders_original = tuple(FitsUtils.MakeXYpoints(datafile, extensions=True, x="wavelength", y="flux", errors="error"))
    #orders_original = tuple(FitsUtils.MakeXYpoints(datafile, errors=2))
    orders_original = orders_original[::-1]

    #Check for command line arguments specifying spectral type endpoints or the logfile name
    logfilename = outfiledir + "summary.dat"
    sensitivity_fcn = lambda x: 1.0
    p_left = 0
    p_right = len(parent_spts)
    s_temp = 5800
    vsini = 10.0
    vel = 0.0
    p_spt = "A0"
    if len(sys.argv) > 2:
        for arg in sys.argv[2:]:
            if "primary" in arg:
                p_spt = arg.split("=")[-1]
            elif "secondary" in arg:
                try:
                    s_temp = float(arg.split("=")[-1])
                except ValueError:
                    s_spt = arg.split("=")[-1]
                    s_temp = MS.Interpolate(MS.Temperature, s_spt)
            elif "vsini" in arg:
                vsini = float(arg.split("=")[-1])
            elif "vel" in arg:
                vel = float(arg.split("=")[-1])


    ############################################
    #Start looping!
    ############################################
    s_spt = MS.GetSpectralType(MS.Temperature, s_temp)
    s_mass = MS.Interpolate(MS.Mass, s_spt)
    radius = MS.Interpolate(MS.Radius, s_spt)
    logg = numpy.log10(Units.G * s_mass * Units.Msun / (radius * Units.Rsun) ** 2)
    best_key = modelfiles.keys()[0]

    for key in modelfiles.keys():
        if numpy.abs(s_temp - key) < numpy.abs(s_temp - best_key):
            best_key = key
            logg_values = [float(fname.split("lte")[-1].split("-")[1]) for fname in modelfiles[key]]

    best_logg = float(modelfiles[best_key][0].split("lte")[-1].split("-")[1])
    logg_index = numpy.argmin(numpy.array(logg_values - logg))
    modelfile = modelfiles[best_key][logg_index]

    #Read in model
    print "Model file: %s" % modelfile
    x, y = numpy.loadtxt(modelfile, usecols=(0, 1), unpack=True)
    x = x * Units.nm / Units.angstrom
    y = 10 ** y
    model = DataStructures.xypoint(x=x, y=y)

    #Rotationally broaden
    model = RotBroad.Broaden(model, vsini * Units.cm / Units.km, findcont=True)

    orders = Add(list(orders_original), model, p_spt, s_spt, vel=vel * Units.cm / Units.km, SNR=1e6,
                 sensitivity=sensitivity_fcn)

    """
    for i, order in enumerate(orders[2:-1]):
      original = orders_original[i]
      original.cont = FindContinuum.Continuum(original.x, original.y, lowreject=2, highreject=2)
      std = numpy.std(order.y/order.cont - original.y/original.cont)
      shift = (std + numpy.std(original.y/original.cont))*3
      pylab.plot(original.x, original.y/original.cont, 'k-')
      pylab.plot(order.x, order.y/order.cont - original.y/original.cont + 1.0 + shift, 'r-')
    pylab.show()
    """
