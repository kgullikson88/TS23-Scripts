"""
This is a script to measure the equivalent width of interstellar sodium.
The main function is Measure(filename), which returns the equivalent
widths of both sodium D lines, as well as the reddening E(B-V).
The relationship from EW to E(B-V) is taken from 
Poznanski et al. 2012 MNRAS, 426, 1465
"""

import sys
import os
import FittingUtilities

from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
import numpy as np
import DataStructures
from astropy.io import fits
import HelperFunctions
import pyspeckit



# ordernum = 30  #The order number with the Sodium lines in it
D1_wave = 589.592
D2_wave = 588.995

#Make matplotlib interactive
plt.ion()


def FindOrderNum(orders, wavelength):
    """
      Given a list of xypoint orders and
      another list of wavelengths, this
      finds the order numbers with the
      requested wavelengths
    """
    num = 0
    for i, order in enumerate(orders):
        if order.x[0] < wavelength and order.x[-1] > wavelength:
            num = i
            break
    return num


def GetReddening(D1_ew, D2_ew):
    #Use the sum of the lines, since it has the lowest error
    EB_V = 10 ** (1.17 * (D1_ew + D2_ew) - 1.85)
    return EB_V


def Measure3(filename, figure, title):
    orders = HelperFunctions.ReadExtensionFits(filename)
    ordernum = FindOrderNum(orders, 589)
    order = orders[ordernum]
    left = np.searchsorted(order.x, 588)
    right = np.searchsorted(order.x, 590.5)
    order = order[left:right]  # Get the line positions
    order.y /= order.cont
    order.err /= order.cont
    xgrid = np.linspace(order.x[0], order.x[-1], order.size())
    order = FittingUtilities.RebinData(order, xgrid)
    lines = FittingUtilities.FindLines(order, tol=0.97, linespacing=0.5)
    print lines
    D1 = 9e9
    D2 = 9e9
    for line in lines:
        if abs(order.x[line] - D1_wave) < abs(D1 - D1_wave):
            D1 = order.x[line]
        if abs(order.x[line] - D2_wave) < abs(D2 - D2_wave):
            D2 = order.x[line]

    #Ask user for where the Na lines start
    first, second = D2 - 0.1, D2 + 0.1
    third, fourth = D1 - 0.1, D1 + 0.1
    ax = figure.add_subplot(111)
    ax.plot(order.x, order.y)
    yrange = ax.get_ylim()

    #Get the range for the Na line
    ax.plot((first, first), yrange, 'r-')
    ax.plot((second, second), yrange, 'r-')
    ax.plot((third, third), yrange, 'r-')
    ax.plot((fourth, fourth), yrange, 'r-')
    plt.show()
    inp = raw_input("Give bounds for Na line (hit enter if okay) ")
    if inp.strip() == "":
        first, second = D2 - 0.1, D2 + 0.1
        third, fourth = D1 - 0.1, D1 + 0.1
    else:
        segments = inp.split()
        first = float(segments[0])
        second = float(segments[1])
        third = float(segments[2])
        fourth = float(segments[3])

    #Make a new array, not including the Na line
    left1 = np.searchsorted(order.x, first)
    right1 = np.searchsorted(order.x, second)
    left2 = np.searchsorted(order.x, third)
    right2 = np.searchsorted(order.x, fourth)
    x = np.r_[order.x[:left1], order.x[right1:left2], order.x[right2:]]
    y = np.r_[order.y[:left1], order.y[right1:left2], order.y[right2:]]
    e = np.r_[order.err[:left1], order.err[right1:left2], order.err[right2:]]
    order2 = DataStructures.xypoint(x=x, y=y, err=e)
    order2.output("Spec.dat")

    done = False
    smooth_value = 3e-4
    while not done:
        done = True
        #Smooth the new array, and use the smoothed value as continuum
        fcn = UnivariateSpline(order2.x, order2.y, s=3e-4)
        order.cont = fcn(order.x)

        #Plot this value to make sure it looks fine
        ax = figure.add_subplot(111)
        plt.cla()
        ax.plot(order.x, order.y)
        ax.plot(order.x, order.cont)
        plt.show()
        inp = raw_input("Continuum fit okay? ")
        if inp.strip() == "":
            done = True
        elif "+" in inp:
            smooth_value += 1e-4
        elif "-" in inp:
            smooth_value -= 1e-4

    #Measure the equivalent width
    dx = np.array([order.x[i + 1] - order.x[i] for i in range(left1, right1)])
    ew1 = np.sum((1.0 - order.y[left1:right1] / order.cont[left1:right1]) * dx)
    dx = np.array([order.x[i + 1] - order.x[i] for i in range(left2, right2)])
    ew2 = np.sum((1.0 - order.y[left2:right2] / order.cont[left2:right2]) * dx)

    return [ew1, ew2], GetReddening(ew1, ew2)


def Measure2(filename, figure, title):
    """
    Given a filename, this will plot the appropriate order
    and help the user measure the equivalent width of
    both sodium D lines. Unlike the 'Measure' function,
    this one just performs a simpson rule interpolation
    to integrate under the curve
    """
    orders = HelperFunctions.ReadExtensionFits(filename)
    ordernum = FindOrderNum(orders, 589)
    order = orders[ordernum]
    left = np.searchsorted(order.x, 588)
    right = np.searchsorted(order.x, 590.5)
    order = order[left:right]  # Get the Na line positions
    order.y /= order.cont
    xgrid = np.linspace(order.x[0], order.x[-1], order.size())
    order = FittingUtilities.RebinData(order, xgrid)
    lines = FittingUtilities.FindLines(order, tol=0.97, linespacing=0.5)
    print lines
    D1 = 9e9
    D2 = 9e9
    for line in lines:
        if abs(order.x[line] - D1_wave) < abs(D1 - D1_wave):
            D1 = order.x[line]
        if abs(order.x[line] - D2_wave) < abs(D2 - D2_wave):
            D2 = order.x[line]

    # We need to fit any large lines as continuum
    # Make a new array, not including the Na line
    first, second = D2 - 0.1, D2 + 0.1
    third, fourth = D1 - 0.1, D1 + 0.1
    left1 = np.searchsorted(order.x, first)
    right1 = np.searchsorted(order.x, second)
    left2 = np.searchsorted(order.x, third)
    right2 = np.searchsorted(order.x, fourth)
    x = np.r_[order.x[:left1], order.x[right1:left2], order.x[right2:]]
    y = np.r_[order.y[:left1], order.y[right1:left2], order.y[right2:]]
    e = np.r_[order.err[:left1], order.err[right1:left2], order.err[right2:]]
    order2 = DataStructures.xypoint(x=x, y=y, err=e)

    # Let user decide if the fit is okay
    done = False
    smooth_value = 3e-4
    while not done:
        done = True
        #Smooth the new array, and use the smoothed value as continuum
        fcn = UnivariateSpline(order2.x, order2.y, s=smooth_value)
        order.cont = fcn(order.x)

        # Plot this value to make sure it looks fine
        ax = figure.add_subplot(111)
        plt.cla()
        ax.plot(order.x, order.y)
        ax.plot(order.x, order.cont)
        plt.show()
        inp = raw_input("Continuum fit okay? (Hit '+', '-', or enter if okay)")
        if inp.strip() == "":
            done = True
        elif "+" in inp:
            smooth_value += 1e-4
            done = False
        elif "-" in inp:
            if smooth_value - 1e-4 <= 0:
                smooth_value /= 2.0
            else:
                smooth_value -= 1e-4
            done = False

    x = pyspeckit.units.SpectroscopicAxis(order.x, units='nm')
    spec = pyspeckit.Spectrum(data=order.y / order.cont, xarr=x, error=order.err / order.cont, units='nm')
    pl = spec.plotter(autorefresh=True, figure=figure)
    plt.title(title)

    # Fit the baseline in this region (continuum)
    spec.baseline(subtract=False, order=10, interactive=True)
    plt.show()
    done = raw_input("hit enter ")

    start1, end1 = spec.baseline.xclicks[1], spec.baseline.xclicks[2]
    start2, end2 = spec.baseline.xclicks[3], spec.baseline.xclicks[4]
    start1 = spec.xarr[start1]
    end1 = spec.xarr[end1]
    start2 = spec.xarr[start2]
    end2 = spec.xarr[end2]
    print spec.baseline.xclicks
    print start1, end1
    print start2, end2

    # Find the equivalent width for both lines
    ews = []
    for line in [[start1, end1], [start2, end2]]:
        xmin = np.searchsorted(order.x, line[0])
        xmax = np.searchsorted(order.x, line[1])
        ew = spec.specfit.EQW(xmin=xmin, xmax=xmax, fitted=False, plot=True)
        ews.append(ew * 10)  #Save the equivalent width in angstroms
    done = raw_input("hit enter ")
    print "\n\n"

    return ews, GetReddening(*ews)


def Measure(filename, figure, title):
    """
    Given a filename, this will plot the appropriate order
    and help the user measure the equivalent width of
    both sodium D lines
    """
    orders = HelperFunctions.ReadExtensionFits(filename)
    ordernum = FindOrderNum(orders, 589)
    order = orders[ordernum]
    left = np.searchsorted(order.x, 588)
    right = np.searchsorted(order.x, 590.5)
    order = order[left:right]

    x = pyspeckit.units.SpectroscopicAxis(order.x, units='nm')
    spec = pyspeckit.Spectrum(data=order.y / order.cont, xarr=x, error=order.err / order.cont, units='nm')
    pl = spec.plotter(autorefresh=True, figure=figure)
    plt.title(title)

    # Fit the baseline in this region (continuum)
    spec.baseline(subtract=False, order=3, interactive=True)
    plt.show()
    done = raw_input("hit enter ")

    # Fit voigt profiles to the lines
    fitguesses = [-0.5, D1_wave, 0.08, 0.0,
                  -0.5, D2_wave, 0.08, 0.0]
    spec.specfit(fittype='voigt', multifit=True, guesses=fitguesses)
    plt.draw()

    # Grab the fitted parameters and their errors
    fitpars = spec.specfit.modelpars
    fiterrs = spec.specfit.modelerrs

    # Sometimes, the lines 'switch' places, so check that
    if fitpars[1] > fitpars[5]:
        fitpars = fitpars[4:] + fitpars[:4]
        fiterrs = fiterrs[4:] + fiterrs[:4]

    #Determine the start and end of each line
    start1 = fitpars[1] - 7 * np.sqrt(fitpars[2] ** 2 + fitpars[3] ** 2)
    end1 = fitpars[1] + 7 * np.sqrt(fitpars[2] ** 2 + fitpars[3] ** 2)
    start2 = fitpars[5] - 7 * np.sqrt(fitpars[6] ** 2 + fitpars[7] ** 2)
    end2 = fitpars[5] + 7 * np.sqrt(fitpars[6] ** 2 + fitpars[7] ** 2)

    #Find the equivalent width for both lines
    ews = []
    for line in [[start1, end1], [start2, end2]]:
        xmin = np.searchsorted(order.x, line[0])
        xmax = np.searchsorted(order.x, line[1])
        ew = spec.specfit.EQW(xmin=xmin, xmax=xmax, fitted=True, continuum=1.0, plot=True)
        ews.append(ew * 10)  #Save the equivalent width in angstroms
    done = raw_input("hit enter ")
    print "\n\n"

    return ews, GetReddening(*ews)


if __name__ == "__main__":
    # Read in the lines in the current log file
    logfilename = "%s/Dropbox/School/Research/AstarStuff/TargetLists/Na_Params.csv" % os.environ[
        'HOME']  #Parse command line arguments
    fileList = []
    for arg in sys.argv[1:]:
        fileList.append(arg)

# For each file, determine the EW and E(B-V)
fig = plt.figure()
for fname in fileList:
    #First, check to make sure that this target was not already measured
    logfile = open(logfilename, "r")
    lines = logfile.readlines()
    logfile.close()
    measure = True
    overwrite = 9999999
    starname = fits.getheader(fname)["object"]
    print starname
    #Check if this star already has an entry in the logfile
    for linenum, line in enumerate(lines):
        star = line.split("|")[0].strip()
        if star.lower() == starname.lower():
            measure = False
            print "%s found in the master logfile." % starname
            inp = raw_input("Do you want to overwrite (y/n)? ")
            if inp.lower() == "y":
                measure = True
                overwrite = linenum

            if measure:
                ews, reddening = Measure2(fname, fig, starname)
                print ews, reddening

                outline = "%s | %.4f | %.4f | %.4f\n" % (starname, ews[0], ews[1], reddening)
                if overwrite > len(lines):
                    lines.append(outline)
                else:
                    lines[overwrite] = outline

        # Finally, output the new logfile
        logfile = open(logfilename, "w")
        logfile.writelines(lines)
        logfile.close()