import sys
from collections import defaultdict
from operator import itemgetter

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import SpectralTypeRelations
import PlotBlackbodies
import Units


"""
  Program to analyze the output of SensitivityAnalysis, and make some pretty plots!
  Very similar to MakeSensitivityPlot.py, but takes output files with S/N in them (2nd column)
     Makes separate figures for the different S/N values
  Command line arguments:
     -combine: will combine several output (say as generated by xgrid) NOT YET IMPLEMENTED
     -xaxis: specifies the variable to use as the x axis. Choices are as follows
         SecondarySpectralType
         SecondaryMass
         MassRatio
         DetectionRate
         AverageSignificance
         MagnitudeDifference
     -yaxis: specifies the variable to use for the y axis. Choices are the same as for -xaxis
     -infile: specifies the input filename (default is Sensitivity/summary.dat).
         If combine is True, the input filename should be a list of comma-separated 
         filenames
     -oneplot: Uses only the highest S/N observation of a given parent spectral type, and
         makes a plot similar to MakeSensitivityPlot.py
     -monochrome: Makes plot black and white
"""


# Set up thing to cycle through matplotlib linestyles
from itertools import cycle

lines = ["-", "--", "-.", ":"]
linecycler = cycle(lines)

if __name__ == "__main__":
    #Defaults
    combine = False
    xaxis = "SecondarySpectralType"
    yaxis = "DetectionRate"
    infilename = "Sensitivity/summary.dat"
    oneplot = False
    #palette = cm.hsv

    #Command-line overrides
    for arg in sys.argv:
        if "combine" in arg:
            combine = True
        elif "xaxis" in arg:
            xaxis = arg.split("=")[-1]
        elif "yaxis" in arg:
            yaxis = arg.split("=")[-1]
        elif "infile" in arg:
            infilename = arg.split("=")[-1]
        elif "oneplot" in arg:
            oneplot = True
            #elif "mono" in arg:
            #palette = cm.Greys
    if combine and "," in infilename:
        infiles = infilename.split(",")
    else:
        infiles = [infilename, ]

    #Set up dictionaries/lists
    s_spt = defaultdict(lambda: defaultdict(list))  #Secondary spectral type
    p_mass = defaultdict(lambda: defaultdict(list))  #Primary mass
    s_mass = defaultdict(lambda: defaultdict(list))  #Secondary mass
    q = defaultdict(lambda: defaultdict(list))  #Mass ratio
    det_rate = defaultdict(lambda: defaultdict(list))  #Detection rate
    sig = defaultdict(lambda: defaultdict(list))  #Average detection significance
    magdiff = defaultdict(lambda: defaultdict(list))  #Magnitude difference
    namedict = {"SecondarySpectralType": s_spt,
                "SecondaryMass": s_mass,
                "MassRatio": q,
                "DetectionRate": det_rate,
                "AverageSignificance": sig,
                "MagnitudeDifference": magdiff}
    labeldict = {"SecondarySpectralType": "Secondary Spectral Type",
                 "SecondaryMass": "SecondaryMass (Solar Masses)",
                 "MassRatio": "Mass Ratio",
                 "DetectionRate": "Detection Rate (Percent)",
                 "AverageSignificance": "Average Significance",
                 "MagnitudeDifference": "Magnitude Difference"}

    if xaxis not in namedict.keys() or yaxis not in namedict:
        print "Error! axis keywords must be one of the following:"
        for key in namedict.keys():
            print key
        sys.exit()

    MS = SpectralTypeRelations.MainSequence()
    vband = np.arange(500, 600, 0.1) * Units.cm / Units.nm
    #Read in file/files  WARNING! ASSUMES A CERTAIN FORMAT. MUST CHANGE THIS IF THE FORMAT CHANGES!

    for infilename in infiles:
        starname = infilename.split("_")[-1].split(".")[0]
        infile = open(infilename)
        lines = infile.readlines()
        for line in lines[1:]:
            segments = line.split()
            p_spt = segments[0]
            snr = float(segments[1])
            T1 = MS.Interpolate(MS.Temperature, p_spt)
            R1 = MS.Interpolate(MS.Radius, p_spt)
            T2 = MS.Interpolate(MS.Temperature, segments[2])
            R2 = MS.Interpolate(MS.Radius, segments[2])
            fluxratio = (PlotBlackbodies.Planck(vband, T1) / PlotBlackbodies.Planck(vband, T2)).mean() * (R1 / R2) ** 2
            s_spt[p_spt][snr].append(segments[2])
            p_mass[p_spt][snr].append(float(segments[3]))
            s_mass[p_spt][snr].append(float(segments[4]))
            q[p_spt][snr].append(float(segments[5]))
            det_rate[p_spt][snr].append(float(segments[6]))
            sig[p_spt][snr].append(float(segments[7]))
            magdiff[p_spt][snr].append(2.5 * np.log10(fluxratio))


        #plot
        index = 0
        spt_sorter = {"O": 1, "B": 2, "A": 3, "F": 4, "G": 5, "K": 6, "M": 7}
        fcn = lambda s: (spt_sorter[itemgetter(0)(s)], itemgetter(1)(s))
        for p_spt in sorted(s_spt.keys(), key=fcn):
            if oneplot:
                #Find highest s/n in this spectral type
                highestsnr = sorted(s_spt[p_spt].keys())[-1]
                x = namedict[xaxis][p_spt][highestsnr]
                y = namedict[yaxis][p_spt][highestsnr]
                if "SpectralType" in xaxis:
                    plt.plot(range(len(x)), y, linestyle=next(linecycler), linewidth=2,
                             label="%s (%s)" % (starname, p_spt))
                    plt.xticks(range(len(x)), x, size="small")
                elif "SpectralType" in yaxis:
                    plt.plot(x, range(len(y)), linestyle=next(linecycler), linewidth=2,
                             label="%s (%s)" % (starname, p_spt))
                    plt.yticks(range(len(y)), y, size="small")
                else:
                    plt.plot(x, y, linestyle=next(linecycler), linewidth=2, label="%s (%s)" % (starname, p_spt))


            else:
                plt.figure(index)
                index += 1
                for snr in sorted(s_spt[p_spt].keys()):
                    x = namedict[xaxis][p_spt][snr]
                    y = namedict[yaxis][p_spt][snr]
                    if "SpectralType" in xaxis:
                        plt.plot(range(len(x)), y, linestyle=next(linecycler), linewidth=2, label="S/N = %.0f" % snr)
                        plt.xticks(range(len(x)), x, size="small")
                    elif "SpectralType" in yaxis:
                        plt.plot(x, range(len(y)), linestyle=next(linecycler), linewidth=2, label="S/N = %.0f" % snr)
                        plt.yticks(range(len(y)), y, size="small")
                    else:
                        plt.plot(x, y, linestyle=next(linecycler), linewidth=2, label="S/N = %.0f" % snr)
                        if "Magnitude" in xaxis:
                            ax = plt.gca()
                            ax.set_xlim(ax.get_xlim()[::-1])
                        elif "Magnitude" in yaxis:
                            ax = plt.gca()
                            ax.set_ylim(ax.get_ylim()[::-1])
                plt.legend(loc='best')
                plt.xlabel(labeldict[xaxis])
                plt.ylabel(labeldict[yaxis])
                plt.title("Sensitivity for " + p_spt + " Primary")
        if oneplot:
            plt.legend(loc='best')
            plt.xlabel(labeldict[xaxis])
            plt.ylabel(labeldict[yaxis])
            if "Magnitude" in xaxis:
                ax = plt.gca()
                ax.set_xlim(ax.get_xlim()[::-1])
            elif "Magnitude" in yaxis:
                ax = plt.gca()
                ax.set_ylim(ax.get_ylim()[::-1])
    plt.show()



