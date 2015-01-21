"""
This script does a chi-squared search to find the best model fit to the given data
"""

import sys
import FittingUtilities

import numpy as np
import StellarModel
import HelperFunctions
from scipy.interpolate import InterpolatedUnivariateSpline as spline
import matplotlib.pyplot as plt
import Broaden
from collections import defaultdict
import pandas as pd


if "darwin" in sys.platform:
    modeldir = "/Volumes/DATADRIVE/Stellar_Models/Sorted/Stellar/Vband/"
elif "linux" in sys.platform:
    modeldir = "/media/FreeAgent_Drive/SyntheticSpectra/Sorted/Stellar/Vband/"
else:
    modeldir = raw_input("sys.platform not recognized. Please enter model directory below: ")
    if not modeldir.endswith("/"):
        modeldir = modeldir + "/"


def residual(orders, model, plot=False):
    fcn = spline(model.x, model.y)
    resid = []
    for order in orders:
        m = fcn(order.x)
        c = FittingUtilities.Continuum(order.x, m, fitorder=3, lowreject=1.5, highreject=10)
        resid.append(order.y / order.cont - m / c)
        if plot:
            plt.plot(order.x, order.y / order.cont, 'k-', alpha=0.4)
            plt.plot(order.x, m / c, 'r-', alpha=0.5)
    if plot:
        plt.xlabel('Wavelength')
        plt.ylabel('Normalized Flux')
        plt.xlim((600, 700))
        plt.ylim((-0.5, 2))
        plt.show()
    return np.hstack(resid)


if __name__ == '__main__':
    file_list = []
    rv = None
    c = 3e5
    vsini = None
    plotflg=False
    for arg in sys.argv[1:]:
        if '--rv' in arg:
            rv_vals = arg.split('=')[-1].split(',')
            rv = [float(v) for v in rv_vals]
            #rv = float(arg.split('=')[-1])
        elif '--vsini' in arg:
            vsini = float(arg.split('=')[-1])
        else:
            file_list.append(arg)

    if len(file_list) == 0:
        sys.exit('Must give at least one file to fit as a command-line argument!')

    if rv is None:
        rv = [0.0]*len(file_list)
    elif len(rv) != len(file_list):
        sys.exit('Must give the same number of rv values as files!')

    # Read in the stellar models
    model_list = StellarModel.GetModelList(metal=[0], alpha=[0], model_directory=modeldir,
                                           temperature=range(5000, 5200, 100))
    model_dict = StellarModel.MakeModelDicts(model_list, vsini_values=[1], logspace=True)[0]

    summary_dict = defaultdict(list)
    for n, fname in enumerate(file_list):
        print 'FILE {}'.format(fname)
        orders = HelperFunctions.ReadExtensionFits(fname)

        # Re-measure continuum
        for i, order in enumerate(orders):
            orders[i].cont = FittingUtilities.Continuum(order.x, order.y, fitorder=3, lowreject=2, highreject=5)

        # Loop through the models
        chi2 = []
        Tvalues = []
        rv_values = []
        for T in sorted(model_dict.keys()):
            for vel in np.arange(rv[n]-1.0, rv[n]+1.1, 0.1):
                model = model_dict[T][4.50][0.0][1.0].copy()
                if vsini is not None:
                    model = Broaden.RotBroad(model, vsini * 1e5)
                print ('RV = {} km/s'.format(vel))
                model.x *= (1.0 + vel / c)
                print 'T = {}'.format(T)
                chi2.append(np.sum(residual(orders, model, plot=plotflg) ** 2))
                Tvalues.append(float(T))
                print '{}  {}'.format(chi2[-1], Tvalues[-1])

        # Find the best T
        print chi2
        print Tvalues
        print rv_values
        idx = np.argmin(chi2)
        print idx
        T = Tvalues[idx]
        vel = rv_values[idx]
        print 'Best T = {} K'.format(T)
        print 'Best rv = {} km/s'.format(vel)

        if plotflg:
            plt.plot(Tvalues, chi2, 'ro')
            plt.xlabel('Temperature')
            plt.ylabel(r'$\chi^2$')
            plt.show()

        summary_dict['file'].append(fname)
        summary_dict['T'].append(T)
        summary_dict['rv'].append(vel)

    summary = pd.DataFrame(data=summary_dict)
    summary.to_csv('Chisquared_Summary.csv', index=False)

