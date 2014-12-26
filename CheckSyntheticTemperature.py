import os
import pandas
import re
from scipy.interpolate import InterpolatedUnivariateSpline as spline

import matplotlib.pyplot as plt
import numpy as np


def classify_filename(fname):
    """
    Given a CCF filename, it classifies the star combination, temperature, metallicity, and vsini
    :param fname:
    :return:
    """
    # First, remove any leading directories
    fname = fname.split('/')[-1]

    # Star combination
    m1 = re.search('[\.[0-9]+kps', fname)
    stars = fname[:m1.start()]
    star1 = stars.split('+')[0].replace('_', ' ')
    star2 = stars.split('+')[1].replace('_', ' ')

    # secondary star vsini
    vsini = float(fname[m1.start():].split('kps')[0])

    # Temperature
    m2 = re.search('[0-9]+\.0K', fname)
    temp = float(m2.group()[:-1])

    # logg
    m3 = re.search('K\+[0-9]\.[0-9]', fname)
    logg = float(m3.group()[1:])

    # metallicity
    metal = float(fname.split(str(logg))[-1])

    return star1, star2, vsini, temp, logg, metal


def get_ccf_data(basedir):
    """
    Searches the given directory for CCF files, and classifies
    by star, temperature, metallicity, and vsini
    :param basedir:
    :return:
    """
    if not basedir.endswith('/'):
        basedir += '/'
    all_files = ['{}{}'.format(basedir, f) for f in os.listdir(basedir)]
    primary = []
    secondary = []
    vsini_values = []
    temperature = []
    gravity = []
    metallicity = []
    ccf = []
    x = np.arange(-900.0, 900.1, 0.1)
    for fname in all_files:
        star1, star2, vsini, temp, logg, metal = classify_filename(fname)
        vel, corr = np.loadtxt(fname, unpack=True)
        fcn = spline(vel, corr)
        ccf.append(fcn(x))
        primary.append(star1)
        secondary.append(star2)
        vsini_values.append(vsini)
        temperature.append(temp)
        gravity.append(logg)
        metallicity.append(metal)

    # Make a pandas dataframe with all this data
    df = pandas.DataFrame(data={'Primary': primary, 'Secondary': secondary, 'Temperature': temperature,
                                'vsini': vsini_values, 'logg': gravity, '[Fe/H]': metallicity, 'CCF': ccf})
    return df


def temperature_plot_star(df, starname1, starname2, pixel='highest'):
    """
    Takes a pandas dataframe such as created by get_ccf_data, and the name of a primary/secondary star.
    It plots the ccf value as a function of temperature
    :param df: pandas DataFrame with the appropriate keys (use get_ccf_data)
    :param starname1: the name of the primary star
    :param starname2: the name of the secondary star
    :keyword pixel: The pixel to measure the CCF at. If 'highest', it uses the maximum of the ccf
    :return:
    """
    good = df.loc[(df.Primary == starname1) & (df.Secondary == starname2)]
    good = good.sort(columns='Temperature')
    Tvalues = good.Temperature.values
    if isinstance(pixel, str) and pixel == 'highest':
        corr_vals = [max(v) for v in good.CCF.values]
    elif isinstance(pixel, int) or isinstance(pixel, float):
        pixel = int(pixel)
        corr_vals = [v[pixel] for v in good.CCF.values]
    else:
        raise TypeError('Bad value entered for pixel ({})'.format(pixel))

    maxidx = np.argmax(corr_vals)
    maxT = Tvalues[maxidx]
    print('Maximum Temperature at {}'.format(maxT))

    plt.plot(Tvalues, corr_vals, 'r-o')
    plt.show()


def plot_temperature_distribution(df, secondary, pixel='highest'):
    """
    Goes through all entries with the given secondary star, finds the best temperature for each one, and plots
    the distribution of temperatures (one for each primary star). Returns the distribution as well.
    :param df: pandas DataFrame with the appropriate keys (use get_ccf_data)
    :param secondary: The name of the secondary star
    :param pixel: The pixel to measure the CCF at. If 'highest', it uses the maximum of the ccf
    :return: numpy.ndarray with the best temperature for each primary star
    """
    good = df.loc[df.Secondary == secondary]
    primary_names = pandas.unique(good.Primary)
    bestT = []
    for primary in primary_names:
        print primary
        prim_df = good.loc[good.Primary == primary]
        Tvalues = prim_df.Temperature.values
        if isinstance(pixel, str) and pixel == 'highest':
            corr_vals = [max(v) for v in prim_df.CCF.values]
        elif isinstance(pixel, int) or isinstance(pixel, float):
            pixel = int(pixel)
            corr_vals = [v[pixel] for v in prim_df.CCF.values]
        else:
            raise TypeError('Bad value entered for pixel ({})'.format(pixel))
        maxidx = np.argmax(corr_vals)
        bestT.append(Tvalues[maxidx])

    plt.hist(bestT)
    plt.show()

    return bestT


if __name__ == '__main__':
    df = get_ccf_data('Cross_correlations/')
    temperature_plot_star(df, 'HIP 42313', 'Gam Ser')