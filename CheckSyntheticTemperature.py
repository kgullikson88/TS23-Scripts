import os
import pandas
from scipy.interpolate import InterpolatedUnivariateSpline as spline
import re

import matplotlib.pyplot as plt
import numpy as np

import StarData
import SpectralTypeRelations


def classify_filename(fname):
    """
    Given a CCF filename, it classifies the star combination, temperature, metallicity, and vsini
    :param fname:
    :return:
    """
    # First, remove any leading directories
    fname = fname.split('/')[-1]

    # Star combination
    m1 = re.search('\.[0-9]+kps', fname)
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


def get_ccf_data(basedir, primary_name=None, secondary_name=None, vel_arr=np.arange(-900.0, 900.0, 0.1), type=None):
    """
    Searches the given directory for CCF files, and classifies
    by star, temperature, metallicity, and vsini
    :param basedir: The directory to search for CCF files
    :keyword primary_name: Optional keyword. If given, it will only get the requested primary star data
    :keyword secondary_name: Same as primary_name, but only reads ccfs for the given secondary
    :keyword vel_arr: The velocities to interpolate each ccf at
    :return: pandas DataFrame
    """
    if not basedir.endswith('/'):
        basedir += '/'
    if type is None:
        all_files = ['{}{}'.format(basedir, f) for f in os.listdir(basedir) if 'MS_scale' not in f]
    elif 'ms' in type.lower():
        all_files = ['{}{}'.format(basedir, f) for f in os.listdir(basedir) if 'MS_scale' in f]
    primary = []
    secondary = []
    vsini_values = []
    temperature = []
    gravity = []
    metallicity = []
    ccf = []
    for fname in all_files:
        star1, star2, vsini, temp, logg, metal = classify_filename(fname)
        if temp == 5000:
            print star1, star2, secondary_name
        if primary_name is not None and star1.lower() != primary_name.lower():
            continue
        if secondary_name is not None and star2.lower() != secondary_name.lower():
            continue
        vel, corr = np.loadtxt(fname, unpack=True)
        fcn = spline(vel, corr)
        ccf.append(fcn(vel_arr))
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


def temperature_plot_star(df, starname1, starname2, velocity='highest', vel_arr=np.arange(-900.0, 900.0, 0.1)):
    """
    Takes a pandas dataframe such as created by get_ccf_data, and the name of a primary/secondary star.
    It plots the ccf value as a function of temperature
    :param df: pandas DataFrame with the appropriate keys (use get_ccf_data)
    :param starname1: the name of the primary star
    :param starname2: the name of the secondary star
    :keyword velocity: The velocity to measure the CCF at. The default is 'highest', and uses the maximum of the ccf
    :keyword vel_arr: The velocities to interpolate each ccf at
    :return:
    """
    good = df.loc[(df.Primary == starname1) & (df.Secondary == starname2)]
    good = good.sort(columns='Temperature')
    Tvalues = good.Temperature.values
    if isinstance(velocity, str) and velocity == 'highest':
        corr_vals = [max(v) for v in good.CCF.values]
    elif isinstance(velocity, int) or isinstance(velocity, float):
        pixel = int(velocity)
        corr_vals = [v[pixel] for v in good.CCF.values]
    else:
        raise TypeError('Bad value entered for pixel ({})'.format(velocity))

    maxidx = np.argmax(corr_vals)
    maxT = Tvalues[maxidx]
    print('Maximum Temperature at {}'.format(maxT))

    plt.plot(Tvalues, corr_vals, 'r-o')
    plt.show()


def get_best_temperature(df, velocity='highest', vel_arr=np.arange(-900.0, 900.0, 0.1)):
    """
    Given a DataFrame with just primary and secondary star combination, find the best temperature
    :param df: DataFrame to search
    :keyword velocity: The velocity to measure the CCF at. The default is 'highest', and uses the maximum of the ccf
    :keyword vel_arr: The velocities to interpolate each ccf at
    :return: best temperature (as a float object)
    """
    Tvalues = df.Temperature.values
    if isinstance(velocity, str) and velocity == 'highest':
        corr_vals = [max(v) for v in df.CCF.values]
    elif isinstance(velocity, int) or isinstance(velocity, float):
        pixel = np.argmin(abs(vel_arr - velocity))
        corr_vals = [v[pixel] for v in df.CCF.values]
    else:
        raise TypeError('Bad value entered for velocity ({})'.format(velocity))
    maxidx = np.argmax(corr_vals)
    return Tvalues[maxidx]


def plot_temperature_distribution(df, secondary, velocity='highest', vel_arr=np.arange(-900.0, 900.0, 0.1)):
    """
    Goes through all entries with the given secondary star, finds the best temperature for each one, and plots
    the distribution of temperatures (one for each primary star). Returns the distribution as well.
    :param df: pandas DataFrame with the appropriate keys (use get_ccf_data)
    :param secondary: The name of the secondary star
    :param velocity: The velocity to measure the CCF at. If 'highest', it uses the maximum of the ccf
    :keyword vel_arr: The velocities to interpolate each ccf at
    :return: numpy.ndarray with the best temperature for each primary star
    """
    good = df.loc[df.Secondary.str.upper() == secondary.upper()]
    primary_names = pandas.unique(good.Primary)
    bestT = []
    for primary in primary_names:
        print primary
        prim_df = good.loc[good.Primary == primary]
        bestT.append(get_best_temperature(prim_df, vel_arr=vel_arr, velocity=velocity))

    plt.hist(bestT)
    plt.show()

    return bestT


def plot_temperature_accuracy(df, velocity='highest', vel_arr=np.arange(-900.0, 900.0, 0.1)):
    """
     Plots the best temperature for each combination of primary and secondary
    :param df:
    :param velocity: The velocity to measure the CCF at. If 'highest', it uses the maximum of the ccf
    :keyword vel_arr: The velocities to interpolate each ccf at
    :return:
    """
    primary_names = pandas.unique(df.Primary)
    secondary_names = pandas.unique(df.Secondary)

    MS = SpectralTypeRelations.MainSequence()

    Tactual = []  # Temperature, as obtained from the spectral type of the secondary star
    Tmeas = []
    for secondary in secondary_names:
        star_data = StarData.GetData(secondary)
        spt = star_data.spectype[0] + re.search('[0-9]\.*[0-9]*', star_data.spectype).group()
        T_sec = MS.Interpolate(MS.Temperature, spt)
        for primary in primary_names:
            good = df.loc[(df.Primary == primary) & (df.Secondary == secondary)]
            Tmeas.append(get_best_temperature(good, velocity=velocity, vel_arr=vel_arr))
            Tactual.append(T_sec)

    plt.scatter(Tactual, Tmeas)
    plt.plot(Tactual, Tactual, 'r--')
    plt.savefig('Temperature_Accuracy.svg')
    # plt.show()

    return Tactual, Tmeas


if __name__ == '__main__':
    secondary = '16 Cyg B'
    df = get_ccf_data('Cross_correlations/')
    # temperature_plot_star(df, 'HIP 42313', 'Gam Ser')
    #plot_temperature_distribution(df, secondary)
    plot_temperature_accuracy(df)