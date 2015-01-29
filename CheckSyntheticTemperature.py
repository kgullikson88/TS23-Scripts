import os
import re
import sys
from george import kernels
from collections import defaultdict

import pandas
from scipy.interpolate import InterpolatedUnivariateSpline as spline
import matplotlib.pyplot as plt
import numpy as np
import StarData
import SpectralTypeRelations
import george
import emcee


def classify_filename(fname, type='bright'):
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
    star2 = stars.split('+')[1].split('_{}'.format(type))[0].replace('_', ' ')

    # secondary star vsini
    vsini = float(fname[m1.start() + 1:].split('kps')[0])

    # Temperature
    m2 = re.search('[0-9]+\.0K', fname)
    temp = float(m2.group()[:-1])

    # logg
    m3 = re.search('K\+[0-9]\.[0-9]', fname)
    logg = float(m3.group()[1:])

    # metallicity
    metal = float(fname.split(str(logg))[-1])

    return star1, star2, vsini, temp, logg, metal


def get_ccf_data(basedir, primary_name=None, secondary_name=None, vel_arr=np.arange(-900.0, 900.0, 0.1), type='bright'):
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
    all_files = ['{}{}'.format(basedir, f) for f in os.listdir(basedir) if type in f.lower()]
    primary = []
    secondary = []
    vsini_values = []
    temperature = []
    gravity = []
    metallicity = []
    ccf = []
    for fname in all_files:
        star1, star2, vsini, temp, logg, metal = classify_filename(fname, type=type)
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


def find_best_pars(df, velocity='highest', vel_arr=np.arange(-900.0, 900.0, 0.1)):
    """
    Find the 'best-fit' parameters for each combination of primary and secondary star
    :param df: the dataframe to search in
    :keyword velocity: The velocity to measure the CCF at. The default is 'highest', and uses the maximum of the ccf
    :keyword vel_arr: The velocities to interpolate each ccf at
    :return: a dataframe with keys of primary, secondary, and the parameters
    """
    # Get the names of the primary and secondary stars
    primary_names = pandas.unique(df.Primary)
    secondary_names = pandas.unique(df.Secondary)

    # Find the ccf value at the given velocity
    if velocity == 'highest':
        df['ccf_max'] = df['CCF'].map(np.max)
    else:
        df['ccf_max'] = df['CCF'].map(lambda arr: arr[np.argmin(np.abs(vel_arr - velocity))])

    # Find the best parameter for each combination
    d = defaultdict(list)
    for primary in primary_names:
        for secondary in secondary_names:
            good = df.loc[(df.Primary == primary) & (df.Secondary == secondary)]
            best = good.loc[good.ccf_max == good.ccf_max.max()]
            d['Primary'].append(primary)
            d['Secondary'].append(secondary)
            d['Temperature'].append(best['Temperature'].item())
            d['vsini'].append(best['vsini'].item())
            d['logg'].append(best['logg'].item())
            d['[Fe/H]'].append(best['[Fe/H]'].item())

    return pandas.DataFrame(data=d)

def add_actual_temperature(df, method='spt'):
    """
    Add the actual temperature to a given summary dataframe
    :param df: The dataframe to which we will add the actual secondary star temperature
    :param method: How to get the actual temperature. Options are:
                   - 'spt': Use main-sequence relationships to go from spectral type --> temperature
                   - 'excel': Use tabulated data, available in the file 'SecondaryStar_Temperatures.xls'
    :return: copy of the original dataframe, with an extra column for the secondary star temperature
    """
    # First, get a list of the secondary stars in the data
    secondary_names = pandas.unique(df.Secondary)
    secondary_to_temperature = defaultdict(float)

    if method.lower() == 'spt':
        MS = SpectralTypeRelations.MainSequence()
        for secondary in secondary_names:
            star_data = StarData.GetData(secondary)
            spt = star_data.spectype[0] + re.search('[0-9]\.*[0-9]*', star_data.spectype).group()
            T_sec = MS.Interpolate(MS.Temperature, spt)
            secondary_to_temperature[secondary] = T_sec

    elif method.lower() == 'excel':
        table = pandas.read_excel('SecondaryStar_Temperatures.xls')
        for secondary in secondary_names:
            T_sec = table.loc[table.Star.str.lower().str.contains(secondary.strip().lower())]['Literature_Temp'].item()
            secondary_to_temperature[secondary] = T_sec

    df['Tactual'] = df['Secondary'].map(lambda s: secondary_to_temperature[s])
    return

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
    truth = []
    medians = []
    low = []
    high = []
    N = max([len(s) for s in secondary_names])
    for i, secondary in enumerate(secondary_names):
        sys.stdout.write('\r' + ' '*(N+34))
        sys.stdout.flush()
        sys.stdout.write("\rCompiling stats for {}, star {}/{}".format(secondary, i+1, len(secondary_names)))
        sys.stdout.flush()
        star_data = StarData.GetData(secondary)
        spt = star_data.spectype[0] + re.search('[0-9]\.*[0-9]*', star_data.spectype).group()
        T_sec = MS.Interpolate(MS.Temperature, spt)
        measured_values = []
        for primary in primary_names:
            good = df.loc[(df.Primary == primary) & (df.Secondary == secondary)]
            bestT = get_best_temperature(good, velocity=velocity, vel_arr=vel_arr)
            Tmeas.append(bestT)
            Tactual.append(T_sec)
            measured_values.append(bestT)

        # Compile statistics for this secondary
        medians.append(np.percentile(measured_values, 50.0))
        low.append(np.percentile(measured_values, 5.0))
        high.append(np.percentile(measured_values, 95.0))
        truth.append(T_sec)

    # Sort the stats so we get a decent plot
    sorter = np.argsort(truth)
    medians = np.array(medians)[sorter]
    low = np.array(low)[sorter]
    high = np.array(high)[sorter]
    truth = np.array(truth)[sorter]

    # plt.scatter(Tactual, Tmeas)
    # Plot the statistics
    plt.plot(truth, medians, 'k-', lw=2, label='Median Measured Temperature')
    plt.fill_between(truth, high, low, facecolor='green', alpha=0.4)
    plt.plot([], [], color='green', alpha=0.4, label=r'$2\sigma$ region', lw=10)  #Dummy region for the legend
    plt.plot(Tactual, Tactual, 'r--', label='Actual Temperature')
    leg = plt.legend(loc='best', fancybox=True)
    plt.xlabel('Actual Temperature', fontsize=15)
    plt.ylabel('Measured Temperature', fontsize=15)
    leg.get_frame().set_alpha(0.4)
    plt.savefig('Temperature_Accuracy.svg')
    plt.show()

    return Tactual, Tmeas


def make_gaussian_process_samples_2(df):
    """
    Make a gaussian process fitting the Tactual-Tmeasured relationship
    :param df: pandas DataFrame with columns 'Temperature' (with the measured temperature)
                 and 'Tactual' (for the actual temperature)
    :return: emcee sampler instance
    """
    # First, find the uncertainties at each actual temperature
    temperature_groups = df.groupby('Temperature')
    Tactual = []
    Tmeasured = []
    error = []
    for T in temperature_groups.groups.keys():
        measurements = temperature_groups.get_group(T)['Tactual']
        if measurements.size < 10:
            continue
        l, m, h = np.percentile(measurements, [5.0, 50.0, 95.0])
        Tactual.append(m)
        Tmeasured.append(T)
        error.append(max(190.0, np.sqrt(((m-l)/2.0)**2 + ((h-m)/2.0)**2)))
    Tactual = np.array(Tactual)
    Tmeasured = np.array(Tmeasured)
    error = np.array(error)
    for Tm, Ta, e in zip(Tmeasured, Tactual, error):
        print Tm, Ta, e
    plt.figure(1)
    plt.errorbar(Tmeasured, Tactual, yerr=error, fmt='.k', capsize=0)
    plt.plot(Tmeasured, Tmeasured, 'r--')
    plt.xlim((min(Tmeasured) - 100, max(Tmeasured) + 100))
    plt.xlabel('Measured Temperature')
    plt.ylabel('Actual Temperature')
    plt.show(block=False)

    # Define some functions to use in the GP fit
    def model(pars, T):
        #polypars = pars[2:]
        #return np.poly1d(polypars)(T)
        return T

    def lnlike(pars, Tact, Tmeas, Terr):
        a, tau = np.exp(pars[:2])
        gp = george.GP(a * kernels.ExpSquaredKernel(tau))
        gp.compute(Tmeas, Terr)
        return gp.lnlikelihood(Tact - model(pars, Tmeas))

    def lnprior(pars):
        lna, lntau = pars[:2]
        polypars = pars[2:]
        if -20 < lna < 20 and 4 < lntau < 20:
            return 0.0
        return -np.inf

    def lnprob(pars, x, y, yerr):
        lp = lnprior(pars)
        return lp + lnlike(pars, x, y, yerr) if np.isfinite(lp) else -np.inf

    # Set up the emcee fitter
    initial = np.array([0, 6])#, 1.0, 0.0])
    ndim = len(initial)
    nwalkers = 100
    p0 = [np.array(initial) + 1e-8 * np.random.randn(ndim) for i in xrange(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(Tactual, Tmeasured, error))

    print 'Running first burn-in'
    p1, lnp, _ = sampler.run_mcmc(p0, 500)
    sampler.reset()

    print "Running second burn-in..."
    p_best = p1[np.argmax(lnp)]
    p2 = [p_best + 1e-8 * np.random.randn(ndim) for i in xrange(nwalkers)]
    p3, _, _ = sampler.run_mcmc(p2, 250)
    sampler.reset()

    print "Running production..."
    sampler.run_mcmc(p3, 1000)

    # Plot a bunch of the fits
    N = 100
    Tvalues = np.arange(3300, 7000, 20)
    for i in range(N):
        pars = sampler.flatchain[i]
        a, tau = np.exp(pars[:2])
        gp = george.GP(a * kernels.ExpSquaredKernel(tau))
        gp.compute(Tmeasured, error)
        s = gp.sample_conditional(Tactual - model(pars, Tmeasured), Tvalues) + model(pars, Tvalues)
        plt.plot(Tvalues, s, 'b-', alpha=0.1)
    plt.draw()



    return sampler

def make_gaussian_process_samples(df):
    """
    Make a gaussian process fitting the Tactual-Tmeasured relationship
    :param df: pandas DataFrame with columns 'Temperature' (with the measured temperature)
                 and 'Tactual' (for the actual temperature)
    :return: FILL IN ONCE I FIGURE THAT OUT!
    """
    # First, find the uncertainties at each actual temperature
    temperature_groups = df.groupby('Tactual')
    Tactual = []
    Tmeasured = []
    error = []
    for T in temperature_groups.groups.keys():
        measurements = temperature_groups.get_group(T)['Temperature']
        l, m, h = np.percentile(measurements, [16.0, 50.0, 84.0])
        Tactual.append(T)
        Tmeasured.append(m)
        error.append(np.sqrt((m-l)**2 + (h-m)**2))
    Tactual = np.array(Tactual)
    Tmeasured = np.array(Tmeasured)
    error = np.array(error)
    plt.figure(1)
    plt.errorbar(Tactual, Tmeasured, yerr=error, fmt='.k', capsize=0)
    plt.plot(Tactual, Tactual, 'r--')
    plt.xlim((min(Tactual) - 100, max(Tactual) + 100))
    plt.show(block=False)

    # Define some functions to use in the GP fit
    def model(pars, T):
        #polypars = pars[2:]
        #return np.poly1d(polypars)(T)
        return T

    def lnlike(pars, Tact, Tmeas, Terr):
        a, tau = np.exp(pars[:2])
        gp = george.GP(a * kernels.ExpSquaredKernel(tau))
        gp.compute(Tact, Terr)
        return gp.lnlikelihood(Tmeas - model(pars, Tact))

    def lnprior(pars):
        lna, lntau = pars[:2]
        polypars = pars[2:]
        if -20 < lna < 20 and 4 < lntau < 20:
            return 0.0
        return -np.inf

    def lnprob(pars, x, y, yerr):
        lp = lnprior(pars)
        return lp + lnlike(pars, x, y, yerr) if np.isfinite(lp) else -np.inf

    # Set up the emcee fitter
    initial = np.array([0, 6])#, 1.0, 0.0])
    ndim = len(initial)
    nwalkers = 100
    p0 = [np.array(initial) + 1e-8 * np.random.randn(ndim) for i in xrange(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(Tactual, Tmeasured, error))

    print 'Running first burn-in'
    p1, lnp, _ = sampler.run_mcmc(p0, 500)
    sampler.reset()

    print "Running second burn-in..."
    p_best = p1[np.argmax(lnp)]
    p2 = [p_best + 1e-8 * np.random.randn(ndim) for i in xrange(nwalkers)]
    p3, _, _ = sampler.run_mcmc(p2, 250)
    sampler.reset()

    print "Running production..."
    sampler.run_mcmc(p3, 1000)

    # Plot a bunch of the fits
    N = 100
    Tvalues = np.arange(3300, 7000, 20)
    for i in range(N):
        pars = sampler.flatchain[i]
        a, tau = np.exp(pars[:2])
        gp = george.GP(a * kernels.ExpSquaredKernel(tau))
        gp.compute(Tactual, error)
        s = gp.sample_conditional(Tmeasured - model(pars, Tactual), Tvalues) + model(pars, Tvalues)
        plt.plot(Tvalues, s, 'b-', alpha=0.1)
    plt.draw()



    return sampler #, p0, p1, p2, p3, truths, median , error


#def



if __name__ == '__main__':
    # Run these in an ipython shell for interactivity. Takes about 5 GB of RAM though!
    df = get_ccf_data('GeneratedObservations/Cross_correlations/')
    Tactual, Tmeas = plot_temperature_accuracy(df)
    Tactual = np.array(Tactual)
    Tmeas = np.array(Tmeas)
    sampler, p0, p1, p2, p3, truths, median , error = make_gaussian_process_samples(Tactual, Tmeas)

    #TODO: Sample the GP hyperparameters to get the measured temperature and uncertainty for each actual temperature. Then, reverse the relationship!
