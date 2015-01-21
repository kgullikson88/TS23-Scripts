"""
Use the outputs of Search_slow to find the best temperature of the secondary star.
"""

import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.interpolate import InterpolatedUnivariateSpline as spline


VEL_GRID = np.arange(-900, 900, 0.1)
COLS = ['original', 'method', 'T', 'logg', 'metallicity', 'vsini']


def classify_file(fname):
    """
    Classify a filename
    :param fname: the filename to classify
    :return: the temperature, logg, metallicity, vsini, ccf method, and original file
    """
    fname = fname.split('/')[-1]
    original_file = fname[:7]
    method = fname.split('-method')[0].split('_')[-1]
    vsini = float(fname.split('kps')[0].split('.')[-1])
    T = float(fname.split('K+')[0].split('_')[-1])
    logg = float(fname[-8:-4])
    metallicity = float(fname[-4:])
    return T, logg, metallicity, vsini, method, original_file


def get_ccf_data(ccf_dir):
    """
    Reads in the ccf files in ccf_dir, and classifies them
    :param ccf_dir: directory containing the ccf files
    :return: pandas DataFrame with the information
    """
    if not ccf_dir.endswith('/'):
        ccf_dir += '/'
    file_list = ['{}{}'.format(ccf_dir, f) for f in os.listdir(ccf_dir)]

    data = defaultdict(list)
    for fname in file_list:
        T, logg, metallicity, vsini, method, original_file = classify_file(fname)
        vel, corr = np.loadtxt(fname, unpack=True)
        ccf = spline(vel, corr)(VEL_GRID)
        data['T'].append(T)
        data['logg'].append(logg)
        data['metallicity'].append(metallicity)
        data['vsini'].append(vsini)
        data['method'].append(method)
        data['original'].append(original_file)
        data['ccf'].append(ccf)

    return pd.DataFrame(data=data)


def find_best_pars(df):
    """
    Finds the 'best-fit' parameters for each original file and method
    :param df:
    :return:
    """
    # First, get the maximum value of the ccf
    df['max_ccf'] = df['ccf'].map(np.max)

    methods = pd.unique(df.method)
    original_files = pd.unique(df.original)
    best_info = defaultdict(list)
    for original_filename in original_files:
        for method in methods:
            good = df.loc[(df.method == method) & (df.original == original_filename)]
            best = good.loc[good['max_ccf'] == good['max_ccf'].max()]
            # print 'File: {}\n\tmethod = {}\n\tT = {}\n\tlogg = {}\n\t[Fe/H] = {}'.format(original_filename,
            #                                                                             method,
            #                                                                             best['T'].item(),
            #                                                                             best['logg'].item(),
            #                                                                             best['metallicity'].item())
            #print '\tvsini = {}'.format(best['vsini'].item())
            best_info['original'].append(original_filename)
            best_info['method'].append(method)
            best_info['T'].append(best['T'].item())
            best_info['logg'].append(best['logg'].item())
            best_info['metallicity'].append(best['metallicity'].item())
            best_info['vsini'].append(best['vsini'].item())

    return pd.DataFrame(data=best_info)


if __name__ == '__main__':
    data = get_ccf_data(sys.argv[1])
    best = find_best_pars(data)
    print best[COLS]
    best.to_csv('CCF_fit.csv', index=False)