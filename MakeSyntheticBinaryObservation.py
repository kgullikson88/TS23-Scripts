"""
This script takes several of my program stars, and several late-type stars, and creates
several synthetic spectra by combining them. Use these data to test my CCF technique and
the parameters it derives!

Usage:
python MakeSyntheticBinaryObservation.py myfile1 myfile2 ... -late latefile1 latefile2...

myfile: a filename for my data. Should be a KG* file
latefile: a filename for a late-type star with known properties, taken with the same instrument as above

"""

import sys
import os
import logging
from scipy.interpolate import InterpolatedUnivariateSpline as spline

from astropy.io import fits
import numpy as np
from astroquery.simbad import Simbad

import HelperFunctions
import ConvertToExtensions


def parse_input(argv):
    early_files = []
    late_files = []
    late = False #Flag for when to start filling the late file list
    for arg in argv:
        if arg.lower() == '-late':
            late = True
        elif late:
            late_files.append(arg)
        else:
            early_files.append(arg)
    if len(early_files) < 1 or len(late_files) < 1:
        logging.warning('Need to give at least one early-type (mine) and one late-type spectrum, '
                        'separated by the "-late" flag')

    return early_files, late_files


def combine(early_filename, late_filename):
    """
    This function does the bulk of the work.
    :param early_filename: the filename of the early-type star
    :param late_filename: the filename of the late-type star. Should be from the same instrument and setup as my data!
    :return: A list of orders, after combination
    """
    # First, we will classify the files by their fits header and some simbad lookups
    early_dict = classify_file(early_filename)
    late_dict = classify_file(late_filename)
    print early_dict
    print late_dict

    # Now, we need to work out how to scale the data as if they were in the same system
    # This ASSUMES that all atmospheric affects are the same for both observations, which is not really true!
    scale = (early_dict['exptime'] / late_dict['exptime']) * (early_dict['plx'] / late_dict['plx']) ** 2
    print scale

    # Get the blazefile for my data
    header = fits.getheader(early_filename)
    try:
        blazefile = "%s.fits" % header['BLAZE']
    except KeyError:
        allfiles = os.listdir("./")
        blazefile = [f for f in allfiles if "BLAZE" in f][0]

    # Read in the orders for both files
    early_orders = ConvertToExtensions.read_orders(early_filename, blazefile)
    late_orders = ConvertToExtensions.read_orders(late_filename)

    # Finally, combine:
    order_list = []
    for order in early_orders:
        num = HelperFunctions.FindOrderNums(late_orders, [np.median(order.x)])[0]
        late = late_orders[num]
        late.y *= scale
        combined = add_order(order, late)
        order_list.append(combined)
        # plt.plot(order.x, order.y)
        #plt.plot(combined.x, combined.y)
    # plt.show()

    return order_list, early_dict, late_dict


def add_order(early, late):
    """
    Adds two xypoint instances. This is only really necessary because the wavelengths
    and/or order lengths are not always the same...
    :param early: the early type star
    :param late: the late type star
    :return: xypoint instance of the combined 'observation'
    """
    left = np.searchsorted(early.x, late.x[0])
    right = np.searchsorted(early.x, late.x[-1])
    late_fcn = spline(late.x, late.y, k=1)
    combined = early[left:right].copy()
    combined.y += late_fcn(combined.x)
    return combined


def classify_file(filename):
    """
    This function uses the fits header information and the Simbad database to classify the object
    :param filename: The filename of the observation to be classified
    :return:
    """
    # Read in the header and get the object name
    header = fits.getheader(filename)
    object = header['object']

    # Make a Simbad object
    sim = Simbad()
    sim.add_votable_fields('plx', 'sp')
    data = sim.query_object(object)
    plx = data['PLX_VALUE'].item()
    spt = data['SP_TYPE'].item()

    d = {'Object': object,
         'plx': plx,
         'SpT': spt,
         'exptime': header['exptime']}
    return d



if __name__ == '__main__':
    early, late = parse_input(sys.argv[1:])

    # Add each late file to all of the early-type files
    for late_file in late:
        for early_file in early:
            total, early_dict, late_dict = combine(early_file, late_file)

            # Prepare for output
            HelperFunctions.ensure_dir('GeneratedObservations')
            outfilename = 'GeneratedObservations/{}_{}.fits'.format(early_file.split('.fits')[0], 1)
            column_list = []
            for order in total:
                column = {'wavelength': order.x,
                          'flux': order.y,
                          'continuum': order.cont,
                          'error': order.err}
                column_list.append(column)
            newheader = {'object1': early_dict['Object'],
                         'object2': late_dict['Object'],
                         'SpT1': early_dict['SpT'],
                         'SpT2': late_dict['SpT'],
                         'file1': early_file,
                         'file2': late_file}
            HelperFunctions.OutputFitsFileExtensions(column_list, early_file, outfilename,
                                                     mode='new', primary_header=newheader)

