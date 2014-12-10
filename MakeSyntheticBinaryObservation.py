"""
This script takes several of my program stars, and several late-type stars, and creates
several synthetic spectra by combining them. Use these data to test my CCF technique and
the parameters it derives!

Usage:
python MakeSyntheticBinaryObservation.py myfile1 myfile2 ... -late latefile1 latefile2...

myfile: a filename for my data. Should be a KG* file
latefile: a filename for a late-type star with known properties, taken with the same instrument as above

"""

import HelperFunctions
import sys
import os

import FittingUtilities
from astropy.io import fits
import numpy as np
from collections import defaultdict
import logging
from astroquery.simbad import Simbad


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
            total = combine(early_file, late_file)

