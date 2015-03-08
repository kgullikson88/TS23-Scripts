import sys
import logging
import os

import pandas
from astropy.io import fits

import GenericSearch
import StarData







# Define regions contaminated by telluric residuals or other defects. We will not use those regions in the cross-correlation
badregions = [[567.5, 575.5],
              [588.5, 598.5],
              [627, 632],
              [647, 655],
              [686, 706],
              [716, 734],
              [759, 9e9],
              # [655, 657],  # H alpha
              # [485, 487],  #H beta
              # [433, 435],  #H gamma
              #[409, 411],  #H delta
              #[396, 398],  #H epsilon
              #[388, 390],  #H zeta
]
interp_regions = []
trimsize = 10

if "darwin" in sys.platform:
    modeldir = "/Volumes/DATADRIVE/Stellar_Models/Sorted/Stellar/Vband/"
elif "linux" in sys.platform:
    modeldir = "/media/FreeAgent_Drive/SyntheticSpectra/Sorted/Stellar/Vband/"
else:
    modeldir = raw_input("sys.platform not recognized. Please enter model directory below: ")
    if not modeldir.endswith("/"):
        modeldir = modeldir + "/"


def get_primary_vsini(file_list):
    homedir = os.environ['HOME']
    vsini = pandas.read_csv("{}/School/Research/Useful_Datafiles/Vsini.csv".format(homedir), sep='|', skiprows=8)[1:]
    vsini_dict = {}
    prim_vsini = []
    for fname in file_list:
        root = fname.split('/')[-1][:9]
        if root in vsini_dict:
            prim_vsini.append(vsini_dict[root])
        else:
            header = fits.getheader(fname)
            star = header['OBJECT']
            try:
                v = vsini.loc[vsini.Identifier.str.strip() == star]['vsini(km/s)'].values[0]
                prim_vsini.append(float(v) * 0.8)
                vsini_dict[root] = float(v) * 0.8
            except IndexError:
                logging.warn('No vsini found for star {}! No primary star removal will be attempted!'.format(star))
                prim_vsini.append(None)
    for fname, vsini in zip(file_list, prim_vsini):
        print fname, vsini
    return prim_vsini


if __name__ == '__main__':
    # Parse command line arguments:
    fileList = []
    for arg in sys.argv[1:]:
        if 1:
            fileList.append(arg)

    # Get the primary star vsini values
    prim_vsini = StarData.get_vsini(fileList)

    GenericSearch.slow_companion_search(fileList, prim_vsini,
                                        hdf5_file='/media/ExtraSpace/PhoenixGrid/TS23_Grid.hdf5',
                                        extensions=True,
                                        resolution=None,
                                        trimsize=trimsize,
                                        modeldir=modeldir,
                                        badregions=badregions,
                                        metal_values=(0.0, -0.5, 0.5),
                                        vsini_values=(1, 5, 10, 20, 30),
                                        Tvalues=range(3000, 7500, 100),
                                        observatory='McDonald',
                                        debug=False,
                                        vbary_correct=True,
                                        addmode='simple',
                                        obstype='real',
                                        output_mode='hdf5')

    """   Use this for the synthetic companion search
    GenericSearch.slow_companion_search(fileList, prim_vsini,
                                        hdf5_file='/media/ExtraSpace/PhoenixGrid/TS23_Grid.hdf5',
                                        extensions=True,
                                        resolution=None,
                                        trimsize=trimsize,
                                        modeldir=modeldir,
                                        badregions=badregions,
                                        metal_values=(0.0, -0.5, 0.5),
                                        vsini_values=(1, 5, 10, 20, 30),
                                        Tvalues=range(3000, 6900, 100),
                                        observatory='McDonald',
                                        debug=False,
                                        vbary_correct=False,
                                        addmode='T-weighted',
                                        obstype='synthetic',
                                        output_mode='hdf5')

    """
