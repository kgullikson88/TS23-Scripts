"""
Sensitivity analysis, using the new search method.
"""
import sys
import logging

import matplotlib.pyplot as plt

import Sensitivity
import StarData
import SpectralTypeRelations
import Search_slow
from HelperFunctions import ensure_dir


logging.basicConfig(level='INFO')

MS = SpectralTypeRelations.MainSequence()


def check_sensitivity():
    fileList = []
    for arg in sys.argv[1:]:
        if 1:
            fileList.append(arg)

    badregions = Search_slow.badregions
    interp_regions = Search_slow.interp_regions
    trimsize = Search_slow.trimsize
    prim_vsini = StarData.get_vsini(fileList)

    Sensitivity.Analyze(fileList, prim_vsini,
                        hdf5_file='/media/ExtraSpace/PhoenixGrid/TS23_Grid.hdf5',
                        extensions=True,
                        resolution=None,
                        trimsize=trimsize,
                        badregions=badregions, interp_regions=interp_regions,
                        metal_values=(0.0,),
                        vsini_values=(100, 150, 200, 250),
                        Tvalues=range(7000, 12100, 1000),
                        debug=False,
                        addmode='all',
                        output_mode='hdf5',
                        output_file='Sensitivity.hdf5')


if __name__ == '__main__':
    if '--analyze' in sys.argv[1]:
        # Make the 2d plots
        df = Sensitivity.analyze_sensitivity(hdf5_file='Sensitivity.hdf5', interactive=False, update=False)

    elif '--marginalize' in sys.argv[1]:
        fig, ax = Sensitivity.marginalize_sensitivity(infilename='Sensitivity_Dataframe.csv')
        # plt.show()
        ensure_dir('Figures/')
        plt.savefig('Figures/Sensitivity_Marginalized.pdf')


    else:
        check_sensitivity()
