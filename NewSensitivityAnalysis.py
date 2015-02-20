"""
Sensitivity analysis, using the new search method.
"""
import sys

import Sensitivity
import StarData

import Search_slow


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
                        vsini_values=(5,),
                        Tvalues=range(3000, 6000, 100),
                        debug=False,
                        addmode='simple',
                        output_mode='hdf5')


def analyze_sensitivity():
    """
    This uses the output of a previous run of check_sensitivity, and makes plots
    :return:
    """


if __name__ == '__main__':
    if 'analyze' in sys.argv[1]:
        analyze_sensitivity()
    else:
        check_sensitivity()