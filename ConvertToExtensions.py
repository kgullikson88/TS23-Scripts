import sys
import os

import FittingUtilities
from astropy.io import fits as pyfits
import numpy as np

import HelperFunctions


left_trim = 0
right_trim = 50
bad_regions = {12: [438, 9e9],
               13: [442.5, 9e9],
               14: [447.5, 9e9],
               15: [455, 9e9],
               16: [1, 9e9],
               17: [1, 9e9],
               18: [1, 9e9],
               58: [1, 9e9],
               59: [1, 9e9]
}


def read_orders(fname, blazefile=None):
    # Read in the data
    try:
        orders = HelperFunctions.ReadFits(fname)
    except ValueError:
        orders = HelperFunctions.ReadFits(fname, errors=2)
    orders = orders[::-1]  # Reverse order so the bluest order is first
    print len(orders)

    # Read in the blaze, if applicable
    if blazefile is not None:
        try:
            blaze = HelperFunctions.ReadFits(blazefile)
        except ValueError:
            blaze = HelperFunctions.ReadFits(blazefile, errors=2)
        except IOError:
            print "Error! blaze file %s does not exist!" % blazefile
            print "Not converting file %s" % fname
            raise IOError
        blaze = blaze[::-1]
        print len(blaze)

    # Trim the data
    order_list = []
    for i, order in enumerate(orders):

        left, right = left_trim, order.size() - right_trim
        if i in bad_regions.keys():
            region = bad_regions[i]
            left = np.searchsorted(order.x, region[0])
            right = np.searchsorted(order.x, region[1])
            if left == 0 or right == order.size():
                order.x = np.delete(order.x, np.arange(left, right))
                order.y = np.delete(order.y, np.arange(left, right))
                order.cont = np.delete(order.cont, np.arange(left, right))
                order.err = np.delete(order.err, np.arange(left, right))
                if blazecorrect:
                    blaze[i].y = np.delete(blaze[i].y, np.arange(left, right))
            else:
                print "Warning! Bad region covers the middle of order %i" % i
                print "Interpolating rather than removing"
                order.y[left:right] = order.cont[left:right]
                order.err[left:right] = 9e9
        else:
            order = order[left:right]
            if blazecorrect:
                blaze[i].y = blaze[i].y[left:right]
        if order.size() < 10:
            continue


        # Blaze correction
        if blazecorrect:
            order.y /= blaze[i].y
            order.err /= blaze[i].y

        order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=3, lowreject=1.5, highreject=5)

        order_list.append(order.copy())
    return order_list


if __name__ == "__main__":
    fileList = []
    blazecorrect = True
    for arg in sys.argv[1:]:
        if "noblaze" in arg:
            blazecorrect = False
        else:
            fileList.append(arg)
    for fname in fileList:
        outfilename = "%s-0.fits" % (fname.split(".fits")[0])

        if blazecorrect:
            header = pyfits.getheader(fname)
            try:
                blazefile = "%s.fits" % header['BLAZE']
            except KeyError:
                allfiles = os.listdir("./")
                blazefile = [f for f in allfiles if "BLAZE" in f][0]
        else:
            blazefile = None
        orders = read_orders(fname, blazefile)

        column_list = []
        for i, order in enumerate(orders):
            columns = {"wavelength": order.x,
                       "flux": order.y,
                       "continuum": order.cont,
                       "error": order.err}
            column_list.append(columns)
        HelperFunctions.OutputFitsFileExtensions(column_list, fname, outfilename, mode="new")
      
      
