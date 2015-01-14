import sys

from astropy.io import fits as pyfits


if __name__ == "__main__":
    for fname in sys.argv[1:]:
        hdulist = pyfits.open(fname)
        header = hdulist[0].header
        # header = pyfits.getheader(fname)
        print fname, header["OBJECT"], header['DATE-OBS'], header['UT']

