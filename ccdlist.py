import sys

from astropy.io import fits as pyfits


if __name__ == "__main__":
    for fname in sys.argv[1:]:
        header = pyfits.getheader(fname)
        print fname, header["OBJECT"], header["date-obs"], header['ut'], header['EXPTIME']

