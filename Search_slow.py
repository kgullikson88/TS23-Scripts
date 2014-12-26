import sys

import GenericSearch
import pandas
from astropy.io import fits


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
              #[433, 435],  #H gamma
              #[409, 411],  #H delta
              #[396, 398],  #H epsilon
              #[388, 390],  #H zeta
]

if "darwin" in sys.platform:
    modeldir = "/Volumes/DATADRIVE/Stellar_Models/PHOENIX/Stellar/Vband/"
elif "linux" in sys.platform:
    modeldir = "/media/FreeAgent_Drive/SyntheticSpectra/Sorted/Stellar/Vband/"
else:
    modeldir = raw_input("sys.platform not recognized. Please enter model directory below: ")
    if not modeldir.endswith("/"):
        modeldir = modeldir + "/"

if __name__ == '__main__':
    # Parse command line arguments:
    fileList = []
    trimsize = 10
    for arg in sys.argv[1:]:
        if 1:
            fileList.append(arg)
    prim_vsini = [100.0]*len(fileList)

    # Get the primary star vsini values
    prim_vsini = []
    vsini = pandas.read_csv("../../Useful_Datafiles/Vsini.csv", sep='|', skiprows=8)[1:]
    vsini_dict = {}
    for fname in fileList:
        root = fname.split('/')[-1][:9]
	if root in vsini_dict:
	    prim_vsini.append(vsini_dict[root])
        else:
            header = fits.getheader(fname)
            star = header['OBJECT']
            v = vsini.loc[vsini.Identifier.str.strip() == star]['vsini(km/s)'].values[0]
            prim_vsini.append(float(v) * 0.8)
            vsini_dict[root] = float(v) * 0.8
    for fname, vsini in zip(fileList, prim_vsini):
        print fname, vsini
    
    GenericSearch.slow_companion_search(fileList, prim_vsini,
	                                extensions=True,
					resolution=60000.0,
					trimsize=trimsize,
					modeldir=modeldir,
					badregions=badregions,
					metal_values=(0.0,),
					vsini_values=(1,),
					#Tvalues=range(3000, 6600, 100),
					observatory='McDonald',
					debug=False,
					vbary_correct=False,
					addmode='simple')

