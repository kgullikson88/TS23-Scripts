import sys

import GenericSearch


# Define regions contaminated by telluric residuals or other defects. We will not use those regions in the cross-correlation
badregions = [[567.5, 575.5],
              [588.5, 598.5],
              [627, 632],
              [647, 655],
              [686, 706],
              [716, 734],
              [759, 9e9],
              [655, 657],  # H alpha
              [485, 487],  #H beta
              [433, 435],  #H gamma
              [409, 411],  #H delta
              [396, 398],  #H epsilon
              [388, 390],  #H zeta
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
    extensions = True
    tellurics = False
    trimsize = 10
    for arg in sys.argv[1:]:
        if '-e' in arg:
            extensions = False
        if '-t' in arg:
            tellurics = True  #telluric lines modeled but not removed
        else:
            fileList.append(arg)

    GenericSearch.CompanionSearch(fileList,
                                  extensions=extensions,
                                  resolution=60000.0,
                                  trimsize=trimsize,
                                  modeldir=modeldir)


