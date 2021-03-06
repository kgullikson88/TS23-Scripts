import sys
import os
from scipy.interpolate import InterpolatedUnivariateSpline as interp

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
import TelluricFitter
import DataStructures

import FitsUtils
import Units
import FindContinuum


homedir = os.environ["HOME"]
weather_file = homedir + "/School/Research/Useful_Datafiles/Weather.dat"
linelist = homedir + "/School/Research/Useful_Datafiles/Linelist_visible.dat"
telluric_orders = [3, 4, 5, 6, 8, 9, 10, 11, 13, 14, 15, 16, 17, 19, 20, 24, 25]

"""
  Function to output a fits file
  column_dict is a dictionary where the key is the name of the column
     and the value is a np array with the data. Example of a column
     would be the wavelength or flux at each pixel
  template is the filename of the template fits file. The header will
     be taken from this file and used as the main header
  mode determines how the outputted file is made. Append will just add
     a fits extension to the existing file (and then save it as outfilename)
     "new" mode will create a new fits file. 
     header_info takes a list of lists. Each sub-list should have size 2 where the first element is the name of the new keyword, and the second element is the corresponding value. A 3rd element may be added as a comment
"""


def EditFitsFile(column_dict, filename, extension, header_info=[]):
    print "Editing extension number %i" % extension
    columns = []
    for key in column_dict.keys():
        columns.append(pyfits.Column(name=key, format="D", array=column_dict[key]))
    cols = pyfits.ColDefs(columns)
    tablehdu = pyfits.new_table(cols)

    # Add keywords to extension header
    num_keywords = len(header_info)
    header = tablehdu.header
    for i in range(num_keywords):
        info = header_info[i]
        if len(info) > 2:
            header.update(info[0], info[1], info[2])
        elif len(info) == 2:
            header.update(info[0], info[1])

    #Open file and update the appropriate extension
    hdulist = pyfits.open(filename, mode='update', save_backup=True)
    if extension < len(hdulist):
        hdulist[extension] = tablehdu
    else:
        hdulist.append(tablehdu)
    hdulist.flush()
    hdulist.close()

    return


if __name__ == "__main__":
    # Initialize fitter
    fitter = TelluricFitter.TelluricFitter()
    fitter.SetTelluricLineListFile(linelist)
    LineList = np.loadtxt(linelist)
    logfile = open("fitlog.txt", "w")

    #Find and read in blaze function (MUST BE IN CURRENT DIRECTORY!)
    files = os.listdir("./")
    blazefile = [fname for fname in files if fname.startswith("BLAZE")][0]
    blaze_orders = FitsUtils.MakeXYpoints(blazefile, errors=2)
    blaze_functions = []
    blaze_errors = []
    for order in blaze_orders:
        blaze_functions.append(interp(order.x, order.y))
        blaze_errors.append(interp(order.x, order.err))

    #filename is the name of the file to re-correct, and should be given first
    fname = sys.argv[1]
    header = pyfits.getheader(fname)
    orders = FitsUtils.MakeXYpoints(fname, errors=2)

    date = header["DATE-OBS"]
    time = header["UT"]
    t_seg = time.split(":")
    time = 3600 * float(t_seg[0]) + 60 * float(t_seg[1]) + float(t_seg[2])

    #Read in weather information (archived data is downloaded from weather.as.utexas.edu
    infile = open(weather_file)
    lines = infile.readlines()
    infile.close()
    times = []
    RH = []
    P = []
    T = []
    idx = 0
    bestindex = 0
    difference = 9e9
    for line in lines[1:]:
        segments = line.split()
        if date in segments[0]:
            segments = segments[1].split(",")
            t = segments[0]
            t_seg = t.split(":")
            weather_time = 3600 * float(t_seg[0]) + 60 * float(t_seg[1]) + float(t_seg[2])
            if np.abs(time - weather_time) < difference:
                difference = np.abs(time - weather_time)
                bestindex = idx
            times.append(segments[0])
            T.append(float(segments[3]))
            RH.append(float(segments[4]))
            P.append(float(segments[5]))
            idx += 1

    angle = float(header["zd"])
    resolution = 60000.0
    humidity = RH[bestindex]
    humidity = 40.0
    T_fahrenheit = T[bestindex]
    pressure = P[bestindex] * Units.hPa / Units.inch_Hg
    temperature = (T_fahrenheit - 32.0) * 5.0 / 9.0 + 273.15

    #Adjust fitter values
    fitter.AdjustValue({"angle": angle,
                        "temperature": temperature,
                        "pressure": pressure,
                        "resolution": resolution})
    fitter.SetBounds({"h2o": [1.0, 96.0],
                      "o2": [5e4, 1e6],
                      "resolution": [40000, 80000]})

    #Loop over the orders to correct (should be the 2nd through last argument to the call):
    for i in sys.argv[2:]:
        i = int(i)
        print "Re-fitting order %i" % i
        fitter.FitVariable({"h2o": humidity,
                            "o2": 2.12e5})

        order = orders[i]
        fitter.AdjustValue({"wavestart": order.x[0] - 20.0,
                            "waveend": order.x[-1] + 20.0})
        order.y /= blaze_functions[i](order.x)

        #Error propagation not right, but (hopefully) close enough!
        order.err /= blaze_functions[i](order.x)
        #order.err = np.sqrt( (order.err/blaze_functions[i+start](order.x))**2 + (order.y/(blaze_functions[i+start](order.x))**2 * blaze_errors[i+start](order.x))**2 )

        order.cont = FindContinuum.Continuum(order.x, order.y, fitorder=3, lowreject=2, highreject=10)
        primary = DataStructures.xypoint(x=order.x, y=np.ones(order.x.size))

        fitter.ImportData(order)
        fitpars = [fitter.const_pars[j] for j in range(len(fitter.parnames)) if fitter.fitting[j]]
        model = fitter.GenerateModel(fitpars, LineList)
        fitter.ImportData(order)  #Re-initialize to original data before fitting
        model_amplitude = 1.0 - min(model.y)
        if model_amplitude < 0.01 or i > 29:
            print "Skipping order %i" % (i)
            #model = DataStructures.xypoint(x=order.x.copy(), y=np.ones(order.x.size))
        elif model_amplitude >= 0.01 and model_amplitude < 0.1:
            print "Fitting line profiles with gaussian profile"
            model = fitter.Fit(resolution_fit_mode="gauss", fit_primary=False, adjust_wave="model")
        else:
            print "Large model amplitude. Using SVD for line profiles"
            model = fitter.Fit(resolution_fit_mode="SVD", fit_primary=False, adjust_wave="model")
        data = fitter.data

        #Set up data structures for OutputFitsFile
        columns = {"wavelength": data.x,
                   "flux": data.y,
                   "continuum": data.cont,
                   "error": data.err,
                   "model": model.y,
                   "primary": primary.y}
        namedict = {"pressure": ["PRESFIT", "PRESVAL", "Pressure"],
                    "temperature": ["TEMPFIT", "TEMPVAL", "Temperature"],
                    "angle": ["ZD_FIT", "ZD_VAL", "Zenith Distance"],
                    "resolution": ["RESFIT", "RESVAL", "Detector Resolution"],
                    "h2o": ["H2OFIT", "H2OVAL", "H2O abundance"],
                    "co2": ["CO2FIT", "CO2VAL", "CO2 abundance"],
                    "o3": ["O3FIT", "O3VAL", "O3 abundance"],
                    "n2o": ["N2OFIT", "N2OVAL", "N2O abundance"],
                    "co": ["COFIT", "COVAL", "CO abundance"],
                    "ch4": ["CH4FIT", "CH4VAL", "CH4 abundance"],
                    "o2": ["O2FIT", "O2VAL", "O2 abundance"],
                    "no": ["NOFIT", "NOVAL", "NO abundance"],
                    "so2": ["SO2FIT", "SO2VAL", "SO2 abundance"],
                    "no2": ["NO2FIT", "NO2VAL", "NO2 abundance"],
                    "nh3": ["NH3FIT", "NH3VAL", "NH3 abundance"],
                    "hno3": ["HNO3FIT", "HNO3VAL", "HNO3 abundance"]}
        header_info = []
        numpars = len(fitter.const_pars)
        for j in range(numpars):
            try:
                parname = fitter.parnames[j]
                parval = fitter.const_pars[j]
                fitting = fitter.fitting[j]
                header_info.append([namedict[parname][0], fitting, namedict[parname][2]])
                header_info.append([namedict[parname][1], parval, namedict[parname][2]])
            except KeyError:
                print "Not saving the following info: %s" % (fitter.parnames[j])

        EditFitsFile(columns, fname, i - 1, header_info=header_info)
  
