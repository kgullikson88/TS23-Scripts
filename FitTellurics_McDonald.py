import numpy
import sys
import os
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
from scipy.interpolate import InterpolatedUnivariateSpline as interp
import TelluricFitter
import FitsUtils
import DataStructures
import Units
import FindContinuum

homedir = os.environ["HOME"]
weather_file = homedir + "/School/Research/Useful_Datafiles/Weather.dat"
linelist = homedir + "/School/Research/Useful_Datafiles/Linelist_visible.dat"
telluric_orders = [3,4,5,6,8,9,10,11,13,14,15,16,17,19,20,24,25]



if __name__ == "__main__":
  #Initialize fitter
  fitter = TelluricFitter.TelluricFitter()
  fitter.SetTelluricLineListFile(linelist)
  LineList = numpy.loadtxt(linelist)
  logfile = open("fitlog.txt", "w")

  #Find and read in blaze function (MUST BE IN CURRENT DIRECTORY!)
  files = os.listdir("./")
  blazefile = [fname for fname in files if fname.startswith("BLAZE")][0]
  blaze_orders = FitsUtils.MakeXYpoints(blazefile, errors=2)
  blaze_functions = []
  blaze_errors = []
  for order in blaze_orders:
    blaze_functions.append( interp(order.x, numpy.ones(order.size())) ) #order.y) )
    blaze_errors.append( interp(order.x, order.err) )


  #START LOOPING OVER INPUT FILES
  for fname in sys.argv[1:]:
    num = fname[2:].split(".fits")[0]
    outfilename = "Corrected_%s.fits" %num

    orders = FitsUtils.MakeXYpoints(fname, errors="error", extensions=True, x="wavelength", y="flux")
    
    hdulist = pyfits.open(fname)
    header = hdulist[0].header
    hdulist.close()
    
    date = header["DATE-OBS"]
    time = header["UT"]
    t_seg = time.split(":")
    time = 3600*float(t_seg[0]) + 60*float(t_seg[1]) + float(t_seg[2])
    
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
        weather_time = 3600*float(t_seg[0]) + 60*float(t_seg[1]) + float(t_seg[2])
        if numpy.abs(time - weather_time) < difference:
          difference = numpy.abs(time - weather_time)
          bestindex = idx
        times.append(segments[0])
        T.append(float(segments[3]))
        RH.append(float(segments[4]))
        P.append(float(segments[5]))
        idx += 1
    
    
    angle = float(header["zd"])
    resolution = 60000.0
    humidity = RH[bestindex]
    T_fahrenheit = T[bestindex]
    pressure = P[bestindex]*Units.hPa/Units.inch_Hg
    temperature = (T_fahrenheit - 32.0)*5.0/9.0 + 273.15
    
    #Adjust fitter values
    fitter.FitVariable({"h2o": humidity, 
                        "o2": 2.12e5})
    fitter.AdjustValue({"angle": angle,
                        "temperature": temperature,
                        "pressure": pressure,
                        "resolution": resolution})
    fitter.SetBounds({"h2o": [1.0, 96.0],
                      "o2": [5e4, 1e6],
                      "resolution": [30000,90000]})
    models = []

    #START LOOPING OVER ORDERS
    start = 2\
    for i, order in enumerate(orders[start:]):
      print "\n***************************\nFitting order %i: " %(i+start)
      fitter.AdjustValue({"wavestart": order.x[0] - 20.0,
                          "waveend": order.x[-1] + 20.0})
      
      order.cont = FindContinuum.Continuum(order.x, order.y, fitorder=3, lowreject=2, highreject=10)
      primary = DataStructures.xypoint(x=order.x, y=numpy.ones(order.x.size))

      
      fitter.ImportData(order)
      fitpars = [fitter.const_pars[j] for j in range(len(fitter.parnames)) if fitter.fitting[j] ]
      model = fitter.GenerateModel(fitpars, LineList)
      fitter.ImportData(order) #Re-initialize to original data before fitting
      model_amplitude = 1.0 - min(model.y)
      #if i+start not in telluric_orders:
      if model_amplitude < 0.01 or i+start > 29:
        print "Skipping order %i" %(i+start)
        data = order.copy()
        #model = DataStructures.xypoint(x=order.x.copy(), y=numpy.ones(order.x.size))
        primary = model.copy()
      elif model_amplitude >= 0.01 and model_amplitude < 0.1:        
        print "Fitting line profiles with gaussian profile"
        model = fitter.Fit(resolution_fit_mode="gauss", fit_primary=False, adjust_wave="model")
        models.append(model)
        data = fitter.data
      else: 
        print "Large model amplitude. Using SVD for line profiles"
        model = fitter.Fit(resolution_fit_mode="SVD", fit_primary=False, adjust_wave="model")
        models.append(model)
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
          header_info.append([namedict[parname][0], fitting, namedict[parname][2] ])
          header_info.append([namedict[parname][1], parval, namedict[parname][2] ])
        except KeyError:
          print "Not saving the following info: %s" %(fitter.parnames[j])
      
      
      if i == 0:
        OutputFitsFile(columns, fname, outfilename, header_info=header_info)
      else:
        OutputFitsFile(columns, outfilename, outfilename, header_info=header_info)
