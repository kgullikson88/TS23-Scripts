import numpy
import sys
import os
import matplotlib.pyplot as plt
import pyfits
from scipy.interpolate import InterpolatedUnivariateSpline as interp
import TelluricFitter
import FitsUtils
import DataStructures
import Units
import FindContinuum

homedir = os.environ["HOME"]
weather_file = homedir + "/School/Research/Useful_Datafiles/Weather.dat"
linelist = homedir + "/School/Research/Useful_Datafiles/Linelist_visible.dat"
telluric_orders = [3,4,5,6,8,9,10,11,14]


"""
  Function to output a fits file
  column_dict is a dictionary where the key is the name of the column
     and the value is a numpy array with the data. Example of a column
     would be the wavelength or flux at each pixel
  template is the filename of the template fits file. The header will
     be taken from this file and used as the main header
  mode determines how the outputted file is made. Append will just add
     a fits extension to the existing file (and then save it as outfilename)
     "new" mode will create a new fits file. 
"""
def OutputFitsFile(column_dict, template, outfilename, mode="append"):
  #Get header from template. Use this in the new file
  if mode == "new":
    hdulist = pyfits.open(template)
    header = hdulist[0].header
    hdulist.close()

  columns = []
  for key in column_dict.keys():
    columns.append(pyfits.Column(name=key, format="D", array=column_dict[key]))
  cols = pyfits.ColDefs(columns)
  tablehdu = pyfits.new_table(cols)

  if mode == "append":
    hdulist = pyfits.open(template)
    hdulist.append(tablehdu)
  elif mode == "new":
    hdu = pyfits.PrimaryHDU(header=header)
    hdulist = pyfits.HDUList([hdu, tablehdu])
  hdulist.writeto(outfilename, clobber=True, output_verify='ignore')
  hdulist.close()


if __name__ == "__main__":
  #Initialize fitter
  fitter = TelluricFitter.TelluricFitter()
  fitter.SetTelluricLineListFile(linelist)

  #Find and read in blaze function (MUST BE IN CURRENT DIRECTORY!)
  files = os.listdir("./")
  blazefile = [fname for fname in files if fname.startswith("BLAZE")][0]
  blaze_orders = FitsUtils.MakeXYpoints(blazefile, errors=2)
  blaze_functions = []
  blaze_errors = []
  for order in blaze_orders:
    blaze_functions.append( interp(order.x, order.y, k=1) )
    blaze_errors.append( interp(order.x, order.err, k=1) )


  #START LOOPING OVER INPUT FILES
  for fname in sys.argv[1:]:
    num = fname[2:-5]
    outfilename = "Corrected_%s.fits" %num

    orders = FitsUtils.MakeXYpoints(fname, errors=2)
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
    fitter.SetBounds({"h2o": [1.0, 99.0],
                      "o2": [5e4, 1e6],
                      "resolution": [10000,100000]})
    models = []

    #START LOOPING OVER ORDERS
    start = 2
    for i, order in enumerate(orders[start:]):
      fitter.AdjustValue({"wavestart": order.x[0] - 20.0,
                          "waveend": order.x[-1] + 20.0})
      order.y /= blaze_functions[i+start](order.x)

      #Error propagation not right, but (hopefully) close enough!
      order.err /= blaze_functions[i+start](order.x)
      #order.err = numpy.sqrt( (order.err/blaze_functions[i+start](order.x))**2 + (order.y/(blaze_functions[i+start](order.x))**2 * blaze_errors[i+start](order.x))**2 )
      
      order.cont = FindContinuum.Continuum(order.x, order.y, fitorder=3, lowreject=2, highreject=10)

      if i+start not in telluric_orders:
        print "Skipping order %i" %(i+start)
        data = order.copy()
        model = DataStructures.xypoint(x=order.x.copy(), y=numpy.ones(order.x.size))
      else:        
        fitter.ImportData(order)
        model = fitter.Fit(resolution_fit_mode="SVD", fit_primary=True, adjust_wave="model")
        models.append(model)
        data = fitter.data
        
      columns = {"wavelength": data.x,
	         "flux": data.y,
                 "continuum": data.cont,
                 "error": data.err,
		 "model": model.y}
      if i == 0:
        OutputFitsFile(columns, fname, outfilename)
      else:
        OutputFitsFile(columns, outfilename, outfilename)
