import FitsUtils
import numpy
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import DataStructures
import FindContinuum
import MakeModel
import sys
import time
import os
import FitTellurics_McDonald
import subprocess
import pyfits


class LineFitter:
  def __init__(self, infilename, telluricfile = "/Users/kgulliks/School/Research/aerlbl_v12.2/rundir2/OutputModels/transmission-792.30-290.93-45.0-7.4-368.50-4.00-10.00-1.50", extensions=False, blazecorrect=True, telluric=True):

    
    print "Reading data"
    if extensions:
      self.orders = FitsUtils.MakeXYpoints(infilename, errors="error", extensions=True, x="wavelength", y="flux")
    else:
      self.orders = FitsUtils.MakeXYpoints(infilename, errors=2)

    if blazecorrect:
      print "Reading blaze function"
      #Find and read in blaze function (MUST BE IN CURRENT DIRECTORY!)
      files = os.listdir("./")
      blazefile = [fname for fname in files if fname.startswith("BLAZE")][0]
      blaze_orders = FitsUtils.MakeXYpoints(blazefile, errors=2)
      blaze_functions = []
      blaze_errors = []
      for order in blaze_orders:
        blaze_functions.append( UnivariateSpline(order.x, order.y, s=0) )
      for i, order in enumerate(self.orders):
        self.orders[i].y /= blaze_functions[i](order.x)
        self.orders[i].err /= blaze_functions[i](order.x)
    else:
      blaze_functions = []
      for order in self.orders:
        blaze_functions.append(lambda x: 1.0)
      

    if telluric:
      print "Reading telluric model from database"
      x,y = numpy.loadtxt(telluricfile, unpack=True)
      self.model = DataStructures.xypoint(x=x[::-1], y=y[::-1])
    else:
      x = numpy.arange(self.orders[-1].x[0] - 20.0, self.orders[0].x[-1] + 20.0, 0.001)
      y = numpy.ones(x.size)
      print x
      print y
      self.model = DataStructures.xypoint(x=x, y=y)

    #Make outfilename
    if "-" in infilename:
      num = int(infilename.split("-")[-1].split(".fits")[0])
      outfilename = "%s-%i.fits" %(infilename.split("-")[-1], num+1)
    else:
      outfilename = "%s-0.fits" %(infilename.split(".fits")[0])
    cmdstring = "cp %s %s" %(infilename, outfilename)
    command = subprocess.check_call(cmdstring, shell=True)
    
    self.fitmode = False
    self.clicks = []
    self.template = infilename
    self.blaze_functions = blaze_functions
    self.infilename = infilename
    self.outfilename = outfilename

  def Plot(self):
    start = 2
    for i, order in enumerate(self.orders[start:]):
      self.fig = plt.figure(1, figsize=(11,10))
      plt.title("Order # %i" %(i+start+1))
      plotgrid = gridspec.GridSpec(2,1)
      self.mainaxis = plt.subplot(plotgrid[0])
      self.fitaxis = plt.subplot(plotgrid[1])
      cid = self.fig.canvas.mpl_connect('key_press_event', self.keypress)
      order.cont = FindContinuum.Continuum(order.x, order.y, fitorder=2)
      self.current_order = order.copy()
      left = numpy.searchsorted(self.model.x, order.x[0]-10.0)
      right = numpy.searchsorted(self.model.x, order.x[-1]+10.0)
      current_model = DataStructures.xypoint(x=self.model.x[left:right], y=self.model.y[left:right])
      current_model = MakeModel.ReduceResolution(current_model, 60000)
      self.current_model = MakeModel.RebinData(current_model, order.x)
      offset = self.CCImprove(self.current_order, self.current_model)
      self.current_model.x -= offset
      self.current_model.y *= 0.8

      self.PlotArrays(((order.x, order.y), (self.current_model.x, (self.current_model.y)*self.current_order.cont)), self.mainaxis, legend=False)
      plt.show()

      print "Done with order %i" %(i+start)
      self.orders[i+start] = self.current_order.copy()
      self.orders[i+start].y *= self.blaze_functions[i+start](self.orders[i+start].x)
      self.orders[i+start].err *= self.blaze_functions[i+start](self.orders[i+start].x)

      data = self.orders[i+start]
      columns = {"wavelength": data.x,
                 "flux": data.y,
                 "continuum": data.cont,
                 "error": data.err}
      self.EditFitsFile(columns, self.outfilename, i+start)
      
      #FitsUtils.OutputFitsFile(self.template, self.orders, errors=2)


  def keypress(self, event):
    if event.key == "f":
      #Enter fit mode
      if self.fitmode:
        print "Deactivating fit mode"
        self.fitmode = False
        self.fig.canvas.mpl_disconnect(self.clickid)
        self.clicks = []
        return
      else:
        print "Activating fit mode"
        self.smoothing_factor = 5e-6
        self.fitmode = True
        self.clickid = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        return
      
    if self.fitmode and event.key == "-":
      print "Decreasing smoothing factor: s = %g --> %g" %(self.smoothing_factor, self.smoothing_factor*1.1)
      #before = self.SmoothData()
      self.smoothing_factor *= 1.1
      smoothed = self.SmoothData()
      #print before.y - smoothed.y
      self.fitaxis.cla()
      self.PlotArrays(((self.smoothing_data.x, self.smoothing_data.y/self.smoothing_data.cont, "Data"), (smoothed.x, smoothed.y, "Smoothed")), self.fitaxis)
      
    if self.fitmode and event.key == "+":
      print "Increasing smoothing factor: s = %g --> %g" %(self.smoothing_factor, self.smoothing_factor*0.9)
      self.smoothing_factor *= 0.9
      smoothed = self.SmoothData()
      self.fitaxis.cla()
      self.PlotArrays(((self.smoothing_data.x, self.smoothing_data.y/self.smoothing_data.cont, "Data"), (smoothed.x, smoothed.y, "Smoothed")), self.fitaxis)

    if event.key == ":":
      inp = raw_input("Enter special key (right now just s= some number for the smoothing factor): ")
      if inp.startswith("s"):
        try:
          self.smoothing_factor = float(inp.split("=")[-1])
          smoothed = self.SmoothData()
          self.fitaxis.cla()
          self.PlotArrays(((self.smoothing_data.x, self.smoothing_data.y/self.smoothing_data.cont, "Data"), (smoothed.x, smoothed.y, "Smoothed")), self.fitaxis)
        except ValueError:
          print "Woops! You entered something starting with s, but did not have the format s=0.001 (or some other number after the equals sign)"
          return
      
    if self.fitmode and event.key == "q":
      #Divide data by smoothed version
      smoothed = self.SmoothData()
      self.smoothing_data.y /= smoothed.y
      left = numpy.searchsorted(self.current_order.x, self.smoothing_data.x[0])
      right = numpy.searchsorted(self.current_order.x, self.smoothing_data.x[-1])
      if right < self.current_order.x.size:
        right += 1
      
      self.current_order.y[left:right] = self.smoothing_data.y
      self.mainaxis.cla()
      self.PlotArrays(((self.current_order.x, self.current_order.y), (self.current_model.x, self.current_model.y*self.current_order.cont)), self.mainaxis, legend=False)
      self.fitaxis.cla()
      self.fitmode = False
      

  def onclick(self, event):
    print event.xdata, event.ydata
    self.clicks.append((event.xdata, event.ydata))

    if len(self.clicks) < 2:
      return
    else:
      #Perform fit. Try just fitting splines?
      x1, y1 = self.clicks[0]    #Left-hand continuum
      x2, y2 = self.clicks[1]    #Right-hand continuum
      #x3, y3 = self.clicks[2]    #Line depth
      self.clicks = []
      left = numpy.searchsorted(self.current_order.x, x1)
      right = numpy.searchsorted(self.current_order.x, x2)
      y1 = numpy.median(self.current_order.y[max(0, left-2):min(self.current_order.size(), left+2)])
      y2 = numpy.median(self.current_order.y[max(0, right-2):min(self.current_order.size(), right+2)])
      cont = numpy.poly1d(numpy.polyfit((x1, x2), (y1, y2), 1) )
      self.smoothing_data = DataStructures.xypoint(x=self.current_order.x[left:right],
                                                   y=self.current_order.y[left:right],
                                                   cont=cont(self.current_order.x[left:right]) )
      self.smoothing_factor *= self.smoothing_data.size()
      smoothed = self.SmoothData()
      #smoothed = UnivariateSpline(data.x, data.y/data.cont, s=6e-4 ) #numpy.median(data.y)/10000.0)
      #mean = data.x.mean()
      mean = 0.0
      #smoothed = numpy.poly1d(numpy.polyfit(data.x - mean, data.y/data.cont, 7) )

      self.PlotArrays(((self.smoothing_data.x, self.smoothing_data.y/self.smoothing_data.cont, "Data"), (smoothed.x, smoothed.y, "Smoothed")), self.fitaxis)
      #plt.show()
      
      return

    
  def SmoothData(self, numiters=10, lowreject=2, highreject=10):
    done = False
    data = self.smoothing_data.copy()
    iterations = 0
    while not done and iterations < numiters:
      iterations += 1
      done = True
      smoother = UnivariateSpline(data.x, data.y/data.cont, s=self.smoothing_factor)
      smoothed = DataStructures.xypoint(x=data.x)
      smoothed.y = smoother(smoothed.x)
      resid = data.y/data.cont - smoothed.y
      std = numpy.std(resid)
      badindices = numpy.where(numpy.logical_or(resid < -lowreject*std, resid > highreject*std))[0]
      #plt.figure(2)
      #plt.plot(data.x, resid, 'ro')
      #plt.plot(data.x, -lowreject*std*numpy.ones(data.x.size), 'b-')
      #plt.plot(data.x, highreject*std*numpy.ones(data.x.size), 'b-')
      #plt.show()
      if badindices.size > 0 and data.size() - badindices.size > 10:
        done = False
        data.x = numpy.delete(data.x, badindices)
        data.y = numpy.delete(data.y, badindices)
        data.cont = numpy.delete(data.cont, badindices)
      
    return DataStructures.xypoint(x=self.smoothing_data.x, y=smoother(self.smoothing_data.x))

  def PlotArrays(self, arrays, axis, legend=True):
    axis.cla()
    for arr in arrays:
      if len(arr) > 2:
        axis.plot(arr[0], arr[1], label=arr[2])
      else:
        axis.plot(arr[0], arr[1])
    if legend:
      axis.legend(loc=4)
    plt.draw()


  def CCImprove(self, data, model, be_safe=True, tol=0.5):
    ycorr = numpy.correlate(data.y/data.cont-1.0, model.y-1.0, mode="full")
    xcorr = numpy.arange(ycorr.size)
    maxindex = ycorr.argmax()
    lags = xcorr - (data.y.size-1)
    distancePerLag = (data.x[-1] - data.x[0])/float(data.x.size)
    offsets = -lags*distancePerLag
    print "maximum offset: ", offsets[maxindex], " nm"

    if numpy.abs(offsets[maxindex]) < tol or not be_safe:
      #Apply offset
      print "Applying offset"
      return offsets[maxindex]
    else:
      return 0.0


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
       header_info takes a list of lists. Each sub-list should have size 2 where the first element is the name of the new keyword, and the second element is the corresponding value. A 3rd element may be added as a comment
  """
  def EditFitsFile(self, column_dict, filename, extension, header_info=[]):
    print "Editing extension number %i" %extension
    columns = []
    for key in column_dict.keys():
      columns.append(pyfits.Column(name=key, format="D", array=column_dict[key]))
    cols = pyfits.ColDefs(columns)
    tablehdu = pyfits.new_table(cols)
  
    #Add keywords to extension header
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
  #Parse command line arguments
  files = []
  extensions=False
  blazecorrect=True
  telluric=True
  for arg in sys.argv[1:]:
    if "-e" in arg:
      extensions=True
    elif "-b" in arg:
      blazecorrect=False
    elif "-t" in arg:
      telluric=False
    else:
      files.append(arg)

  #Loop over files
  for fname in files:
    if "linux" in sys.platform:
      fitter = LineFitter(fname, telluricfile="/home/kgullikson/School/Research/aerlbl_v12.2/rundir3/OutputModels/transmission-796.23-270.40-27.1-40.8-368.50-3.90-1.80-1.40", extensions=extensions, blazecorrect=blazecorrect, telluric=telluric)
    else:
      fitter = LineFitter(fname, extensions=extensions, blazecorrect=blazecorrect, telluric=telluric)
    fitter.Plot()