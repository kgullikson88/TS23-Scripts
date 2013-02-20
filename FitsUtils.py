import pyfits
import numpy
from numpy.polynomial import chebyshev
import DataStructures
import Units
import FindContinuum
import readmultispec as multispec
import subprocess


#Make a dictionary for converting from standard polynomial to chebyshev
#   For now, only go up to order 5
#standard_to_cheb = {0: [[1,],],
#                    1: [[1,0], [0,1]],
#                    2: [[1,0,0.5], [0,1,0], [0,0,0.5]],
#                    3: [[1,0,0.5,0], [0,1,0,0.75], [0,0,0.5,0], [0,0,0,0.25]], 
#                    4: [[1,0,0.5,0,0.625], [0,1,0,0.75,0], [0,0,0.5,0,-0.5], [0,0,0,0.25,0], [0,0,0,0,0.125]], 
#                    5: [[1,0,0.5,0,0.625, 0], [0,1,0,0.75,0,0.625], [0,0,0.5,0,-0.5,0], [0,0,0,0.25,0,5./16.], [0,0,0,0,0.125,0], [0,0,0,0,0,1./16.]]}

def Chebyshev(pvals, coefficients):
  order = coefficients.size
  pmid = pvals[pvals.size/2]
  prange = pvals[-1] - pvals[0]
  nvals = (pvals - pmid)/(prange/2.0)
  wave = numpy.zeros(pvals.size)
  x = []
  x.append(numpy.ones(pvals.size))
  x.append(nvals)
  for i in range(2,order):
    x.append(2.0*nvals*x[i-1] - x[i-2])
  for i in range(order):
    wave = wave + coefficients[i]*x[i]
  return wave
  
  
  
def GetChebyshevCoeffs(data, pvals, order=5):
  pmid = pvals[pvals.size/2]
  prange = pvals[-1] - pvals[0]
  nvals = (pvals - pmid)/(prange/2.0)
  
  return chebyshev.chebfit(nvals, data.x, order)
  
"""

The following was my attempt. It does not work quite right. Now I use Rick White's implementation (in .PythonModules/GeneralScripts as readmultispec.py)

#Takes a fits file, with nonlinear wavelength calibration given in the header,
#  and generates a DataStructures.xypoint object for each spectral order
#TODO: make it work for functions other than chebyshev, and for linear/log-linear spacing
#      allow for multiple functions. Right now it only works for one type, while fits allows many
def MakeXYpoints(header, data=None):
  #Header can be the filename, too!
  if type(header) == str:
    hdulist = pyfits.open(header)
    header = hdulist[0].header
    data = hdulist[0].data
    hdulist.close()
  
  try:
    ltv = header['ltv1']
  except KeyError:
    ltv = 0.0
  try:
    ltm = header['ltm1_1']
  except KeyError:
    ltm = 1.0

  #Make list of the wavelength data for each order
  string = ""
  wave_factor = 1.0   #factor to multiply wavelengths by to get them in nanometers
  for key in sorted(header.keys()):
    if "WAT2" in key:
      string = string + header[key]
      #print key, header[key]
      for i in range(len(header[key]), 68):
        string = string + " "
    elif "WAT1" in key:
      #Get wavelength units
      if "label=Wavelength"  in header[key] and "units" in header[key]:
        units = header[key].split("units=")[-1]
        if units == "angstroms" or units == "Angstroms":
          wave_factor = Units.nm/Units.angstrom
          print "Wavelength units are Angstroms. Scaling wavelength by ", wave_factor
  orders = string.split(" spec")

  #Make a list of xypoints called DATA
  DATA = []
  index = 0
  for order in orders[1:]:
    #print order
    order = order.split("=")[-1].split('"')[1]
    segments = order.split()
    size = int(float(segments[5]))
    z = float(segments[6])  #Dopplar correction: 1/(1+z)
    
    xypt = DataStructures.xypoint(size)
    pvals = (numpy.arange(size) - ltv)/ltm

    wlen0 = float(segments[10])
    func_type = int(float(segments[11]))
    if func_type == 1:
      #Chebyshev
      func_order = int(float(segments[12]))
      pmin = int(float(segments[13]))
      pmax = int(float(segments[14]))
      coefficients = []
      for segment in segments[15:]:
        coefficients.append(float(segment))
      xypt.x = (wlen0 + Chebyshev(pvals, numpy.array(coefficients)))/(1.0+z)*wave_factor
      xypt.y = data[index]
      xypt.cont = numpy.ones(xypt.x.size)
      xypt.err = numpy.ones(xypt.size())*1e9
      xypt.err[xypt.y > 0] = numpy.sqrt(xypt.y[xypt.y > 0])
    else:
      print "Sorry! This function type is not currently implemented!"

    DATA.append(xypt)
    index += 1

  return DATA
  
"""



"""
  Function to read in a multispec file and return a set of orders (a list of xypoints)
  The 'errors' keyword can be set to an integer to give the band index of the errors
"""
def MakeXYpoints(datafile, errors=False):
  #Call Rick White's script
  retdict = multispec.readmultispec(datafile)
  
  #Check if wavelength units are in angstroms (common, but I like nm)
  hdulist = pyfits.open(datafile)
  header = hdulist[0].header
  hdulist.close()
  wave_factor = 1.0   #factor to multiply wavelengths by to get them in nanometers
  for key in sorted(header.keys()):
    if "WAT1" in key:
      if "label=Wavelength"  in header[key] and "units" in header[key]:
        units = header[key].split("units=")[-1]
        if units == "angstroms" or units == "Angstroms":
          wave_factor = Units.nm/Units.angstrom
          print "Wavelength units are Angstroms. Scaling wavelength by ", wave_factor

  if errors == False:
    numorders = retdict['flux'].shape[0]
  else:
    numorders = retdict['flux'].shape[1]
  orders = []
  for i in range(numorders):
    wave = retdict['wavelen'][i]*wave_factor
    if errors == False:
      flux = retdict['flux'][i]
      err = numpy.ones(flux.size)*1e9
      err[flux > 0] = numpy.sqrt(flux[flux > 0])
    else:
      flux = retdict['flux'][0][i]
      err = retdict['flux'][errors][i]
    cont = FindContinuum.Continuum(wave, flux, lowreject=2, highreject=4)
    orders.append(DataStructures.xypoint(x=wave, y=flux, err=err , cont=cont))
  return orders
  
  
#Function to output a fits file with the same format
#as the template function. Only implementing Chebyshev for now...
def OutputFitsFile(template, orders, func_order=None, outfilename=None):
  hdulist = pyfits.open(template)
  header = hdulist[0].header
  
  #First, get what the current wavelength calibration is from the header
  try:
    ltv = header['ltv1']
  except KeyError:
    ltv = 0.0
  try:
    ltm = header['ltm1_1']
  except KeyError:
    ltm = 1.0

  #Make list of the wavelength data for each order
  string = ""
  wave_factor = 1.0   #factor to multiply wavelengths by to get them in nanometers
  for key in sorted(header.keys()):
    if "WAT2" in key:
      string = string + header[key]
      #print key, header[key]
      for i in range(len(header[key]), 68):
        string = string + " "
    elif "WAT1" in key:
      #Get wavelength units
      if "label=Wavelength"  in header[key] and "units" in header[key]:
        units = header[key].split("units=")[-1]
        if units == "angstroms" or units == "Angstroms":
          header.update(key,"wtype=multispec label=Wavelength units=nanometers")
        
  waveinfo = string.split(" spec")
  output_string = waveinfo[0]
  for info, i in zip(waveinfo[1:], range(1,len(waveinfo))):
    output_string = output_string + " spec" + info.split('"')[0] + '"'
    segments = info.split('"')[1].split()
    for segment in segments[:10]:
      output_string = output_string + segment + " "
    size = int(float(segments[5]))
    z = float(segments[6])  #Dopplar correction: 1/(1+z)
    
    xypt = DataStructures.xypoint(size)
    pvals = (numpy.arange(size) - ltv)/ltm

    wlen0 = 0.0
    func_type = int(float(segments[11]))
    if func_type == 1:
      #Chebyshev
      if func_order == None:
        func_order = int(float(segments[12]))
      pmin = int(float(segments[13]))
      pmax = int(float(segments[14]))
      coefficients = GetChebyshevCoeffs(orders[i-1], pvals, func_order-1)

      #Put stuff in output string
      output_string = output_string + ("%.10g " %wlen0).replace("e","E")
      output_string = output_string + segments[11] + " "
      output_string = output_string + str(func_order) + " "
      output_string = output_string + segments[13] + " "
      output_string = output_string + segments[14] + " "
      for coeff in coefficients:
        output_string = output_string + ("%.14g " %coeff).replace("e","E")
     
    else:
      print "Sorry! This function type is not currently implemented!"

    output_string = output_string[:-1] + '"'
  
  #Split output string so it isn't too many characters
  size = 68
  leftindex=0
  rightindex = size
  keyindex = 1
  while rightindex < len(output_string):
    key = "WAT2_"
    if len(str(keyindex)) == 1:
      key = key + "00"
    elif len(str(keyindex)) == 2:
      key = key + "0"
    key = key+str(keyindex)
    
    header.update(key, output_string[leftindex:rightindex])
    leftindex = rightindex
    rightindex = min(leftindex+size, len(output_string))
    keyindex += 1
  
  #Do the last one outside of the loop
  key = "WAT2_"
  if len(str(keyindex)) == 1:
    key = key + "00"
  elif len(str(keyindex)) == 2:
    key = key + "0"
  key = key+str(keyindex)
  #print key, output_string[leftindex:rightindex]
  header.update(key, output_string[leftindex:rightindex])
  
  #Delete any WAT2 keys that were not used
  for key in sorted(header.keys()):
    if "WAT2" in key:
      if int(key[-3:]) > keyindex:
        del header[key]
    
  hdulist[0].header = header

  #Now, update the data (much easier)
  for i in range(hdulist[0].data.shape[0]):
    hdulist[0].data[i] = orders[i].y

  if outfilename == None:
    if "-" in template:
      i = int(template.split("-")[-1].split(".fits")[0])
      outfilename = template.split("-")[0] + "-" + str(i+1) + ".fits"
    else:
      outfilename = template.split(".fits")[0] + "-0.fits"
  
  print "Outputting to ", outfilename
  try:
    hdulist.writeto(outfilename)
  except IOError:
    cmd = subprocess.check_call("rm "+outfilename, shell=True)
    hdulist.writeto(outfilename)
  hdulist.close()
  
  
def GetOrderWaveCalPars(header):
  string = ""
  for key in sorted(header.keys()):
    if "WAT2" in key:
      string = string + header[key]
      #print key, header[key]
      for i in range(len(header[key]), 68):
        string = string + " "
        
  return string.split(" spec")
  
  
def CopyWaveCal(copyfrom, copyto, order=None, scale=1.0):
  """
    This function will copy the wavelength calibration of one or more orders and put it in another file. 
    copyfrom and copyto should both be filenames
    The order keyword can be:
       1: left as None, which will copy all of the orders
       2: an integer, which will only copy that order (number in fortran style!)
       3: a list of integers, in which the calibration will be copied for each one
     scale is a factor to apply to the wavelengths (useful for converting between units)
  """
  
  #First, read in the files
  hdulistfrom = pyfits.open(copyfrom)
  hdulistto = pyfits.open(copyto)
  headerfrom = hdulistfrom[0].header
  headerto = hdulistto[0].header
  
  numorders = hdulistfrom[0].data.shape[0]
  if numorders != hdulistto[0].data.shape[0]:
    print "Warning in CopyWaveCal! Files have different number of orders!"
    numorders = min(numorders, hdulistto[0].data.shape[0])
  
  if type(order) == int:
    order = [order,]
  elif order == None:
    order = range(1,numorders+1)
    
  #Make sure we don't try to do too many orders
  order = sorted(order)
  while len(order) > 0 and order[-1] > numorders:
    order.pop()
  
  #Now, get all the orders from both files
  ordersfrom = GetOrderWaveCalPars(headerfrom)
  ordersto = GetOrderWaveCalPars(headerto)

  
  #Do the copy
  for i in order:
    print "Copying order ", i
    print numorders
    ordersto[i] = ordersfrom[i]
    
    if scale != 1.0:
      #Now we have to get more complicated... sigh
      #Only tested with chebyshev, should work with legendre. Maybe others?
      segments = ordersto[i].split()
      
      #Coefficients start in segment 16
      for j in range(17,len(segments)):
        coeff = segments[j]
        print coeff
        coeff = "%.14g" %(float(coeff.replace("E","e").strip('"'))*scale)
        print coeff
        segments[j] = coeff
        
      #re-save as a string
      ordersto[i] = segments[0]
      for segment in segments[1:]:
        ordersto[i] = ordersto[i] + " " + segment
      ordersto[i] = ordersto[i] + '"'
  
  #Convert the orders into one long string
  output_string = ordersfrom[0]
  for segment in ordersto[1:]:
    output_string = output_string + " spec" + segment
  
  #Update the header
  size = 68
  leftindex=0
  rightindex = size
  keyindex = 1
  while rightindex < len(output_string):
    key = "WAT2_"
    if len(str(keyindex)) == 1:
      key = key + "00"
    elif len(str(keyindex)) == 2:
      key = key + "0"
    key = key+str(keyindex)
    
    print key, output_string[leftindex:rightindex]
    headerto.update(key, output_string[leftindex:rightindex])
    leftindex = rightindex
    rightindex = min(leftindex+size, len(output_string))
    keyindex += 1
  
  #Do the last one outside of the loop
  key = "WAT2_"
  if len(str(keyindex)) == 1:
    key = key + "00"
  elif len(str(keyindex)) == 2:
    key = key + "0"
  key = key+str(keyindex)
  print key, output_string[leftindex:rightindex]
  headerto.update(key, output_string[leftindex:rightindex])
  
  #Delete any WAT2 keys that were not used
  for key in sorted(headerto.keys()):
    if "WAT2" in key:
      if int(key[-3:]) > keyindex:
        del headerto[key]
    
  hdulistto[0].header = headerto
  
  #Save with new filename
  if "-" in copyto:
    i = int(copyto.split("-")[-1].split(".fits")[0])
    outfilename = copyto.split("-")[0] + "-" + str(i+1) + ".fits"
  else:
    outfilename = copyto.split(".fits")[0] + "-0.fits"
  print "Outputting to ", outfilename
  try:
    hdulistto.writeto(outfilename)
  except IOError:
    cmd = subprocess.check_call("rm "+outfilename, shell=True)
    hdulistto.writeto(outfilename)
  hdulistto.close()
  hdulistfrom.close()
  
  return outfilename
  
  
  
  
if __name__ == "__main__":
  hdulist = pyfits.open("output-1.fits")
  orders = MakeXYpoints(hdulist[0].header, hdulist[0].data)
  OutputFitsFile("output-1.fits", orders)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
