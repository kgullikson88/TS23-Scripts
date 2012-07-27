import numpy
import sys
import subprocess
import ConvertMIPASto_lblrtm_format as convert
import scipy.interpolate
import os
import DataStructures

homedir = os.environ['HOME']
TelluricModelingDir = homedir + "/School/Research/lblrtm/run_examples/MyModel/"
if sys.platform.startswith("linux"):
  ModelDir = "/media/FreeAgent_Drive/TelluricLibrary/"
else:
  ModelDir = TelluricModelingDir + "OutputFiles/"


def Main(pressure, temperature, humidity, lowfreq, highfreq, angle, co2, o3, ch4, co, o2, wavegrid, resolution):
    #First, check to make sure the model does not already exist
    print "\n\n\n\n****************************************************"
    print "O2: %.15f" %o2
    print "H20: %.15f " %humidity
    print "P: %.15f" %pressure
    print "T: %.15f" %temperature
    print "Angle: %.15f" %angle
    print "Resolution: %.15f" %resolution
    print "CO: %.15f" %co
    print "CH4: %.15f" %ch4
    print "O3: %.15f" %o3
 
    model_name = ModelDir + "transmission"+"-%.2f" %pressure + "-%.2f" %temperature + "-%.1f" %humidity + "-%.1f" %angle + "-%.2f" %(co2) + "-%.2f" %(o3*100) + "-%.2f" %ch4 + "-%.2f" %(co*10) + "-.2f" %(o2/1e5)
   
    
    #Read in MIPAS_atmosphere_profile
    print "Generating new atmosphere profile"
    filename = TelluricModelingDir + "MIPAS_atmosphere_profile"
    infile = open(filename)
    lines = infile.readlines()
    outputlines = []
    pressureindex = 0
    temperatureindex = 0
    H2Oindex = 0
    CO2index = 0
    COindex = 0
    CH4index = 0
    O3index = 0
    O2index = 0
    for i in range(len(lines)):
        line = lines[i]
        if line.find("HGT") > 0 and line.find("[") > 0:
            numlevels = int(lines[i-1].split()[0])
        elif line.find("PRE") > 0 and line.find("[") > 0:
            pressureindex = i
        elif line.find("TEM") > 0 and line.find("[") > 0:
            temperatureindex = i
        elif line.find("O2 ") == 1 and line.find("[") > 0:
            O2index = i
        elif line.find("H2O ") == 1 and line.find("[") > 0:
            H2Oindex = i
        elif line.find("CO2 ") == 1 and line.find("[") > 0:
            CO2index = i
        elif line.find("CO ") == 1 and line.find("[") > 0:
            COindex = i
        elif line.find("CH4 ") == 1 and line.find("[") > 0:
            CH4index = i
        elif line.find("O3 ") == 1 and line.find("[") > 0:
            O3index = i
    infile.close()
    #Determine the number of lines that follow each header
    levelsperline = 5.0
    linespersection = int(numlevels/levelsperline + 0.9)
    #assume the order in which parameters appear is:
    #pressure, temperature, O2, CO2, O3, H20, CH4, CO
    for line in lines[:pressureindex+1]:
        outputlines.append(line)
    #Loop over the pressure indices
    oldpressure = float(lines[pressureindex+1].split()[2])
    scalefactor = pressure/oldpressure
    for i in range(pressureindex+1,pressureindex+1+linespersection):
        line = lines[i]
        levels = line.split()
        newline = ""
        for level in levels:
            level = float(level)*scalefactor
            newline = newline + " " + "%.5E" %level
        newline = newline + "\n"
        outputlines.append(newline)
    
    #Loop over everything from pressure to temperature
    for i in range(pressureindex+1+linespersection, temperatureindex+1):
        outputlines.append(lines[i])

    #Loop over temperature indices
    oldtemperature = float(lines[temperatureindex+1].split()[2])
    scalefactor = temperature/oldtemperature
    for i in range(temperatureindex+1, temperatureindex+1+linespersection):
        line = lines[i]
        levels = line.split()
        newline = ""
        for level in levels:
            level = float(level)*scalefactor
            newline = newline + " " + "%.2f" %level
        newline = newline + "\n"
        outputlines.append(newline)
        
    #Loop over everything from temperature to O2
    for i in range(temperatureindex+1+linespersection, O2index+1):
        outputlines.append(lines[i])

    #Loop over O2 indices
    oldo2 = float(lines[O2index+1].split()[2])
    scalefactor = o2/oldo2
    for i in range(O2index+1, O2index+1+linespersection):
        line = lines[i]
        levels = line.split()
        newline = ""
        for level in levels:
            level = float(level)*scalefactor
            newline = newline + " " + "%.3e" %level
        newline = newline + "\n"
        outputlines.append(newline)
  
    #loop over everything from O2 to CO2:
    for i in range (O2index+1+linespersection, CO2index+1):
        outputlines.append(lines[i])   

    #Loop over CO2 indices
    oldno2 = float(lines[CO2index+1].split()[2])
    scalefactor = co2/oldno2
    for i in range(CO2index+1, CO2index+1+linespersection):
        line = lines[i]
        levels = line.split()
        newline = ""
        for level in levels:
            level = float(level)*scalefactor
            newline = newline + " " + "%.3e" %level
        newline = newline + "\n"
        outputlines.append(newline)
  
    #loop over everything from CO2 to O3:
    for i in range (CO2index+1+linespersection, O3index+1):
        outputlines.append(lines[i])

    #Loop over O3 indices
    oldo3 = float(lines[O3index+1].split()[2])
    scalefactor = o3/oldo3
    for i in range(O3index+1, O3index+1+linespersection):
        line = lines[i]
        levels = line.split()
        newline = ""
        for level in levels:
            level = float(level)*scalefactor
            newline = newline + " " + "%.3e" %level
        newline = newline + "\n"
        outputlines.append(newline)

    #loop over everything from O3 to H20:
    for i in range (O3index+1+linespersection, H2Oindex+1):
        outputlines.append(lines[i])

    #Convert from relative humidity to concentration (ppm)
    #formulas and constants come from http://www.scribd.com/doc/53758450/Humidity-Conversion-Formulas-B210973EN-B-Lores
    Psat = 6.1162*10**(7.5892*(temperature-273.15)/(240.71+(temperature-273.15)))
    Pw = Psat*humidity/100.0
    concentration = Pw/(pressure - Pw)*1e6
    
    #Loop over H2O indices
    oldconcentration = float(lines[H2Oindex+1].split()[2])
    scalefactor = concentration/oldconcentration
    for i in range(H2Oindex+1, H2Oindex+1+linespersection):
        line = lines[i]
        levels = line.split()
        newline = ""
        for level in levels:
            level = float(level)*scalefactor
            newline = newline + " " + "%.3e" %level
        newline = newline + "\n"
        outputlines.append(newline)

    #loop over everything from H20 to CH4:
    for i in range (H2Oindex+1+linespersection, CH4index+1):
        outputlines.append(lines[i])

    #Loop over CH4 indices
    oldch4 = float(lines[CH4index+1].split()[2])
    scalefactor = ch4/oldch4
    index = 0
    for i in range(CH4index+1, CH4index+1+linespersection):
        line = lines[i]
        levels = line.split()
        newline = ""
        for i in range(len(levels)):
            #level = CH4_atmosphere_profile(ch4, breakpoint, decay, float(index))
            #print index, "\t", level
            level = float(levels[i])*scalefactor
            newline = newline + " " + "%.3e" %level
            index = index + 1
        newline = newline + "\n"
        outputlines.append(newline)
        
    #loop over everything from CH4 to CO:
    for i in range (CH4index+1+linespersection, COindex+1):
        outputlines.append(lines[i])  

    #Loop over CO indices
    oldno = float(lines[COindex+1].split()[2])
    scalefactor = co/oldno
    for i in range(COindex+1, COindex+1+linespersection):
        line = lines[i]
        levels = line.split()
        newline = ""
        for level in levels:
            level = float(level)*scalefactor
            newline = newline + " " + "%.3e" %level
        newline = newline + "\n"
        outputlines.append(newline)

    #Loop over the rest
    for i in range(COindex + 1 + linespersection, len(lines)):
        outputlines.append(lines[i])
        
    
    ###########################################################
    #  We are done generating the new atmosphere. Now, 
    #  reformat it so that lblrtm can use it
    ###########################################################
    lines = convert.Main(outputlines)
    
    
    #Read in the new file and ParameterFile, and replace the atmosphere section of ParameterFile
    #Note: this method assumes that ParameterFile has the same number of atmosphere layers as the new atmosphere. May not necessarily be true!
    infile = open(TelluricModelingDir + "ParameterFile")
    outputlines = infile.readlines()
    infile.close()
    #infile = open(outfilename + ".reformatted")
    #lines = infile.readlines()
    #infile.close()
    firstindex = 0
    lastindex = 0
    for i in range(len(outputlines)):
        line = outputlines[i]
        if line.find("Record3.5") > 0:
            firstindex = i
        elif line.find("Record3.2H") > 0:
            lastindex=i
        elif line.find("ANGLE") != -1:
            #change angle to user-given value
            string = line.split()[0] + "    " + str(angle) + "\n"
            outputlines[i] = string
 
    for i in range(firstindex+1, lastindex):
        outputlines[i] = lines[i-(firstindex+1)] + "\n"

    if (highfreq - lowfreq > 2000.0):
        while lowfreq + 2000.0 <= highfreq:
            print lowfreq, "\t", lowfreq+2000
            for i in range(len(outputlines)):
                line = outputlines[i]
                if line.find("V1....") != -1:
                    #change frequency to lowfreq
                    string = line.split()[0] + "    %.2f" %(lowfreq) + "\n"
                    outputlines[i] = string
                elif line.find("V2...") != -1:
                    #change frequency to highfreq
                    string = line.split()[0] + "    %.2f" %(lowfreq + 2000.0) + "\n"
                    outputlines[i] = string
            lowfreq = lowfreq + 2000.0
            
             #Output new parameter file
            outfile = open(TelluricModelingDir + "ParameterFile", "w")
            for line in outputlines:
                outfile.write(line)
            outfile.close()
            
            #Run lblrtm
            cmd = "cd " + TelluricModelingDir + ";sh runlblrtm.sh"
            command = subprocess.check_call(cmd, shell=True)
    else:
        for i in range(len(outputlines)):
            line = outputlines[i]
            if line.find("V1....") != -1:
                #change frequency to lowfreq
                string = line.split()[0] + "    %.2f" %(lowfreq) + "\n"
                outputlines[i] = string
            elif line.find("V2...") != -1:
                #change frequency to highfreq
                string = line.split()[0] + "    %.2f" %(highfreq) + "\n"
                outputlines[i] = string
    
        #Output new parameter file
        outfile = open(TelluricModelingDir + "ParameterFile", "w")
        for line in outputlines:
            outfile.write(line)
        outfile.close()

        #Run lblrtm
        cmd = "cd " + TelluricModelingDir + ";sh runlblrtm.sh"
        command = subprocess.check_call(cmd, shell=True)

    #Convert from frequency to wavelength units
    wavelength, transmission = FixTelluric(TelluricModelingDir + "FullSpectrum.freq")
    
    if "FullSpectrum.freq" in os.listdir(TelluricModelingDir):
      cmd = "rm " + TelluricModelingDir + "FullSpectrum.freq"
      command = subprocess.check_call(cmd, shell=True)
    
    #Append model name to library_list
    outfile = open("library_list", "a")
    outfile.write(model_name+"\n")
    outfile.close()
    
    #Interpolate model to a constant wavelength grid
    #For now, just interpolate to the right spacing
    #Also, only keep the section near the chip wavelengths
    wavelength = wavelength[::-1]
    transmission = transmission[::-1]
    xspacing = (wavelength[-1] - wavelength[0])/float(wavelength.size)
    tol = 10  #Go 10 nm on either side of the chip
    left = numpy.searchsorted(wavelength, wavegrid[0]-tol)
    right = numpy.searchsorted(wavelength, wavegrid[-1]+tol) 
      
    Model = scipy.interpolate.UnivariateSpline(wavelength, transmission, s=0)
    model = DataStructures.xypoint(right-left+1)
    model.x = numpy.arange(wavelength[left], wavelength[right], xspacing)
    model.y = Model(model.x)
    
    print "All done! Output Transmission spectrum is located in the file below:"
    print model_name
    return model
    
  


def FixTelluric(filename):
  wavenumber, transmission = numpy.loadtxt(filename,unpack=True)
  wavelength = 1e4/wavenumber
  outfile = open(TelluricModelingDir + "FullSpectrum.wave", "w")
  for i in range(wavelength.size):
    outfile.write(str(wavelength[i]) + "\t" + str(transmission[i]) + "\n")
  return wavelength*1000.0, transmission

#This function rebins (x,y) data onto the grid given by the array xgrid
def RebinData(data,xgrid):
  Model = scipy.interpolate.UnivariateSpline(data.x, data.y, s=0)
  newdata = DataStructures.xypoint(xgrid.size)
  newdata.x = numpy.copy(xgrid)
  newdata.y = Model(newdata.x)
  
  left = numpy.searchsorted(data.x, (3*xgrid[0]-xgrid[1])/2.0)
  for i in range(xgrid.size-1):
    right = numpy.searchsorted(data.x, (xgrid[i]+xgrid[i+1])/2.0)
    newdata.y[i] = numpy.mean(data.y[left:right])
    left = right
  right = numpy.searchsorted(data.x, (3*xgrid[-1]-xgrid[-2])/2.0)
  newdata.y[xgrid.size-1] = numpy.mean(data.y[left:right])
  
  return newdata

#This function reduces the resolution by convolving with a gaussian kernel
def ReduceResolution(data,resolution, cont_fcn=None, extend=True):
  centralwavelength = (data.x[0] + data.x[-1])/2.0
  xspacing = data.x[1] - data.x[0]   #NOTE: this assumes constant x spacing!
  FWHM = centralwavelength/resolution;
  sigma = FWHM/(2.0*numpy.sqrt(2.0*numpy.log(2.0)))
  left = 0
  right = numpy.searchsorted(data.x, 10*sigma)
  x = numpy.arange(0,10*sigma, xspacing)
  gaussian = numpy.exp(-(x-5*sigma)**2/(2*sigma**2))
  if extend:
    #Extend array to try to remove edge effects (do so circularly)
    before = data.y[-gaussian.size/2+1:]
    after = data.y[:gaussian.size/2]
    extended = numpy.append(numpy.append(before, data.y), after)

    first = data.x[0] - float(int(gaussian.size/2.0+0.5))*xspacing
    last = data.x[-1] + float(int(gaussian.size/2.0+0.5))*xspacing
    x2 = numpy.linspace(first, last, extended.size) 
    
    conv_mode = "valid"

  else:
    extended = data.y.copy()
    x2 = data.x.copy()
    conv_mode = "same"

  newdata = DataStructures.xypoint(data.x.size)
  newdata.x = numpy.copy(data.x)
  if cont_fcn != None:
    cont1 = cont_fcn(newdata.x)
    cont2 = cont_fcn(x2)
    cont1[cont1 < 0.01] = 1
  
    newdata.y = numpy.convolve(extended*cont2, gaussian/gaussian.sum(), mode=conv_mode)/cont1

  else:
    newdata.y = numpy.convolve(extended, gaussian/gaussian.sum(), mode=conv_mode)
    
  return newdata

#Just a convenince fcn which combines the above two
def ReduceResolutionAndRebinData(data,resolution,xgrid):
  data = ReduceResolution(data,resolution)
  return RebinData(data,xgrid)
  
  
def CH4_atmosphere_profile(ch4, breakpoint, decay, z):
  if z < breakpoint:
    return ch4
  else:
    return ch4*numpy.exp(-(z-breakpoint)**2/(2*decay**2))
  

"""
if __name__ == "__main__":
    if len(sys.argv) > 10:
        pressure = float(sys.argv[1])
        temperature = float(sys.argv[2])
        humidity = float(sys.argv[3])
        lowfreq = float(sys.argv[4])
        highfreq = float(sys.argv[5])
        angle = float(sys.argv[6])
        co2 = float(sys.argv[7])
        o3 = float(sys.argv[8])
        ch4 = float(sys.argv[9])
        co = float(sys.argv[10])
        Main(pressure, temperature, humidity, lowfreq, highfreq, angle, co2, o3, ch4, co)
        #Main(pressure, temperature, humidity, lowfreq, highfreq, angle)
        if "FullSpectrum.freq" in os.listdir("./"):
          cmd = "rm FullSpectrum.freq"
          command = subprocess.check_call(cmd, shell=True)
          
          
"""
