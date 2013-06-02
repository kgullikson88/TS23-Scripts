import numpy
import matplotlib.pyplot as plt
import sys
import os
from collections import defaultdict

Corr_dir = "./Cross_correlations/"


if __name__ == "__main__":
  basename = sys.argv[1]
  velocity = float(sys.argv[2])  #Turn this into an optional argument later...
  
  files = os.listdir(Corr_dir)
  fileList = defaultdict(list)
  for fname in files:
    if fname.startswith(basename) and "10kps" in fname:
      vsini = float(fname.split("kps")[0].split(".")[-1])
      fileList[vsini].append(fname)
  for val in sorted(fileList.keys()):
    files = fileList[val]
    #Temperatures = []
    #Significances = []
    Temperatures = defaultdict(list)
    Significances = defaultdict(list)
    for fname in files:
      vels, corr = numpy.loadtxt(Corr_dir + fname, unpack=True)
      fit = numpy.poly1d( numpy.polyfit(vels, corr, 2) )
      std = numpy.std(corr - fit(vels))
      index = numpy.searchsorted(vels, velocity)
      metallicity = fname[-4:]
      #Temperatures.append(float(fname.split("_")[-1].split("K")[0]))
      Temperatures[metallicity].append(float(fname.split("_")[-1].split("K")[0]))
      #Significances.append((corr[index] - fit(vels[index]))/std)
      Significances[metallicity].append((corr[index] - fit(vels[index]))/std)
    Temperatures2 = defaultdict(numpy.array)
    sorter = defaultdict(list)
    for metal in Temperatures:
      Temperatures2[metal] = numpy.array(Temperatures[metal])
      sorter[metal] = sorted(range(len(Temperatures2[metal])),key=lambda x:Temperatures2[metal][x])
    Significances2 = defaultdict(numpy.array)
    for metal in Significances:
      Significances2[metal] = numpy.array(Significances[metal])
    #Temperatures = numpy.array(Temperatures)
    #Significances = numpy.array(Significances)
    #sorter = sorted(range(len(Temperatures)),key=lambda x:Temperatures[x])
    #T = Temperatures[sorter]
    #S = Significances[sorter]
    for Z in Significances2.keys():
      plt.plot(Temperatures2[Z][sorter[Z]], Significances2[Z][sorter[Z]], 'o', label="[Fe/H] = %s" %Z)
    #plt.plot(T, S, '-o', color='black')
    plt.xlabel("Secondary Temperature (K)")
    plt.ylabel("CCF Value at v = %g" %velocity)
    plt.legend(loc='best')
    plt.title("vsini = %g" %val)
    plt.show()
    
