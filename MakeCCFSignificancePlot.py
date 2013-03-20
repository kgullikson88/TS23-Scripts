import numpy
import matplotlib.pyplot as plt
import sys
import os
from collections import defaultdict

Corr_dir = "Cross_Correlations/"


if __name__ == "__main__":
  basename = sys.argv[1]
  velocity = float(sys.argv[2])  #Turn this into an optional argument later...
  
  files = os.listdir(Corr_dir)
  fileList = defaultdict(list)
  for fname in files:
    if fname.startswith(basename):
      vsini = float(fname.split("kps")[0].split(".")[-1])
      fileList[vsini].append(fname)
  for val in sorted(fileList.keys()):
    files = fileList[val]
    Temperatures = []
    Significances = []
    for fname in files:
      vels, corr = numpy.loadtxt(Corr_dir + fname, unpack=True)
      index = numpy.searchsorted(vels, velocity)
      Temperatures.append(float(fname.split("_")[-1][:-1]))
      Significances.append(corr[index])
    Temperatures = numpy.array(Temperatures)
    Significances = numpy.array(Significances)
    sorter = sorted(range(len(Temperatures)),key=lambda x:Temperatures[x])
    T = Temperatures[sorter]
    S = Significances[sorter]
    plt.plot(T, S, '-o')
    plt.xlabel("Secondary Temperature (K)")
    plt.ylabel("CCF Value at v = %g" %velocity)
    plt.title("vsini = %g" %val)
    plt.show()
    
