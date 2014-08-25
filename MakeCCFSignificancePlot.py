import sys
import os
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt


Corr_dir = "./Cross_correlations/"

if __name__ == "__main__":
    basename = sys.argv[1]
    velocity = float(sys.argv[2])  # Turn this into an optional argument later...

    files = os.listdir(Corr_dir)
    fileList = defaultdict(list)
    for fname in files:
        if fname.startswith(basename) and "10kps" in fname:
            vsini = float(fname.split("kps")[0].split(".")[-1])
            fileList[vsini].append(fname)
    for val in sorted(fileList.keys()):
        files = fileList[val]
        # Temperatures = []
        #Significances = []
        Temperatures = defaultdict(list)
        Significances = defaultdict(list)
        for fname in files:
            vels, corr = np.loadtxt(Corr_dir + fname, unpack=True)
            fit = np.poly1d(np.polyfit(vels, corr, 2))
            std = np.std(corr - fit(vels))
            index = np.searchsorted(vels, velocity)
            metallicity = fname[-4:]
            #Temperatures.append(float(fname.split("_")[-1].split("K")[0]))
            Temperatures[metallicity].append(float(fname.split("_")[-1].split("K")[0]))
            #Significances.append((corr[index] - fit(vels[index]))/std)
            Significances[metallicity].append((corr[index] - fit(vels[index])) / std)
        Temperatures2 = defaultdict(np.array)
        sorter = defaultdict(list)
        for metal in Temperatures:
            Temperatures2[metal] = np.array(Temperatures[metal])
            sorter[metal] = sorted(range(len(Temperatures2[metal])), key=lambda x: Temperatures2[metal][x])
        Significances2 = defaultdict(np.array)
        for metal in Significances:
            Significances2[metal] = np.array(Significances[metal])
        #Temperatures = np.array(Temperatures)
        #Significances = np.array(Significances)
        #sorter = sorted(range(len(Temperatures)),key=lambda x:Temperatures[x])
        #T = Temperatures[sorter]
        #S = Significances[sorter]
        for Z in Significances2.keys():
            plt.plot(Temperatures2[Z][sorter[Z]], Significances2[Z][sorter[Z]], 'o', label="[Fe/H] = %s" % Z)
        #plt.plot(T, S, '-o', color='black')
        plt.xlabel("Secondary Temperature (K)")
        plt.ylabel("CCF Value at v = %g" % velocity)
        plt.legend(loc='best')
        plt.title("vsini = %g" % val)
        plt.show()
    
