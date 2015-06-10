import sys
import itertools

import matplotlib.pyplot as plt
import numpy as np

import FitsUtils


if __name__ == "__main__":
    fileList = []
    tellurics = False
    normalize = False
    byorder = False  # Plots one order at a time
    pixelscale = False
    oneplot = False
    for arg in sys.argv[1:]:
        if "tellcorr" in arg:
            tellurics = True
        elif "-norm" in arg:
            normalize = True
        elif "-order" in arg:
            byorder = True
        elif "-pix" in arg:
            pixelscale = True
            # byorder = True
        elif "-one" in arg:
            oneplot = True
        else:
            fileList.append(arg)

    # linestyles = ['k-', 'r-', 'b-', 'g-']
    linestyles = itertools.cycle(('-', '--', ':', '-.'))
    colors = itertools.cycle(('r', 'g', 'b', 'c', 'm', 'y', 'k'))

    for fnum, fname in enumerate(fileList):
        #ls = linestyles[fnum % len(linestyles)]
        col = colors.next()
        if fnum % 7 == 0:
            style = linestyles.next()
        ls = '{}{}'.format(col, style)
        orders = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", cont="continuum",
                                        errors="error")
        print fname, len(orders)
        if not oneplot:
            plt.figure(fnum)
            plt.title(fname)
        if tellurics:
            model = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="model")
        for i, order in enumerate(orders):

            # order.cont = FindContinuum.Continuum(order.x, order.y, lowreject=3, highreject=3)
            if pixelscale:
                order.x = np.arange(order.size())
            if tellurics:
                plt.plot(order.x, order.y / order.cont, 'k-')
                plt.plot(order.x, model[i].y, 'r-')
            else:
                if normalize:
                    plt.plot(order.x, order.y / order.cont, ls, label=fname)
                    plt.text(order.x.mean(), 1.1, str(i + 1))
                else:
                    if i == 0:
                        plt.plot(order.x, order.y, ls, label=fname)
                    else:
                        plt.plot(order.x, order.y, ls)
                        #plt.plot(order.x, order.cont)
            plt.xlabel("Wavelength (nm)")
            plt.ylabel("Flux")
            plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
            if byorder:
                plt.title("Order %i" % i)
                plt.show()
    if not byorder:
        if 'oneplot':
            leg = plt.legend(loc='best', fancybox=True)
            leg.get_frame().set_alpha(0.4)
        plt.show()
