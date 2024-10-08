import numpy as np
import LT.box as B
# Random Number Generator
import scipy.special as sp
# Parameter Module
import LT_Fit.parameters as P
# Genfit module.
import LT_Fit.gen_fit as G
import matplotlib.pyplot as plt


def set_range(first_bin = 0, bin_width = 1., Nbins = 5):

    rmin = first_bin - bin_width/2.0
    rmax = rmin + Nbins*bin_width
    return (rmin, rmax)


fileName = 'PartB.1sec.txt'
f = B.get_file(fileName)
Counts = B.get_data(f, 'C')
CountsMean = Counts.mean()
CountsStd = Counts.std()


h = B.histo(Counts, range = set_range(first_bin = 0, bin_width = 1.,
Nbins =5), bins =5, title = 'Counts for 1 event per Time Interval', xlabel = 'Number of Counts per Interval',
ylabel = 'Frequency of Count' )


# Putting the Bin Centers and Bin Contents into two Arrays
hx = h.bin_center
hy = h.bin_content
dy = np.sqrt(hy) # Error of Bin Contents


# Poisson Function Fitting
muP = P.Parameter(1.0, 'muP') # Initial Guess for mu
normP = P.Parameter(10.0, 'normP') # Initial Guess for normalization factor
def Poisson(x):
    value = normP()*muP()**(x)*np.exp(-muP())
    return value
fitPoisson = G.genfit(Poisson, [muP, normP], x = hx, y = hy)
stdP = np.sqrt(muP.value) # Variance for Poisson Fitting

h.plot()
h.plot_exp(color = 'black')
B.plot_line(fitPoisson.xpl, fitPoisson.ypl, color = 'red')

plt.text(3.,40.0, r'$P(x) = \frac{N_P}{x!}\mu^x\cdot e^{-\mu}$',
fontsize=9)
plt.text(3.,30.0, r'$\mu =${:.2f}'.format(muP.value),fontsize=8)

plt.text(3.,20.0, r'$\sigma =${:.2f}'.format(stdP),fontsize=8)

plt.text(3.,11.5, r'$N_P =${:}'.format(int(normP.value)),fontsize=8)
B.pl.show()
