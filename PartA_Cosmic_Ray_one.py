import numpy as np
import LT.box as B
import matplotlib.pyplot as plt
import math
# Parameter Module
import LT_Fit.parameters as P
# Genfit module.
import LT_Fit.gen_fit as G
import scipy.special as sp

def set_range(first_bin = math.pi/4, bin_width = math.pi/10, Nbins = 10):
    
    """
Function used to set the range and bin width of a histrogram
Input: first_bin = bin center of the first bin,
bin_width = bin width,
Nbins = total number of bins
Returns: A tuple used in the range key word defining
a histogram.
NOTE: For the histogram use the same number of bins:
example: h = B.histo( r, range = set_range(-5.0, 1, 11), bins = 11)
This created a histogram where the first bin is centered at -5.0,
the next at -4.0 etc. a total of 11 bins are created and the bin
center of the last one is at first_bin + (Nbins - 1)*bin_width = 5.0
"""

    rmin = first_bin - bin_width/2.0
    rmax = rmin + Nbins*bin_width
    return (rmin, rmax)

Nbins = 40
first_bin = -math.pi/2
bin_width = math.pi/Nbins
rmin = first_bin + bin_width/2.0
rmax = rmin + Nbins*bin_width

q = []
for m in range(Nbins):
    q.append(rmin)
    rmin += bin_width
N = []
A= 1
C= 100

for i in range(len(q)):
    Counts = A + C*(math.cos(q[i]))**2
    Counts =int(Counts)
    print(Counts)
    for j in range(Counts):
        N.append(q[i])

print(q)
h = B.histo(N, range = set_range(first_bin = -math.pi/2, bin_width = math.pi/Nbins,
Nbins = Nbins), bins = Nbins, title = 'Cosmic Ray Simulation', xlabel = 'Azimuth Angle',
ylabel = 'Perfect Number of Counts' )

# Putting the Bin Centers and Bin Contents into two Arrays
hx = h.bin_center
hy = h.bin_content
dy = np.sqrt(hy) # Error' of Bin Contents

h.plot()
h.plot_exp(color = 'black')
# Helper function

# Poisson Function Fitting
muP = P.Parameter(100.0, 'muP') # Initial Guess for mu
normP = P.Parameter(100.0, 'normP') # Initial Guess for normalization factor
def Poisson(x):
    value = normP()*muP()**(x)*np.exp(-muP())/sp.factorial(x)
    return value
fitPoisson = G.genfit(Poisson, [muP, normP], x = hx, y = hy)
stdP = np.sqrt(muP.value)# Variance for Poisson Fitting

h.plot()
h.plot_exp(color = 'black')
B.plot_line(fitPoisson.xpl, fitPoisson.ypl, color = 'red')

plt.text(1.0,100.6, r'$P(x) = \frac{N_P}{x!}\mu^x\cdot e^{-\mu}$',
fontsize=9)
plt.text(1.0,80.2, r'$\mu =${:.2f}'.format(muP.value),fontsize=8)
plt.text(1.0,60.3, r'$\sigma =${:.2f}'.format(stdP),fontsize=8)
plt.text(1.0,40.4, r'$N_P =${:}'.format(int(normP.value)),fontsize=8)
B.pl.show()
