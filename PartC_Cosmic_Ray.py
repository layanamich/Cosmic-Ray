import numpy as np
import LT.box as B
import matplotlib.pyplot as plt



fileName = 'PartC_1.txt'
f = B.get_file(fileName)



dV = f.par.get_value('dV') # Error in Voltage (V)
dA = f.par.get_value('dA') # Error in Zenith Angle (Degrees)
dt = f.par.get_value('dt') # Error in Time(s)
l = f.par.get_value('l') # Length of PMT Panel (cm)
dl = f.par.get_value('dl')
w = f.par.get_value('w') # Width of PMT Panel (cm)
dw = f.par.get_value('dw')
d = f.par.get_value('d') # Distance Between PMT Panels (cm)
dd = f.par.get_value('dd')


A = B.get_data(f, 'A') # Zenith Angle (Degrees)
n = B.get_data(f, 'n') # Counts
dn = np.sqrt(n)
t = B.get_data(f, 't') # Time (min.)


# Counts per min. (Experimental)
r = n/t
Pr_Pn = 1./t
Pr_Pt = -n/t**2
dr = np.sqrt((Pr_Pn*dn)**2+(Pr_Pt*dt)**2)


# Counts per min. (Theoretical)
def convDeg(Degrees):
    radian = Degrees*np.pi/180.0
    return radian
solidAngle = l*w/d**2


# Partial Derivatives of Solid Angle with Respect to:
PsA_Pl = w/d**2 # Length
PsA_Pw = l/d**2 # Width
PsA_Pd = -2.*l*w/d**3 # Distance


# Error of Solid Angle
dsA = np.sqrt((PsA_Pl*dl)**2+(PsA_Pw*dw)**2+(PsA_Pd*dd)**2)
x=np.linspace(0.0, 190.0, 89)
rTheo = 0.66*(np.cos(convDeg(x)))**2*solidAngle*l*w # Theoretical Rate


# Partial Derivatives of Theoretical Rate with Respect to:
PrTheo_PsA = 0.66*(np.cos(convDeg(x)))**2*l*w # Solid Angle
PrTheo_Pl = 0.66*(np.cos(convDeg(x)))**2*solidAngle*w # Length
PrTheo_Pw = 0.66*(np.cos(convDeg(x)))**2*solidAngle*l # Width

# Error of Theoretical Rate
drTheo = np.sqrt((PrTheo_PsA*dsA)**2+(PrTheo_Pl*dl)**2+(PrTheo_Pw*dw)**2)

# Plotting Sequence in Python 2.7
B.plot_exp(A, r, dr)

#B.plot_exp(x,rTheo,drTheo)
B.pl.title('Cosmic Ray Count Rate vs. Zenith Angle')
B.pl.xlabel('Zenith Angle (Degrees)')
B.pl.ylabel('Count Rate (Counts/min)', rotation = 'vertical')
B.pl.xlim(-10.0, 190.0)
B.pl.ylim(-.05, 1.1*(np.max(r)+np.max(dr)))
B.pl.grid(True)
plt.plot(A,r, label='Experimental Rate')
plt.plot(x,rTheo, label='Theoretical Rate')
plt.legend()
B.pl.show()

