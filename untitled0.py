# Lab 03 - Millikan Oil Drop
import nump as np
import LT.box as B
import matplotlib.pyplot as plt


fileName = 'Lab03.data'
f = B.get_file(fileName)

d_c = f.par.get_value('d_c')/100 # Plate Separation (cm)->(m)
dd_c = f.par.get_value('dd_c')/100.
g = f.par.get_value('g') # Gravitational Acceleration (m/s^2)
rho = f.par.get_value('rho') # Oil Density (kg/m^3)
b = f.par.get_value('b') # Constant (N/m)
dt = f.par.get_value('dt') # Error in Time (s)
T_C = f.par.get_value('T_C') # Chamber Temperature (Celsius)
dT_C = 0.2
P = f.par.get_value('P') # Barometric Pressure (Pa)
P0 = 101325. # Atmospheric Pressure (Pa)
dP = (P-P0)
V = f.par.get_value('V') # Capacitor Voltage [V]
dV = 1.

d_f = B.get_data(f, 'd_f')/1000. # Falling Distance (mm)->(m)
d_r = B.get_data(f, 'd_r')/1000. # Rising Distance (mm)->(m)
dd_f, dd_r = 0.1/1000., 0.1/1000.
t_f = B.get_data(f, 't_f') # Falling Time (s)
dt_f = dt
t_r = B.get_data(f, 't_r') # Rising Time (s)
dt_r = dt

elec = 1.60217662 *10.**(-19 ) # Elementary Charge Accepted Value

v_f = d_f/t_f # Falling Velocity
Pvf_Pdf = 1./t_f
Pvf_Ptf = -d_f/t_f**2
dv_f = np.sqrt((Pvf_Pdf*dd_f)**2+(Pvf_Ptf*dt_f)**2)

v_r = d_r/t_r # Rising Velocity
Pvr_Pdr = 1./t_r
Pvr_Ptr = -d_r/t_r**2
dv_r = np.sqrt((Pvr_Pdr*dd_r)**2+(Pvr_Ptr*dt_r)**2)

T_1 = 15.0 # (Celsius)
dT_1 = 0.2
Vis_1 = 1.8000*10.**(-5) # (N*s/m^2)
dVis_1 = 0.001*10.**(-5)
T_2 = 20.0 # (Celsius)
dT_2 = 0.2
1
Vis_2 = 1.8240*10.**(-5) # (N*s/m^2)
dVis_2 = 0.001*10.**(-5)

slope = (Vis_2 - Vis_1)/(T_2 - T_1)
Pslope_PT1 = (Vis_2-Vis_1)/(T_2-T_1)**2
Pslope_PVis1 = -1./(T_2-T_1)
Pslope_PT2 = -(Vis_2-Vis_1)/(T_2-T_1)**2
Pslope_PVis2 = 1./(T_2-T_1)
dslope = np.sqrt((Pslope_PT1*dT_1)**2+(Pslope_PT2*dT_2)**2+(Pslope_PVis1*dVis_1)**2+(Pslope_PVis2*dVis_2)**2)

eta = slope*(T_C-T_1)+Vis_1 # Viscosity of Air (N*s/m^2)
Peta_Pslope = (T_C-T_1)
Peta_PTC = slope
Peta_PT1 = -slope
Peta_PVis1 = 1.
deta = np.sqrt((Peta_Pslope*dslope)**2+(Peta_PTC*dT_C)**2+(Peta_PT1*dT_1)**2+(Peta_PVis1*dVis_1)**2)

# Charge on Droplet (C)
q = (4.*np.pi*rho*g*d_c/3.)*\
(np.sqrt((b/(2.*P))**2+9.*eta*v_f/(2.*g*rho))-b/(2.*P))**3*\
(v_f+v_r)/(v_f*V)

# Partial Derivatives with respect to:
# Capacitor Separation (d_c)
Pq_Pdc = (4.*np.pi*rho*g/3.)*\
(np.sqrt((b/(2.*P))**2+9.*eta*v_f/(2.*g*rho))-b/(2.*P))**3*\
(v_f+v_r)/(v_f*V)

# Pressure (P)
Pq_PP = (4.*d_c*g*np.pi*rho)*\
(b/(2.*P**2)-b**2/(4.*P**3*np.sqrt(b**2/(4.*P**2)+9.*eta*v_f/(2.*g*rho))))*\
(-b/(2.*P)+np.sqrt(b**2/(4.*P**2)+9.*eta*v_f/(2.*g*rho)))**2*\
(v_f+v_r)/(v_f*V)

# Eta
Pq_Peta = (9.*d_c*np.pi)*(-b/(2.*P)+np.sqrt(b**2/(4.*P**2)+9.*eta*v_f/(2.*g*rho)))**2*\
(v_f+v_r)/(V*np.sqrt(b**2/(4.*P**2)+9.*eta*v_f/(2.*g*rho)))

# Falling Velocity (v_f)
Pq_Pvf = (4.*d_c*g*np.pi*rho)*\
(-b/(2.*P)+np.sqrt(b**2/(4.*P**2)+9.*eta*v_f/(2.*g*rho)))**3/(3.*V*v_f)+\
9.*d_c*np.pi*eta*(-b/(2.*P)+np.sqrt(-b**2/(4.*P**2)+9.*eta*v_f/(2.*g*rho)))**2*\
(v_f+v_r)/(V*v_f*np.sqrt(b**2/(4.*P**2)+9.*eta*v_f/(2.*g*rho)))-\
4.*d_c*g*np.pi*rho*(-b/(2.*P)+np.sqrt(b**2/(4.*P**2)+9.*eta*v_f/(2.*g*rho)))**3*\
(v_f+v_r)/(3.*V*v_f**2)

# Rising Velocity (v_r)
Pq_Pvr = 4.*d_c*g*np.pi*rho*\
(-b/(2.*P)+np.sqrt(b**2/(4.*P**2)+9.*eta*v_f/(2.*g*rho)))**3/\
(3.*V*v_f)

# Voltage (V)
Pq_PV = (-4.*d_c*g*np.pi*rho)*\
(np.sqrt((b/(2.*P))**2+9.*eta*v_f/(2.*g*rho))-b/(2.*P))**3*\
2
(v_f+v_r)/(3.*V**2*v_f)
dq = np.sqrt((Pq_Pdc*dd_c)**2+(Pq_PP*dP)**2+(Pq_Peta*deta)**2+\
(Pq_Pvf*dv_f)**2+(Pq_Pvr*dv_r)**2+(Pq_PV*dV)**2)

    # Sorting Routine
qSort = np.zeros((np.size(q),2))
i_qSort = 0 # Sorting Index
while i_qSort < np.size(q): # Pair charges with error
    qSort[i_qSort,0]=q[i_qSort]
    qSort[i_qSort,1]=dq[i_qSort]
    i_qSort +=1
qSort = np.sort(qSort, axis=0) # Sort by charge q (sorting rows)
i_qSort = 0 # Reinitialize index
while i_qSort < np.size(q): # Charges Sorted with Error
    q[i_qSort] =qSort[i_qSort,0]
    dq[i_qSort]=qSort[i_qSort,1]
    i_qSort += 1

#print("The charge of the droplet is:")
#print(q)
print('\ n The minimum charge is: {:4.2e}'.format(np.min(q)))
# Plotting Sequence 'Charge on Droplet per Trial'
#fig = plt.figure()

#ax = fig.add_subplot(111)
x = np.arange(0,np.size(q),1) # Domain for Plot
#ax.plot(x,q, 'k.')
B.plot_exp(x, q, dq,elinewidth=0.5, ecolor='black',color='b',marker='.')
#plt.xlabel(r'Droplet Trial')
B.pl.xlabel(r'Droplet Trial')
#plt.ylabel(r'Charge of Droplet [C]', rotation = 'vertical')
B.pl.ylabel(r'Charge of Droplet [C]', rotation = 'vertical')
#plt.title('Charge on Droplet per Trial')
B.pl.title('Charge on Droplet per Trial')
#plt.axis([0.8*np.min(x), 1.2*np.max(x), 0.8*np.min(q), 1.2*np.max(q)])
B.pl.xlim(0.8*np.min(x), 1.2*np.max(x))
B.pl.ylim(0.8*np.min(q), 1.2*np.max(q))
#plt.grid(True)
B.pl.grid(True)
#plt.show()
B.pl.show()



# Determining Average Minimum Charge
Points = 17
q1 = np.zeros(Points) # Minimum Charge Array
dq1 = np.zeros(Points)
i = 00 # Start Index For q1
j = 00 # Start Index for q

while i < Points:
    q1[i]=q[j]
    dq1[i]=dq[j]
    i+=1
    j+=1
    
def WeightedAverage(a, da): # Defining Weighted Average Formula
        i = 0
        N = 0.
        D = 0.
        for i in range(0,np.size(a)):
            N+=(a[i]/da[i]**2.)
            D+=(1./da[i]**2.)
    
        x_bar = N/D
        dx_bar = np.sqrt(1./D)
        return x_bar, dx_bar

q1Avg, dq1Avg = WeightedAverage(q1,dq1)
print ('The weighted average of the minimum charge is: {:4.2e} +/- {:4.2e}'.format(q1Avg, dq1Avg))
i_n = 0 # Index for Integer Array n
n = np.zeros(np.size(q))

while i_n < np.size(q):
    n[i_n] = round(q[i_n]/q1Avg) # Round to Nearest Integer
    i_n += 1


# Plotting Sequence'Integer Multiple of Average Minimum Charge per Drop'
fig = plt.figure()
ax = fig.add_subplot(111)
x = np.arange(0,np.size(n),1) # Domain for Plot
ax.plot(x,n, 'k.')
plt.xlabel(r'Droplet Trial')
plt.ylabel(r'Integer Multiple', rotation = 'vertical')
plt.title('Integer Multiple of Average Minimum Charge per Drop')
plt.axis([0.8*np.min(x), 1.2*np.max(x), 0.8*np.min(n), 1.2*np.max(n)])
plt.grid(True)
plt.show()

# Determining the Elementary Charge
elecExp = q/n
PelecExp_Pq = 1./n
PelecExp_Pn = -q/n**2
delecExp = np.sqrt((PelecExp_Pq*dq)**2+(PelecExp_Pn*np.sqrt(n))**2)
delecExp=np.abs(delecExp)
# Plotting Sequence 'Elementary Charge (q/n) per Drop' without Error Bars
fig = plt.figure()
ax = fig.add_subplot(111)
4
x = np.arange(0,np.size(elecExp),1) # Domain for Plot
ax.plot(x,elecExp, 'k.')
plt.xlabel(r'Droplet Trial')
plt.ylabel(r'Charge [C]', rotation = 'vertical')
plt.title('Elementary Charge (q/n) without Error Bars')
plt.axis([0.8*np.min(x), 1.2*np.max(x), 0.3*np.min(elecExp)-np.min(delecExp), 1.2*np.max(elecExp)+np.max(delecExp)])
plt.grid(True)
plt.show()



''' Plotting Sequence 'Elementary Charge (q/n) per Drop' with Error Bars
B.plot_exp(x, elecExp, delecExp,elinewidth=0.5, ecolor='black',color='b',marker='.')
B.pl.xlabel(r'Droplet Trial')
B.pl.ylabel(r'Charge [C]', rotation = 'vertical')
B.pl.title('Elementary Charge (q/n) with Error Bars')
B.pl.xlim(0.8*np.min(x), 1.2*np.max(x))
B.pl.ylim(0.3*np.min(elecExp)-np.min(delecExp), 1.2*np.max(elecExp)+np.max(delecExp))
B.pl.grid(True)
B.pl.show()
#'''


elecExp, delecExp = WeightedAverage(elecExp, delecExp)
print ('Elementary Charge (q/n) based on weighted average: {:4.2e} +/- {:4.2e}'.format(elecExp, delecExp))
print ('Relative uncertainty is: {:4.2f} %'.format(delecExp/elecExp*100.))
Deviation = (elec-elecExp)*100./elec
print ('Deviation from the established value is: {:4.2f} %'.format(Deviation))
print ('A total of {:} trials were conducted.'.format(np.size(q)))