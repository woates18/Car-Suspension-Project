#!/usr/bin/env python
# coding: utf-8

# In[816]:


import matplotlib.pyplot as plt
import numpy as np
import math
import sympy as sym
from scipy.integrate import odeint
from numpy.fft import fft, ifft
from sympy import DiracDelta, symbols, diff, sin
from scipy.integrate import solve_ivp
from scipy.optimize import minimize

# Parameters

mc = 1500 # mass of car
mcs =20 # mass of car seat
mp = 80 # mass of person

kcs = 21065000         # spring constant of car seat
ccs = 1000        # damper value of car seat
ks = 555000               # spring constant of carand ground
cs = 1000 # damper value of car and ground


'''Equation of motion'''
#EOM in State space form
def dx_dt(x, t):
    u = 0.05 * np.sin(2 * math.pi * t)
    ud = 0.05 * np.cos(2 * math.pi * t) 
    return [x[1], 
            ((-2*kcs-2*ks) / mc)*x[0] + (-2*cs - 2*ccs)/mc*x[1] + kcs/mc*x[2] + ccs/mc*x[3]
            + kcs/mc * x[4] + ccs/mc*x[5] + (2*ks/mc)*(u) + (2*cs/mc)*ud,
            x[3],
            kcs/(mp+mcs)*x[0] + (ccs / (mp + mcs))*x[1] - kcs / (mp + mcs)*x[2] - ccs / (mp + mcs)*x[3],
            x[5],
            kcs/mcs*x[0] + ccs/mcs*x[1] - kcs/mcs*x[4] - ccs/mcs*x[5]]
    

def dx2_dt(x, t):
    if t < 0.01:
        u = 1
    else:
        u = 0
    ud = 0
    ud = 0
    return [x[1], 
            ((-2*kcs-2*ks) / mc)*x[0] + (-2*cs - 2*ccs)/mc*x[1] + kcs/mc*x[2] + ccs/mc*x[3]
            + kcs/mc * x[4] + ccs/mc*x[5] + (2*ks/mc)*(u) + (2*cs/mc)*ud,
            x[3],
            kcs/(mp+mcs)*x[0] + (ccs / (mp + mcs))*x[1] - kcs / (mp + mcs)*x[2] - ccs / (mp + mcs)*x[3],
            x[5],
            kcs/mcs*x[0] + ccs/mcs*x[1] - kcs/mcs*x[4] - ccs/mcs*x[5]]

# Initial Conditions
x0 = [0, 0, 0, 0, 0, 0]


# Sampling Rate
sr = 100

# Sampling Interval (step size)
ts = 1/sr

t = np.arange(0,10,ts) # Time constant


# Integrate to solve for state vector
zt = odeint(dx_dt, x0, t)
zt2 = odeint(dx2_dt, x0, t)

#Position
x1 = zt[:,0]
x2 = zt[:,2]
x3 = zt[:,4]

#Position
x12 = zt2[:,0]
x22 = zt2[:,2]
x32 = zt2[:,4]
# Velocity
x1dot = zt[:,1]
x2dot = zt[:,3]
x3dot = zt[:,5]


'''Plotting for bumpy road'''
# Plot Results
u,ud = bumpy_road(t)
plt.figure(figsize=(8, 4))
plt.plot(t,x1)
plt.xlabel('Time (s)')
plt.ylabel('x1 displacment (m)')
plt.title("Car chasis displacment graph on a Modeled Bumpy road")
plt.xlim([0, 10]);
#plt.ylim([-0.06, 0.06]);

plt.figure(figsize=(8, 4))
plt.plot(t,x2)
plt.xlabel('Time (s)')
plt.ylabel('x2 displacment (m)')
plt.title("Drivers Seat displacment graph on a Modeled Bumpy road")
plt.xlim([0, 10]);
#plt.ylim([-0.06, 0.06]);

plt.figure(figsize=(8, 4))
plt.plot(t,x3)
plt.xlabel('time (s)')
plt.ylabel('x3 displacment (m)')
plt.title("Passenger Seat displacment graph on a Modeled Bumpy road")
plt.xlim([0, 10]);
#plt.ylim([-0.06, 0.06]);
'''plotting for impulse'''

# Plot Results
plt.figure(figsize=(8, 4))
plt.plot(t,x12)
plt.xlabel('Time (s)')
plt.ylabel('x1 displacment (m)')
plt.title("Car chasis displacment graph on a Modeled Speed Bump")
plt.xlim([0, 10]);
#plt.ylim([-0.25, 0.25]);

plt.figure(figsize=(8, 4))
plt.plot(t,x22)
plt.xlabel('Time (s)')
plt.ylabel('x2 displacment (m)')
plt.title("Drivers Seat displacment graph on a Modeled Speed Bump")
plt.xlim([0, 10]);
#plt.ylim([-0.25, 0.25]);

plt.figure(figsize=(8, 4))
plt.plot(t,x32)
plt.xlabel('Time (s)')
plt.ylabel('x3 displacment (m)')
plt.title("Passenger Seat displacment graph on a Modeled Speed Bump")
plt.xlim([0, 10]);
#plt.ylim([-0.25, 0.3]);

'''Input displacment plots'''
#impulse plot
u,ud = speed_bump(t)
plt.figure(figsize=(8, 4))
plt.plot(t,u)
plt.xlabel('Time (s)')
plt.ylabel('Displacment (m)')
plt.title("Speed Bump Displacment Modeled as an Impulse")
plt.xlim([-0.3, 0.3]);
# plt.ylim([-0.03, 0.03]);

#bumpy road plot
u = 0.05 * np.sin(2 * math.pi * t)
plt.figure(figsize=(8, 4))
plt.plot(t,u)
plt.xlabel('Time (s)')
plt.ylabel('Displacment (m)')
plt.title("Bumpy Road Modeled as a Sin Wave")
plt.xlim([0, 10]);
# plt.ylim([-0.03, 0.03]);



# In[ ]:







# In[ ]:





# In[ ]:




