
# coding: utf-8

# In[11]:


###    Cholera Virus Model
###    
###  *Header Stuff Goes Here*
###
###
###


### Initializations ###

import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D
#import pymc3 as pm


### Initial Conditions ###
N = 10000                  # (t) may not be needed
S0 = 7000
A0 = 0
I0 = 0
R0 = 0
V0 = 0
BH0 = 1
BL0 = 1

V0 = np.array([S0, I0, A0, R0, V0, BH0, BL0])


### Parameter Definitions ###  

alpha = 0.02                  # [day^-1 person^-1], Rate of contaminated water consumption From G.G. Kolaye et al. Paper
kL = 10**5                     # [cells], LI V. cholera infectious dose
kH = kL / 50                   # [cells], HI V. cholera infectious dose 
chi = 0.8                      # [(24h)^-1] --> [day^-1], Rate of decay of HI to LI V. Cholera
delta = 0.03333                # [(30d)^-1] --> [day^-1], Death rate of V. Cholera in the environment 
pi = 4.491 * 10 ** (-5)        # [(61yr)^-1] --> [day^-1], Natural birth rate
mu = 4.491 * 10 ** (-5)        # [(61yr)^-1] --> [day^-1], Natural death rate
muC = 0.046                     # [day^-1], Mortality rate, symptomatic cholera
gamma = 0.2                    # [(5d)^-1], Recovery rate cholera                        
p = 0.79                       # [Unitless Preportion], Proportion of cases asymptomatic 
xiS = 1.3 * 10 ** 11           # [cells*d^-1], Rate of excretion of V. cholera, symptomatic patient 
xiA = 1.3 * 10 ** 8            # [cells*d^-1], Rate of excretion of V. cholera, asymptomatic patient 
omega = 0.00342                # [(0.8yr)^-1] --> [day^-1], Rate of waning of natural immunity
nu = 10 ** 1                   # [vaccines*d^-1], Rate of vaccination per day, guestimation for DC-Metro 9mill pop  
tau = 0.67                     # [Unitless Preportion], Vaccine efficacy
epsilon = 0.00137              # [(2yr)^-1] --> [day^-1], Rate of waning, vaccine induced immunity
theta = 0.08                   # [Unitless Preportion], Prop. symptmtic indvduals receiving antbtics
phi = 0.52                     # [Unitless Preportion], Relative rate of shedding, receiving antbtics
delta = 2.3                    # [Unitless Preportion], Relative rate of recovery, receiving antbtics
W = 15*N*365                   # [Liters], Size of water reservoir


### System of Diff Eqs ###

def dX_dt(X, t=0):
    
    S = X[0]
    I = X[1]
    A = X[2]
    R = X[3]
    V = X[4]
    BH = X[5]
    BL = X[6]
    
    dS_dt = ((pi*N)) +(omega*R) + (epsilon*V) - (alpha*S*((BL)/(kL + BL))) - (alpha*S*((BH)/(kH + BH))) - (mu*S) - (tau*nu)
    dI_dt = ((1-p)*alpha*S*(BL/(kL+BL))) + ((1-p)*alpha*S*((BH)/(kH+BH))) - (muC+mu+((1-theta)*gamma)+(theta*gamma*delta))*I
    dA_dt = (p*alpha*S*((BL)/(kL+BL))) + (p*alpha*S*((BH)/(kH+BH))) - (mu+gamma) * A
    dR_dt = (((1-theta)*gamma) + (theta*gamma*delta))*I + (gamma*A) - (mu+omega) * R
    dV_dt = (tau*nu) - (epsilon*V) - (mu*V)
    dBH_dt = ((phi*theta+(1-theta))*(xiS/W)*I) + ((xiA/W)*A) - (chi*BH)
    dBL_dt = (chi*BH) - (delta*BL)
    
    return np.array([dS_dt, dI_dt, dA_dt, dR_dt, dV_dt, dBH_dt, dBL_dt])


### System Integrations ###

TS = (10 ** 4)                    
t = np.linspace(0, 300, TS)          
ODEsSolved = odeint(dX_dt, V0, t)


### Figures and Plots ###

plt.figure(figsize=[17,8])
#plt.plot(t, ODEsSolved[:,0], 'b-', label='Susceptibles')
plt.plot(t, ODEsSolved[:,1]+ODEsSolved[:,2], 'g-', label='Total Infected')
#plt.plot(t, ODEsSolved[:,2], 'r-', label='Asymptomatic')
#plt.plot(t, ODEsSolved[:,3], 'c-', label='Recovered')
#plt.plot(t, ODEsSolved[:,4], 'm-', label='Vaccinated')
#plt.plot(t, ODEsSolved[:,5], 'y-', label='Hyper-infectious V. Cholerae')
#plt.plot(t, ODEsSolved[:,6], 'k-', label='Low-infectious V. Cholerae')
plt.xlim([0, 75])
#plt.ylim([0,10])
plt.tick_params(axis = 'both', direction = 'in', top = 1, right = 1)
plt.xlabel('Time (Days)', fontsize = '15')
plt.ylabel('Population' , fontsize = '15')
plt.title('Total Infected Populations vs Time' , fontsize = '20')
plt.legend(loc='best', prop={'size': 17})
plt.show()

plt.figure(figsize=[17,8])
#plt.plot(t, ODEsSolved[:,0], 'b-', label='Susceptibles')
plt.plot(t, ODEsSolved[:,1], 'g-', label='Symptomatic')
plt.plot(t, ODEsSolved[:,2], 'r-', label='Asymptomatic')
#plt.plot(t, ODEsSolved[:,3], 'c-', label='Recovered')
#plt.plot(t, ODEsSolved[:,4], 'm-', label='Vaccinated')
#plt.plot(t, ODEsSolved[:,5], 'y-', label='Hyper-infectious V. Cholerae')
#plt.plot(t, ODEsSolved[:,6], 'k-', label='Low-infectious V. Cholerae')
plt.xlim([0, 75])
#plt.ylim([0,10])
plt.tick_params(axis = 'both', direction = 'in', top = 1, right = 1)
plt.xlabel('Time (Days)', fontsize = '15')
plt.ylabel('Population' , fontsize = '15')
plt.title('Symptomatic & Asymptomatic Populations vs Time' , fontsize = '20')
plt.legend(loc='best', prop={'size': 17})
plt.show()

plt.figure(figsize=[17,8])
#plt.plot(t, ODEsSolved[:,0], 'b-', label='Susceptibles')
plt.plot(t, ODEsSolved[:,1], 'g-', label='Symptomatic')
plt.plot(t, ODEsSolved[:,2], 'r-', label='Asymptomatic')
plt.plot(t, ODEsSolved[:,3], 'y-', label='Recovered')
plt.plot(t, ODEsSolved[:,4], 'b-', label='Vaccinated')
#plt.plot(t, ODEsSolved[:,5], 'y-', label='Hyper-infectious V. Cholerae')
#plt.plot(t, ODEsSolved[:,6], 'k-', label='Low-infectious V. Cholerae')
#plt.xlim([0, .5])
#plt.ylim([0,10])
plt.tick_params(axis = 'both', direction = 'in', top = 1, right = 1)
plt.xlabel('Time (Days)', fontsize = '15')
plt.ylabel('Population' , fontsize = '15')
plt.title('Long Term Recovery & Vaccination Behavior' , fontsize = '20')
plt.legend(loc='best', prop={'size': 14})
plt.show()


plt.figure(figsize=[17,8])
plt.plot(t, ODEsSolved[:,0], 'b-', label='Susceptibles')
plt.plot(t, ODEsSolved[:,1], 'g-', label='Symptomatic')
plt.plot(t, ODEsSolved[:,2], 'r-', label='Asymptomatic')
plt.plot(t, ODEsSolved[:,3], 'c-', label='Recovered')
plt.plot(t, ODEsSolved[:,4], 'm-', label='Vaccinated')
plt.plot(t, ODEsSolved[:,5], 'y-', label='Hyper-infectious V. Cholerae')
plt.plot(t, ODEsSolved[:,6], 'k-', label='Low-infectious V. Cholerae')
plt.xlim([0, 120])
#plt.ylim([0,10])
plt.tick_params(axis = 'both', direction = 'in', top = 1, right = 1)
plt.xlabel('Time (Days)', fontsize = '15')
plt.ylabel('Population' , fontsize = '15')
plt.title('All Agents vs Time' , fontsize = '20')
plt.legend(loc='best', prop={'size': 15})
plt.show()


