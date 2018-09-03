
# coding: utf-8

# In[3]:


###   World War III or World War Z? -  The Complex Dynamics of Doom
###   (Urbano França, Gabriela Michel, C.Brandon Ogbunu, Sean Robinson)

### Initialzations ###

import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


# In[10]:


### Fig. 1 Initial Definitions & Conditions ###        

S0 = 1.6 * 10 ** 6                # Number of susceptibles
P0 = 0                            # Number of susceptibles that are panicked
F0 = 0                            # Number of susceptibles that have converted into fighters
Z0 = 1                            # Number of zombies

b_S = (365/1.6) * 10 ** (-8)      # Suseptible Encounters Per Year in Manhattan
b_R = 1/25                        # Birth Rate of S 
omega = 0.9*(3/4)*b_S             # Rate of encounter S and P * probability of panicking
lambd = 0.1*(3/2)*b_S            # Rate of encounter S and F * probability of becoming a fighter
beta_P = 0.25*b_S                 # Rate of encounter S and Z * probability of panicking
beta_F = 0.1*b_S                  # Rate of encounter S and Z * probability of becoming fighter
beta_Z = 0.65*b_S                 # Rate of encounter S and Z * probability of becoming a Zombie
alpha = 0.9*(1/5)*b_S             # Rate of encounter P and Z * probability of becoming a Zombie
V = 0.9*(2)*b_S                   # Rate of encounter F and Z * probability of killing a Zombie

mu_S = 1/55                      # Mortality rate of i (i => S, P, F, Z)
mu_P = 1/55                      #                 ''
mu_F = 2/55                      #                 ''
mu_Z = 1/3                       #                 ''


### System Equations ==> Eq.4, Eq.5, Eq.6, Eq.7 ###

def dX_dt(X, t=0):
    S = X[0]
    P = X[1]
    F = X[2]
    Z = X[3]
    dS_dt = b_R *(1 - mu_S) * S + b_R * (P + F) - ((beta_P + beta_F + beta_Z) * Z + omega * P + lambd * F) * S
    dP_dt = ((omega * P + beta_P * Z) * S) - ((alpha * Z + mu_P) * P)
    dF_dt = ((lambd * F + beta_F * Z) * S) - (mu_F * F)
    dZ_dt = ((alpha * P + beta_Z * S) * Z) - ((V * F + mu_Z) * Z)
    
    return np.array([dS_dt, dP_dt, dF_dt, dZ_dt])

V0 = np.array([S0, P0, F0, Z0])     # Initial Condition Vector


### System Integration ###

N = (1 * 10 ** 5)                   # Number of Time Steps
t = np.linspace(0, 100, N)           # Time Initialization
OD = odeint(dX_dt, V0, t) 


### Fig 1 Ploting ###

plt.figure(figsize=[15,7])
plt.plot(t, OD[:,0], 'b-', label='Susceptible')
plt.plot(t, OD[:,1], 'y-', label='Panicked')
plt.plot(t, OD[:,2], 'r-', label='Fighters')
plt.plot(t, OD[:,3], 'g-', label='Zombies')
plt.plot(t, (1.6 * 10 ** 6)*np.ones(N), 'k--', label = 'Initial Population')
#plt.xlim([0, .5])
#plt.ylim([0,2e6])
plt.xlabel('Time (Days)', fontsize = '15')
plt.ylabel('Population' , fontsize = '15')
plt.title('Zombie Apocalypse - Fig. 1' , fontsize = '20')
plt.legend(loc='best')
plt.show()




# In[5]:


### Fig. 2 Initial Definitions & Conditions ###

S0 = 1.6 * 10 ** 6                # Number of susceptibles
P0 = 0                            # Number of susceptibles that are panicked
F0 = 0                            # Number of susceptibles that have converted into fighters
Z0 = 1                            # Number of zombies

b_S = (365/2.133) * 10 ** (-8)    # Suseptible Encounters Per Year in Manhattan
b_R = 1/25                        # Birth Rate of S 
omega = 0.1*(3/4)*b_S             # Rate of encounter S and P * probability of panicking
lambd = 0.25*(3/2)*b_S            # Rate of encounter S and F * probability of becoming a fighter
beta_P = 0.25*b_S                 # Rate of encounter S and Z * probability of panicking
beta_F = 0.1*b_S                  # Rate of encounter S and Z * probability of becoming fighter
beta_Z = 0.65*b_S                 # Rate of encounter S and Z * probability of becoming a Zombie
alpha = 0.9*(1/5)*b_S             # Rate of encounter P and Z * probability of becoming a Zombie
V = 0.9*(2)*b_S                   # Rate of encounter F and Z * probability of killing a Zombie

mu_S = 1/55                      # Mortality rate of i (i => S, P, F, Z)
mu_P = 1/55                      #                 ''
mu_F = 2/55                      #                 ''
mu_Z = 1/3                       #                 ''


### System Equations ==> Eq.4, Eq.5, Eq.6, Eq.7 ###

def dX_dt(X, t=0):
    S = X[0]
    P = X[1]
    F = X[2]
    Z = X[3]
    dS_dt = b_R *(1 - mu_S) * S + b_R * (P + F) - ((beta_P + beta_F + beta_Z) * Z + omega * P + lambd * F) * S
    dP_dt = ((omega * P + beta_P * Z) * S) - ((alpha * Z + mu_P) * P)
    dF_dt = ((lambd * F + beta_F * Z) * S) - (mu_F * F)
    dZ_dt = ((alpha * P + beta_Z * S) * Z) - ((V * F + mu_Z) * Z)
    
    return np.array([dS_dt, dP_dt, dF_dt, dZ_dt])

V0 = np.array([S0, P0, F0, Z0])     # Initial Condition Vector


### System Integration ###

N = (1 * 10 ** 5)                   # Number of Time Steps
t = np.linspace(0, 25, N)           # Time Initialization
OD = odeint(dX_dt, V0, t) 


### Fig 2 Ploting ###

plt.figure(figsize=[12,7])
plt.plot(t, OD[:,0], 'b-', label='Susceptible')
plt.plot(t, OD[:,1], 'y-', label='Panicked')
plt.plot(t, OD[:,2], 'r-', label='Fighters')
plt.plot(t, OD[:,3], 'g-', label='Zombies')
plt.plot(t, (1.6 * 10 ** 6)*np.ones(N), 'k--', label = 'Initial Population')
#plt.xlim([0, .5])
#plt.ylim([0,2e6])
plt.xlabel('Time (Days)', fontsize = '15')
plt.ylabel('Population' , fontsize = '15')
plt.title('Zombie Apocalypse - Fig. 2' , fontsize = '20')
plt.legend(loc='best')
plt.show()



# In[33]:


### Fig. 3 Initial Definitions & Conditions ###

S0 = 1.6 * 10 ** 6                # Number of susceptibles
P0 = 0                            # Number of susceptibles that are panicked
F0 = 0                            # Number of susceptibles that have converted into fighters
Z0 = 1                            # Number of zombies

b_S = (365/2.13) * 10 ** (-6)     # Suseptible Encounters Per Year in Manhattan
b_R = 1/25                        # Birth Rate of S 
omega = 0.9*(3/4)*b_S             # Rate of encounter S and P * probability of panicking
lambd = 0.25*(3/2)*b_S            # Rate of encounter S and F * probability of becoming a fighter
beta_P = 0.25*b_S                 # Rate of encounter S and Z * probability of panicking
beta_F = 0.1*b_S                  # Rate of encounter S and Z * probability of becoming fighter
beta_Z = 0.65*b_S                 # Rate of encounter S and Z * probability of becoming a Zombie
alpha = 0.9*(1/5)*b_S             # Rate of encounter P and Z * probability of becoming a Zombie
V = 0.9*(2)*b_S                   # Rate of encounter F and Z * probability of killing a Zombie

mu_S = 1/55                       # Mortality rate of i (i => S, P, F, Z)
mu_P = 1/55                       #                 ''
mu_F = 2/55                       #                 ''
mu_Z = 1/3                        #                 ''


### System Equations ==> Eq.4, Eq.5, Eq.6, Eq.7 ###

def dX_dt(X, t=0):
    S = X[0]
    P = X[1]
    F = X[2]
    Z = X[3]
    dS_dt = b_R *(1 - mu_S) * S + b_R * (P + F) - ((beta_P + beta_F + beta_Z) * Z + omega * P + lambd * F) * S
    dP_dt = ((omega * P + beta_P * Z) * S) - ((alpha * Z + mu_P) * P)
    dF_dt = ((lambd * F + beta_F * Z) * S) - (mu_F * F)
    dZ_dt = ((alpha * P + beta_Z * S) * Z) - ((V * F + mu_Z) * Z)
    
    return np.array([dS_dt, dP_dt, dF_dt, dZ_dt])

V0 = np.array([S0, P0, F0, Z0])     # Initial Condition Vector


### System Integration ###

N = (1 * 10 ** 5)                   # Number of Time Steps
t = np.linspace(0, .8, N)           # Time Initialization
OD = odeint(dX_dt, V0, t) 


### Fig. 3 Ploting ###

plt.figure(figsize=[15,7])
plt.plot(t, OD[:,0], 'b-', label='Susceptible')
plt.plot(t, OD[:,1], 'y-', label='Panicked')
plt.plot(t, OD[:,2], 'r-', label='Fighters')
plt.plot(t, OD[:,3], 'g-', label='Zombies')
plt.plot(t, (1.6 * 10 ** 6)*np.ones(N), 'k--', label = 'Initial Population')
plt.xlim([0, .3])
plt.ylim([0,2e6])
plt.xlabel('Time (Days)', fontsize = '15')
plt.ylabel('Population' , fontsize = '15')
plt.title('Zombie Apocalypse - Fig. 3' , fontsize = '20')
plt.legend(loc='best')
plt.show()



# In[37]:


### Fig. 4 Initial Definitions & Conditions ###

S0 = 1.6 * 10 ** 6                # Number of susceptibles
P0 = 0                            # Number of susceptibles that are panicked
F0 = 0                            # Number of susceptibles that have converted into fighters
Z0 = 1                            # Number of zombies

b_S = (365/1.6) * 10 ** (8)     # Suseptible Encounters Per Year in Manhattan
b_R = 1/25                        # Birth Rate of S 
omega = 0.9*(3/4)*b_S             # Rate of encounter S and P * probability of panicking
lambd = 0.1*(3/2)*b_S            # Rate of encounter S and F * probability of becoming a fighter
beta_P = 0.25*b_S                 # Rate of encounter S and Z * probability of panicking
beta_F = 0.1*b_S                  # Rate of encounter S and Z * probability of becoming fighter
beta_Z = 0.65*b_S                 # Rate of encounter S and Z * probability of becoming a Zombie
alpha = 0.9*(1/5)*b_S             # Rate of encounter P and Z * probability of becoming a Zombie
V = 0.9*(2)*b_S                   # Rate of encounter F and Z * probability of killing a Zombie

mu_S = 1/55                       # Mortality rate of i (i => S, P, F, Z)
mu_P = 1/55                       #                 ''
mu_F = 2/55                       #                 ''
mu_Z = 1/3                        #                 ''


### System Equations ==> Eq.4, Eq.5, Eq.6, Eq.7 ###

def dX_dt(X, t=0):
    S = X[0]
    P = X[1]
    F = X[2]
    Z = X[3]
    dS_dt = b_R *(1 - mu_S) * S + b_R * (P + F) - ((beta_P + beta_F + beta_Z) * Z + omega * P + lambd * F) * S
    dP_dt = ((omega * P + beta_P * Z) * S) - ((alpha * Z + mu_P) * P)
    dF_dt = ((lambd * F + beta_F * Z) * S) - (mu_F * F)
    dZ_dt = ((alpha * P + beta_Z * S) * Z) - ((V * F + mu_Z) * Z)
    
    return np.array([dS_dt, dP_dt, dF_dt, dZ_dt])

V0 = np.array([S0, P0, F0, Z0])     # Initial Condition Vector


### System Integration ###

N = (1 * 10 ** 5)                   # Number of Time Steps
t = np.linspace(0, 25, N)           # Time Initialization
OD = odeint(dX_dt, V0, t) 


### Fig. 4 Ploting ###

#plt.figure(figsize=[15,7])
#plt.plot(t, OD[:,0], 'b-', label='Susceptible')
#plt.plot(t, OD[:,1], 'y-', label='Panicked')
#plt.plot(t, OD[:,2], 'r-', label='Fighters')
#plt.plot(t, OD[:,3], 'g-', label='Zombies')
#plt.plot(t, (1.6 * 10 ** 6)*np.ones(N), 'k--', label = 'Initial Population')
#plt.xlim([0, .5])
#plt.ylim([0,2e5])
#plt.xlabel('Time (Days)', fontsize = '15')
#plt.ylabel('Population' , fontsize = '15')
#plt.title('Zombie Apocalypse - Fig. 4' , fontsize = '20')
#plt.legend(loc='best')
#plt.show()



# In[4]:


### Fig. 5 Initial Definitions & Conditions ###

S0 = 1.6 * 10 ** 6                # Number of susceptibles
P0 = 0                            # Number of susceptibles that are panicked
F0 = 0                            # Number of susceptibles that have converted into fighters
Z0 = 1                            # Number of zombies

b_S = (365/1.6) * 10 ** (8)     # Suseptible Encounters Per Year in Manhattan
b_R = 1/25                        # Birth Rate of S 
omega = 0.5*(3/4)*b_S             # Rate of encounter S and P * probability of panicking
lambd = 0.1*(3/2)*b_S            # Rate of encounter S and F * probability of becoming a fighter
beta_P = 0.25*b_S                 # Rate of encounter S and Z * probability of panicking
beta_F = 0.1*b_S                  # Rate of encounter S and Z * probability of becoming fighter
beta_Z = 0.65*b_S                 # Rate of encounter S and Z * probability of becoming a Zombie
alpha = 0.9*(1/5)*b_S             # Rate of encounter P and Z * probability of becoming a Zombie
V = 0.9*(2)*b_S                   # Rate of encounter F and Z * probability of killing a Zombie

mu_S = 1/55                       # Mortality rate of i (i => S, P, F, Z)
mu_P = 1/55                       #                 ''
mu_F = 2/55                       #                 ''
mu_Z = 1/3                        #                 ''


### System Equations ==> Eq.4, Eq.5, Eq.6, Eq.7 ###

def dX_dt(X, t=0):
    S = X[0]
    P = X[1]
    F = X[2]
    Z = X[3]
    dS_dt = b_R *(1 - mu_S) * S + b_R * (P + F) - ((beta_P + beta_F + beta_Z) * Z + omega * P + lambd * F) * S
    dP_dt = ((omega * P + beta_P * Z) * S) - ((alpha * Z + mu_P) * P)
    dF_dt = ((lambd * F + beta_F * Z) * S) - (mu_F * F)
    dZ_dt = ((alpha * P + beta_Z * S) * Z) - ((V * F + mu_Z) * Z)
    
    return np.array([dS_dt, dP_dt, dF_dt, dZ_dt])

V0 = np.array([S0, P0, F0, Z0])     # Initial Condition Vector


### System Integration ###

N = (1 * 10 ** 5)                   # Number of Time Steps
t = np.linspace(0, 25, N)           # Time Initialization
OD = odeint(dX_dt, V0, t) 


### Fig. 5 Ploting ###

#plt.figure(figsize=[15,7])
#plt.plot(t, OD[:,0], 'b-', label='Susceptible')
#plt.plot(t, OD[:,1], 'y-', label='Panicked')
#plt.plot(t, OD[:,2], 'r-', label='Fighters')
#plt.plot(t, OD[:,3], 'g-', label='Zombies')
#plt.plot(t, (1.6 * 10 ** 6)*np.ones(N), 'k--', label = 'Initial Population')
#plt.xlim([0, .5])
#plt.ylim([0,2e5])
#plt.xlabel('Time (Days)', fontsize = '15')
#plt.ylabel('Population' , fontsize = '15')
#plt.title('Zombie Apocalypse - Fig. 5' , fontsize = '20')
#plt.legend(loc='best')
#plt.show()


