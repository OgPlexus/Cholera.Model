
# coding: utf-8

# In[1]:


###    LV Model Equaitons
###
###    dU/dt = aU - bUV
###    dV/dt = -cV + dbU*V
###
###    Assumptions (Wikipedia.Lotka-Volterra)
###      1) The prey population finds ample food at all times.
###      2) The food supply of the predator population depends entirely 
###         on the size of the prey population.
###      3) The rate of change of population is proportional to its size.
###      4) During the process, the environment does not change in favour 
###         of one species, and genetic adaptation is inconsequential.
###         Predators have limitless appetite.


# ## Initial Definitions

# In[2]:


from numpy import *
import pylab as p

### Definition of the parameters ###                       ### For Multiple Initial Values
                                                           ### Define an Array and Loop over It
# U       # Rabits
# V       # Foxes
a = 1.0   # natural growth rate of rabits, without foxes
b = 0.5   # natural dying rate of rabits due to fox predation
c = 1.5   # natural dying rate of foxes, without rabits
d = 0.7   # how many caught rabits let create a new fox

### Definition of the equations ###

def dX_dt(X,t=0):
    """ Return the growth rate of fox and rabbit populations"""
    return array([ a*X[0] -   b*X[0]*X[1] ,
                  -c*X[1] + d*b*X[0]*X[1] ])

### As Array U = X[0], V = X[1]


# ## Aux Calculations

# In[3]:


### Position Equalibrium ###

X_f0 = array([      0.,  0.])
X_f1 = array([ c/(d*b), a/b])

all(dX_dt(X_f0) == zeros(2)) and all (dX_dt(X_f1) == zeros(2))

### Jacobian Matrix ###

def d2X_dt2(X, t=0):
    """ Return the Jacobian matrix evaluated in X. """
    return array([[a -b*X[1],   -b*X[0]     ],
                  [b*d*X[1] ,   -c +b*d*X[0]] ])


A_f0 = d2X_dt2(X_f0)  
A_f1 = d2X_dt2(X_f1)
A_f0, A_f1    

## Eigenvalues & Period ###

lambda1, lambda2 = linalg.eigvals(A_f1)
lambda1, lambda2

T_f1 = 2*pi/abs(lambda1) 
T_f1


# ## System Integration

# In[4]:


from scipy import integrate
t = linspace(0, 15,  1000)              # time
X0 = array([10, 5])                     # initials conditions: 10 rabbits and 5 foxes
X, infodict = integrate.odeint(dX_dt, X0, t, full_output=True)
infodict['message']                     


# ## Plots 

# In[7]:


values  = linspace(0.3, 0.9, 5)                          # position of X0 between X_f0 and X_f1
vcolors = p.cm.autumn_r(linspace(0.3, 1., len(values)))  # colors for each trajectory

f2 = p.figure(figsize=(18,6))

p.subplot(122)
# plot trajectories
for v, col in zip(values, vcolors):
    X0 = v * X_f1                               # starting point
    X = integrate.odeint( dX_dt, X0, t)         # we don't need infodict here
    p.plot( X[:,0], X[:,1], lw=3.5*v, color='k', label='X0=(%.f, %.f)' % ( X0[0], X0[1]) )

# define a grid and compute direction at each point
ymax = p.ylim(ymin=0)[1]                        # get axis limits
xmax = p.xlim(xmin=0)[1]
nb_points   = 20

x = linspace(0, xmax, nb_points)
y = linspace(0, ymax, nb_points)

X1 , Y1  = meshgrid(x, y)                       # create a grid
DX1, DY1 = dX_dt([X1, Y1])                      # compute growth rate on the gridt
M = (hypot(DX1, DY1))                           # Norm of the growth rate 
M[M == 0] = 1.                                  # Avoid zero division errors 
DX1 /= M                                        # Normalize each arrows
DY1 /= M

p.title('Contours in Fox-Rabit Space', fontsize = 20)
#p.quiver(X1, Y1, DX1, DY1, M, pivot='mid', cmap=p.cm.jet) #quiver plots work sometimes
p.xlabel('Number of Rabbits', fontsize = 14)
p.ylabel('Number of Foxes', fontsize = 14)
#p.legend()
#p.grid()
p.xlim(0, xmax)
p.ylim(0, ymax)
#f2.savefig('rabbits_and_foxes_2.png')

rabbits, foxes = X.T
p.subplot(121)
p.plot(t, rabbits, 'r-', label='Rabbits')
p.plot(t, foxes  , 'b-', label='Foxes')
#p.grid()
p.legend(loc='upper right')
p.legend(prop={'size': 15})
p.xlabel('Time', fontsize = 14)
p.ylabel('Population', fontsize = 14)
p.title('Evolution of Fox and Rabbit Populations', fontsize = 20)
#f1.savefig('rabbits_and_foxes_1.jpg')
p.show()

