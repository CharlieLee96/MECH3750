# MECH3750: Semester 2, 2018
# author: Travis Mitchell

# Part II: Parabolic PDE example 
# Solving the diffusion equation with finite difference
#       df/dt=D*d2f/dx^2 by finite difference    
# with the Dirichlet boundary conditions:  f=0 at x=0 and f=0 at x=xb 

import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# the fundamental solution function
def FI(D,x,t) :     
	ff=np.exp(-np.array(x)**2/4.0/D/t)/np.sqrt(4.0*D*np.pi*t)
	return(np.matrix(ff).transpose())  

#### STEP 1: DISCRETISATION ####
# define the system and parameters required to discretise #
n=100    # n+1 is the number of grid points  
D=1.0    # diffusion coefficient
cfl=0.4  # CFL number 


xa=0.0      # domain limits    [xa,xb]
xb=1.0
t0=0.00   # initial time 
tf=0.1   # the final time 

x0=0.6;

# Here we calculate dx and use our CFL criterion to #
#  define a time step that will suit our stability  #
#  requirements.                                    #
dx=xb/n 
dt=dx*dx*cfl/D 
x=np.array(range(0,n+1))*dx + xa  

#### STEP 2: INITIALISATION ####
# You may have a discrete set of points or a fn describing this #
#  Here we will use the a shifted fundamental solution for this #
#  You could code in here though a dirac delta function for e.g #
f0=FI(D,x-x0,t0+0.0025)  
print(f0)

#### STEP 3: RUN ####
### For the diffusion equation we can formulate as: ###
### Implicit, explicit, Crank-Nicolson              ###
### Here I will show you explicit, it is your job   ###
###  to code in the rest                            ###

## The explicit scheme is given by:            ##
##   f(i,jp1) = cfl* (f(ip1,j) - 2f(i,j) + f(im1,j)) + f(i,j)

## Finite difference using matrices ##
S = (1 - 2*cfl) * np.eye(x.size) + np.diag(cfl*np.ones(x.size-1),1) + np.diag(cfl*np.ones(x.size-1),-1)
# Re-construct first and last row for B.C.
S[0,:] = 0; S[0,0] = 1
S[-1,:]= 0; S[-1,-1]=1
print(S)
# Set I.C. and B.C.
t = t0
f = f0
f[0] = 0; f[-1] = 0 # Dirichlet B.C. equal to 0
iter = 0
while t<tf:
	f=S*f                                # finite difference -- step by step
	t=t+dt 
	iter=iter+1 
	if iter%100 == 0:	    
		fig = plt.figure()
		print([iter,t,tf])  		 
		plt.clf()
		plt.plot(x,f0,'r-',x,f,'b-')
		plt.draw()
		plt.waitforbuttonpress(0)
		plt.close(fig)
	
