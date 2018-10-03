# mech3750 2018  Part 2(ii) demo (can be used in assignments - AYK)      
# solving   df/dt=D*d2f/dx^2 by the spectral method    
# with the Dirichlet boundary conditions  f=0 at x=0 and x=xb 
#     using np.array()  // use np.matrix() for matrix multiplications    

import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


# the fundamental solution function (used for initial cond. & comparison)
def FI(x,t) :     
	return(np.exp(-np.array(x)**2/4/D/t)/np.sqrt(4*D*np.pi*t))

tfs=0.0025   # initial time for fundamental solution 
xfs=0.6      # location for the fundamental solution 

##########################################    BASIC PARAMETERS 
n=100    # n+1 is the number of grid points  
D=1      # diffusion coefficient
ns=5     # number of time steps    

xb=1      # domain limits    [0,xb]
t0=0      # initial time
t1=0.05   # the final time 

dt=(t1-t0)/ns    # time step, which does not have to be small
dx=xb/n 

cfl=dt*D/dx/dx   # cfl -- is not important for spectral methods

x=np.array(range(0,n+1))*dx   # x - grid:    n+1 grid points: [0,1,...,n]



################################# INITIAL AND BOUNDARY CONDITIONS
#  initial conditions 
f0=FI(x-xfs,tfs)         # initial condition: FI

#  boundary conditions   
f0[0]=0
f0[-1]=0

t=0     # initial time         
f=f0*1  # setting initial f to f0         



################################### SELECTING THE BASIS FOR EXPANSION 

# number of terms in the expansion is determined by the number of free values f[j]
n1=n-1    # in case of the Dirichlet boundary conditions since f[0] and f[-1] are fixed 


# select the expansion basis   ( sin(...) for Dirichlet b.c.)
ggg=np.zeros((n+1,n1))
for ii in range(0,n1): 
 ggg[:,ii]=(np.sin(x*np.pi*(ii+1)/xb))



################################### SPECTRAL EXPANSION  
# identifying and removing the linear component:
f_01=f0[0]+x*(f0[-1]-f0[0])/x[-1]    
f1=f0-f_01


#===================expand and check that f0==f00 
#           expanding  
A0=f0*0; 
for ii in range(0,n1):
 A0[ii]=np.sum(f1*ggg[:,ii])*2/n 

#           reconstructing f0
f00=f_01*1;
for ii in range(0,n1):
 f00=f00+A0[ii]*ggg[:,ii] 

#           always check that f00=f0
#df0=np.linalg.norm(f0-f00)
df0=np.ndarray.max(np.absolute(f0-f00))
print('>>>>>>>>>>>>>>>> init. cond. check:: df0 =',df0)


###############################  EVALUATION at any t>0 -------------------------
A1=A0*0;

for iii in range(0,ns+1):
 t=t0+dt*iii
 f=f_01                               # initial value is set to the linear component 
 for ii in range(0,n1):
   kk=np.pi*(ii+1)/xb
   lam=D*kk**2                        
   fi=np.exp(-t*lam)      	      # time exponents    
   A1[ii]=A0[ii]*fi                   # current expansion coefficients
   f=f+A1[ii]*ggg[:,ii]               # the sum 
 fs=FI(x-xfs,tfs+t)                   # current fundamental solution        
 plt.clf()
 ax=plt.gca()
 ax.plot(x,f0,'r-',x,f,'b',x,fs,'r--')
 plt.show()
 print('t =',t)

 # www=input("press ENTER")
 # plt.savefig('fff')

 
print('cfl = ',cfl)
print('====END====') 	