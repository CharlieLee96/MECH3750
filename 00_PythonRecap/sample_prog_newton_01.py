# Author: Pr. Tony Roberts
#	Git Maintain: @TravisMitchell
#	    24-07-2018: + Created on git
#			+ Altered from for loop to tolerance

#	Date: 24-07-2018

# TASK:
#	Uses Newton's method to find the roots of
#	a function:
#		x(n+1) = x(n) - f(x(n))/df(x(n))
#	Iterate till convergence

# Import necessary modules
import math
import numpy as np 
import matplotlib.pyplot as plt 

def fun1(x):
    """ 
	We will solve the roots of this function 
    """
    return x**3-x

def deriv1(x):
    """ 
	We can use analytic or numerical derivative here
	but this function returns the derivative of
	fun1(x)	
    """
    return 3*x**2-1


tol = 1e-4

xguess=np.zeros(2)  
x0 = 1.8			#Initial guess

xguess[0]=x0            	
xguess[1]=x0 - fun1(x0)/deriv1(x0)
err = xguess[1] - xguess[0]

i = 0				#Counter
while abs(err) > tol:             
    i += 1
    xn = xguess[i]
    xguess = np.append( xguess, xn - fun1(xn)/deriv1(xn) )
    err = xguess[i+1] - xguess[i]
    print ('Iteration=%i, x=%.4f, dx=%.4f') % (i,xguess[i+1], err)

print ('###### Convergence to Analytic:')
xexact=1                        #Here we know the answer
n = len(xguess)

error = np.zeros((n))  	#Initialize array
ivalue= np.zeros((n))  	#Initialize array of ivalues for plotting
for i in range(n):       
   error[i]=xguess[i]-xexact
   ivalue[i]=i                  
   print ('Iteration=%i: Error=%.4e') % (i,error[i])

# https://matplotlib.org/api/_as_gen/matplotlib.pyplot.html

plt.figure(1)
plt.plot(ivalue,error,'o')
plt.title('Newtons method', fontsize=20)
plt.text(1, 0.6, r'$x_{n+1}=x_n-f(x_n)/\frac{df}{dx}(x_n)$', fontsize=20)  
#The $ signs indicate latex code
plt.xlabel('Iteration', fontsize=15)
plt.ylabel('Error', fontsize=15)
plt.ylim([-0.1, 0.9])
plt.show()
# End of code

