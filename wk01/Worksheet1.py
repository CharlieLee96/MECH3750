# Python as a calculator:
#   Here we can run in a script or just from the cmd line
#e.g:
"""
>> 2 + 2
>> 4.0 - 2
>> 5 / 2    # Be careful with Integer Devision if you're using Python2!
>> a = 5
>> b = 5.0
>> type(a)
>> type(b)

# Strings and lists
>> a = 'Test String'
>> type(a)

>> B = [0.3, 2, 5, 6.7]
>> B[0]     # Indices start from 0! Not like matlab, i.e. starting from 1
"""
# Unlike Matlab, we need to include external libraries
#   this allows us to select only the sources you need
#   i.e. more efficient. For neatness of code, place
#        imports all at top of script file.
import math
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

print("Task 1:")
def myFactorial(n):
    ans = 1
    if (n == 0):
        return ans
    else:
        for i in range(1,n+1):
            ans *= i
        return ans

x = 5
myFactor = myFactorial(x)
exactFac = math.factorial(x)
if abs(myFactor - exactFac) < 1e-6:
    print("Correct! %d factorial = %d") % (x, myFactor)
else:
    print("Try again..... :(")



print("Task 2:") 
n = 6
x = 0.1
total = 0
for j in range (n+1):
    total=total+x**j/math.factorial(j)

exact=math.exp(x)
print ('Exact=%.6f, Approx=%.6f') % (exact, total)

print("Task 3:") 
def exptaylor(n,x):
    total=0
    for j in range (n+1):
        total=total+x**j/math.factorial(j)
    return total

n=3
x=0.1
approx=exptaylor(n,x)
exact=math.exp(x)
print ('Exact=%.6f, Approx=%.6f') % (exact, total)

print("Task 4:") 
xvals = np.arange(-3, 3, 0.1)
myY = [exptaylor(5,xvals[i]) for i in range(len(xvals))]
yvals = np.exp(xvals)

plt.plot(xvals, yvals)
plt.plot(xvals, myY,'r.')
plt.show()

print("Task 5:")
m=4
dx=0.25
n=4
dy=0.25

x=np.zeros( (2*m+1,2*n+1) )
y=np.zeros( (2*m+1,2*n+1) )
f=np.zeros( (2*m+1,2*n+1) )

for i in range(2*m+1):
    for j in range(2*n+1):
        xx=dx*( i-m )
        yy=dy*( j-n )
        ff=xx**2+yy**2
        #print (i,j,xx,yy,ff)
        x[i,j]=xx
        y[i,j]=yy
        f[i,j]=ff

fig = plt.figure(2)
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(x, y, f, rstride=2, cstride=2)
plt.show()

















