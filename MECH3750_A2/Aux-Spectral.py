from math import *
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('data1_auxiliary.txt')
x = data[:,0]
f0 = data[:,1]
f_01=f0[0]+x*(f0[-1]-f0[0])/x[-1]    
f1=f0-f_01
time=2
D=0.5
pi=np.pi

def g(n,x):
    answer=np.cos(((pi+2*n*pi)/28)*x)
    return answer

def f(n,x):
    answer=((exp(-(((pi+2*n*pi)/28)**2)*time*D))*(B[n])*np.cos(((pi+2*n*pi)/28)*x))
    return answer
                  
    
B=np.zeros(71)
gm=np.zeros([71,71])
Bnum=np.zeros(71)
Bdenom=np.zeros(71)
fm=np.zeros([70,71])
fxj=np.zeros(71)

for n in range(0,70):
    gm[n,:]=g(n,x)
    B[n]=(2/71)*((np.sum(f0[1:-2]*(gm[n,1:-2])))+(((f0[0]*gm[0,0])+(f0[70]*gm[70,70]))/2))
    fm[n,:]=f(n,x)

for n in range (0,70):
    fxj[n]=np.sum(fm[:,n])
    



        
plt.figure(4)
plt.plot(x, f0,"b", label='Initial Condition' )
plt.plot(x, fxj, "r", label='After two days')
plt.show()
