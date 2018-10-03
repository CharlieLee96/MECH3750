# -*- coding: utf-8 -*-
"""
Spyder Editor

Work of Chul Lee (43951024)

"""
import numpy as np
from math import *
import matplotlib
import matplotlib.pyplot as plt
from numpy.linalg import inv
import sympy
L=14
D=0.5 #Diffusion coefficient
dx=0.2
#CFL = [1.4,2.1,1]
#dt=0.023 
#c=D*(dt/(dx**2))
#print(c)

c=0.001
print("c = ",c)
dt=(c*(dx**2))/D
#print(dt)
iteration=2/dt + 1
print("number of iterations = ",int(iteration))
def data_read(data):
    x=[] #column 1
    fx = [] #column 2
    with open(str(data+'.txt'),'r') as fl:
        for i in fl:
            x.append(i.split('  ', 1)[0])
            fx.append(i.split('  ', 1)[1])
    return x,fx
x = data_read("data2_investigation")[0]
print(len(x))
for i in x:
    x[x.index(i)]=float(i)
fx2 = data_read("data2_investigation")[1]
#print(f)
for i in fx2:
    fx2[fx2.index(i)]=float(i)
n=len(x) #71
x=np.asarray(x)
fx2=np.asarray(fx2)
 
    
def g(n,x,t):
        g=np.cos(((pi+2*n*pi)/28)*x)
#        print("g at n = ",n," = ",g)
#        print(np.shape(g))
        return g
    
def f(n,x,t,B):
        f=((exp(-(((pi+(2*n*pi))/28)**2)*t*D))*(B[n])*np.cos(((pi+2*n*pi)/28)*x))
        return f
def SpectralSol(x,f0,t,D,N):

    f_01=f0[0]+x*(f0[-1]-f0[0])/x[-1]    
    f1=f0-f_01
    B=np.zeros(71)
    gm=np.zeros([N-1,71]) #n x 71 matrix
#    Bnum=np.zeros(71)
#    Bdenom=np.zeros(71)
    fm=np.zeros([N-1,71]) #n x 71 matrix
    fxj=np.zeros(71)

    for n in range(0,N-1):
        #replace an empty matrix with g values (each row is n values)
        gm[n,:]=g(n,x,t) #g is a list of g values corresponding to all x values
#        print("gm at n=",gm)
#        print(np.shape(gm))
        #B[n] is An
        B[n]=(2/71)*((np.sum(f0[1:-2]*(gm[n,1:-2])))+(((f0[0]*gm[0,0])+(f0[70]*gm[n,70]))/2))
        fm[n,:]=f(n,x,t,B)
#    print("fm",fm)
#    print(np.shape(fm))
    
#    print("gm",gm)
    for n in range (0,len(x)-1):
#        print(fm[:,n])
        fxj[n]=np.sum(fm[:,n])
    try:
        fxj[-2] <= 0.000000000000001
    except ValueError:
        plt.figure(4)
        plt.plot(x, f0,"b", label='provided data' )
        plt.plot(x, fxj, "r", label='Solution data')
        plt.title(str("t= "+str(t)+", N= "+str(N)))
        plt.show()
        print("we have reached your destination at t = ",t)
        sys.exit()
            
#    for i in range(0,len(x)-1):
#        fxj[n]=np.sum(fm[:,i])
    if N == 9: #decide the harmonic you wish to plot
        plt.figure(4)
        plt.plot(x, f0,"b", label='provided data' )
        plt.plot(x, fxj, "r", label='Solution data')
        plt.title(str("t= "+str(t)+", N= "+str(N)))
        plt.show()
    return fxj
#SpectralSol(x,fx2,-8,0.5,7)
#SpectralSol(col1,col2,2,0.5,71)
t=5
while t > -10:
    for i in range(1,12):
        SpectralSol(x,fx2,t,0.5,i)
    t +=-0.5

def Spectral(t,k,x,fx,N,D):
    fmat = []  #matrix consisting of exp and cosine terms (row: exponential term with increasing kn values  column: increasing x values)
    B=np.zeros(71)
    gm=np.zeros([71,71])
    Bnum=np.zeros(71)
    Bdenom=np.zeros(71)
#    fx=np.transpose(fx)
#    print("fx",np.shape(fx))
    for i in x: #for each x value
        fmat_row=[] #creates an empty matrix to append values
        for n in range(N+1): #level of harmonics
#            print(type(n))
            e= np.exp(-(k**2)*D*t)#exponential term
#            print("e",np.shape(e))
            fterm=e*np.cos(((pi+2*n*pi)/28)*i) #fterm is a term multiplying cosine term with exponential term (not considering constant term)
#            print(np.shape(fterm))
#            print(type(fterm))
            fmat_row.append(fterm)
#        print(np.shape(fmat_row))
        fmat.append(fmat_row)
    print(fmat)
#    fmat_inv=inv(fmat)
    Amat=np.linalg.lstsq(fmat,fx)
    print (Amat)
    return Amat[0]
#Spectral(5,(pi+(2*n*pi)),x,fx2,5,D)

#for i in col3:
    