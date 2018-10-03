# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
L=14
D=0.5 #Diffusion coefficient
dx=0.2
CFL = [1.4,2.1,1]
dt=0.02
c=D*(dt/(dx**2))
print(c)

#c=0.3
#dt=(c*(dx**2))/D
#print(dt)
iteration=10
def data_read(data):
    col1=[] #column 1
    col2 = [] #column 2
    with open(str(data+'.txt'),'r') as fl:
        for i in fl:
            col1.append(i.split('  ', 1)[0])
            col2.append(i.split('  ', 1)[1])
    return col1,col2
col1 = data_read("data1_auxiliary")[0]
for i in col1:
    col1[col1.index(i)]=float(i)
col2 = data_read("data1_auxiliary")[1]
for i in col2:
    col2[col2.index(i)]=float(i)
n=len(col1) #71
#print(n)
#Cited work from Travis Mitchell - put reference!
#def FI(D,x,t) :
#    #think of a way to include the range of t and create list of f for each t
#    fxt = []
#        
#    ff=np.exp(-np.array(x)**2/4.0/D/t)/np.sqrt(4.0*D*np.pi*t)
#    return(np.matrix(ff).transpose())  
#FI(D,col1,t)
def Explicit(n,xa,xb,bc,f0,c,xlist,iteration,dt):
    counter = 0
    f0_a=f0
    while counter < iteration:
        Z = np.zeros(n) #zero matrix size 71x71
        S = Z + (1 - 2*c)*np.eye(n) + c*np.eye(n,k=-1)+c*np.eye(n,k=1) #constructs S matrix without BC considered
#        print(S)
        #uses Dirichlet BC to finalise matrix S
#        S[0,:] = 0; S[0,0] = 1
        """Change left BC to Neumann's (1st order)"""
        S[0,:] = S[1,:]
        S[-1,:]= 0; S[-1,-1]=1
        f0 = np.transpose(np.asarray(f0_a))
        f1=np.matmul(S,f0)
        plt.plot(col1,f1,label=counter)
        counter +=1
        f1=f0_a
    plt.xlabel("x")
    plt.xlabel("f(x)")
#    plt.legend()
    plt.title("Explicit")
    plt.show()
    return f1

#def 1DPDE(n,xa,xb,bc,f0,c,iteration):
#    T = np.ones
#1DPDE(n,0,14,[0,0],col2,CFL[0],10)
def Implicit(n,xa,xb,bc,f1,c,xlist,iteration,dt):
    counter = 0
    f1_a=f1
    while counter < iteration:
        Z = np.zeros(n) #zero matrix size 71x71
        S = Z + (1 + 2*c)*np.eye(n) + (-c)*np.eye(n,k=-1)+(-c)*np.eye(n,k=1) #constructs S matrix without BC considered
        
        #uses Dirichlet BC to finalise matrix S
#        S[0,:] = 0; S[0,0] = 1
        """Change left BC to Neumann's (1st order)"""
        S[0,:] = S[1,:]
        S[-1,:]= 0; S[-1,-1]=1
#        print(S)
        f1 = np.transpose(np.asarray(f1_a))
        f0=np.matmul(S,f1)
        plt.plot(col1,f0,label=counter)
        counter +=1
        f0=f1_a
    plt.xlabel("x")
    plt.xlabel("f(x)")
    plt.title("Implicit")
#    plt.legend()
    plt.show()
    return f0    
Explicit(n,0,14,[0,0],col2,c,col1,iteration,dt)
Implicit(n,0,14,[0,0],col2,c,col1,iteration,dt)


"""Q1: How do I know that a solution is unstable?"""
"""Q2: Am I on the right tracK?"""
