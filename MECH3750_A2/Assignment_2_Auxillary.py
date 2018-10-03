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
L=14
D=0.5 #Diffusion coefficient
dx=0.2
#CFL = [1.4,2.1,1]
#dt=0.02
#c=D*(dt/(dx**2))
#print(c)

c=0.5
print("c = ",c)
dt=(c*(dx**2))/D
#print(dt)
iteration=2/dt + 1
print("number of iterations = ",int(iteration))
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
col1=np.asarray(col1)
col2=np.asarray(col2)
#print(n)
#Cited work from Travis Mitchell - put reference!
#def FI(D,x,t) :
#    #think of a way to include the range of t and create list of f for each t
#    fxt = []
#        
#    ff=np.exp(-np.array(x)**2/4.0/D/t)/np.sqrt(4.0*D*np.pi*t)
#    return(np.matrix(ff).transpose())  
#FI(D,col1,t)
def Explicit(n,xa,xb,f0,c,xlist,iteration,dt):
    counter = 0
    f1list=[]
    f0_a=f0
    while counter < iteration:
        Z = np.zeros(n) #zero matrix size 71x71
#        print(np.shape(Z))
        S = Z + (1 - 2*c)*np.eye(n) + c*np.eye(n,k=-1)+c*np.eye(n,k=1) #constructs S matrix without BC considered
        
#        print(S)
        #uses Dirichlet BC to finalise matrix S
#        S[0,:] = 0; S[0,0] = 1
        """Change left BC to Neumann's (1st order)"""
        S[0,:] = S[1,:]
        S[-1,:]= 0; S[-1,-1]=1
#        print(S)
        f0 = np.transpose(np.asarray(f0_a))
#        print(counter)
#        print(f0)
#        f1=np.matmul(S,f0)
        f1=np.dot(S,f0)
#        print(f1)
        f1list.append(f1)
        plt.plot(col1,f1,label=counter)
        f0_a=f1
#        print(f0_a)
        
#        print(counter)
        counter +=1
        
    plt.xlabel("x")
    plt.xlabel("f(x)")
#    plt.legend()
    plt.title("Explicit")
    plt.show()
    plt.plot(col1,f1,label=counter)
    plt.xlabel("x")
    plt.xlabel("f(x)")
#   plt.legend()
    plt.title(str("Explicit - final iteration at c = "+str(c)))
    plt.show()
#    print(f1)
#    print(f0)
#    print(f1list)
    return f1,f1list

#def 1DPDE(n,xa,xb,bc,f0,c,iteration):
#    T = np.ones
#1DPDE(n,0,14,[0,0],col2,CFL[0],10)
def Implicit(n,xa,xb,f1,c,xlist,iteration,dt):
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
        f1_a=f0
    plt.xlabel("x")
    plt.xlabel("f(x)")
    plt.title("Implicit")
#    plt.legend()
    plt.show()
    plt.plot(col1,f0,label=counter)
    plt.xlabel("x")
    plt.xlabel("f(x)")
#   plt.legend()
    plt.title(str("Implicit - final iteration at c = "+str(c)))
    plt.show()
    return f0    

#E1=Explicit(n,0,14,col2,c,col1,iteration,dt)
#I1=Implicit(n,0,14,col2,c,col1,iteration,dt)
#E1
#I1
def Crank_Nicolson(n,xa,xb,f0,c,xlist,iteration,dt):
    counter = 0
    f0_a=f0
    while counter < iteration:
        Z = np.zeros(n) #zero matrix size 71x71
        H= Z + (1 + c)*np.eye(n) + (-c/2)*np.eye(n,k=-1)+(-c/2)*np.eye(n,k=1) #creates H matrix
        #Apply BCs: Left BC is Neumann's and right BC is Dirichlet's
        H[0,:] = 0; H[0,0] = -1; H[0,1]=1
        H[-1,:]= 0; H[-1,-1]=1
        
        S= Z + (1 - c)*np.eye(n) + (c/2)*np.eye(n,k=-1)+(c/2)*np.eye(n,k=1) #creates H matrix
        #Apply BCs: Left BC is Neumann's and right BC is Dirichlet's
#        S[0,:] = 0; S[0,0] = 1
        S[0,:] = 0
#        S[0,:] = S[1,:]
#        S[-1,:]= 0; S[-1,-1]=1
        S[-1,:] = 0
        f0 = np.transpose(np.asarray(f0_a))
#        print("H=",H)
#        print("S=",S)
#        print(inv(H))
        Hinv_S=np.dot(inv(H),S)
        f1=np.dot(Hinv_S,f0)
        plt.plot(col1,f1,label=counter)
        counter +=1
        f0_a=f1
    plt.xlabel("x")
    plt.xlabel("f(x)")
    plt.title("Crank-Nicolson")
#    plt.legend()
    plt.show()
    plt.plot(col1,f1,label=counter)
    plt.xlabel("x")
    plt.xlabel("f(x)")
#   plt.legend()
    plt.title(str("Crank-Nicolson - final iteration at c = "+str(c)))
    plt.show()
    return f0    
C1=Crank_Nicolson(n,0,14,col2,c,col1,iteration,-0.1)       

def trap_rule(f,a,b,n):

    """From below, trapezoid rule to funcion f within bounaries up to given harmonics, n is applied.
    Parameters:
        f: function to integrate
        a: lower bounds
        b: upper bounds
        n: mode of harmonics
    Output:
    res: result : result of applying trapezoid rule
    """
    h = (b-a)/n #h = (b-a)/n which b=10,a=0 and n=100
    s = 0.5 * (f(a) + f(b))
    for i in range(1, n+1, 1):
        s += + f(a + i * h)
    # print("result using numpy is ", np.trapz(xx1, 'x', b-a, 'x'))
    print (h*s)
    return h*s

#def sep_var(N,J,x,F_0,xb):
#    fjlist=[] #list containing f(x) at each x value (x values are elements of col1)
#    for i in col1:
#        fj=0
#        for n in range(len(col2)):
#            An=0
#            for j in range(len(col2)-1):
##                print("j at n= ",n," is ",j)
#                g0 = np.cos(col1[0]*(pi/28)*(1+2*0))
#                gn = np.cos(col1[-1]*(pi/28)*(1+2*71))
#                glist=[]
#                for i in col1:
#                    
#                    glist.append(np.cos(i*(pi/28)*(1+2*n)))
#                glist=np.asarray(glist)
#                An = ((np.sum(np.asarray(col2)*glist)) + (col2[0]*g0 + col2[-1]*gn)/2)*(2/71)
##                print(An)
#            fj += An*np.exp(-((pi*n/xb)**2)*D*2)*cos(pi*i*n/xb)
##        print("fj at x= ",i , "is ", fj)
#        fjlist.append(fj)
#    plt.plot(col1,fjlist,label="spectral")
#    plt.plot(col1,col2,label="actual")
#    plt.xlabel("x")
#    plt.ylabel("f(x)")
#    plt.legend()
#    plt.title("Spectral Solution")
#    plt.show()
#sep_var(71,71,col1,col2,14)
def g(n,x,t,B):
        g=np.cos(((pi+2*n*pi)/28)*x)
        return g
    
def f(n,x,t,B):
        f=((exp(-(((pi+(2*n*pi))/28)**2)*t*D))*(B[n])*np.cos(((pi+2*n*pi)/28)*x))
        return f


col3 = data_read("data2_investigation")[0]
for i in col3:
    col3[col3.index(i)]=float(i)
col4 = data_read("data2_investigation")[1]
for i in col4:
    col4[col4.index(i)]=float(i)
n=len(col3) #71
col3=np.asarray(col3)
col4=np.asarray(col4)

def SpectralSol(x,f0,t,D,N):

    f_01=f0[0]+x*(f0[-1]-f0[0])/x[-1]    
    f1=f0-f_01
    B=np.zeros(71)
    gm=np.zeros([71,71])
    Bnum=np.zeros(71)
    Bdenom=np.zeros(71)
    fm=np.zeros([70,71])
    fxj=np.zeros(71)

    for n in range(0,N-1):
        gm[n,:]=g(n,x,t,B)
        B[n]=(2/71)*((np.sum(f0[1:-2]*(gm[n,1:-2])))+(((f0[0]*gm[0,0])+(f0[70]*gm[70,70]))/2))
        fm[n,:]=f(n,x,t,B)    
    for n in range (0,N-1):
        fxj[n]=np.sum(fm[:,n])
           
    plt.figure(4)
    plt.plot(x, f0,"b", label='provided data' )
    plt.plot(x, fxj, "r", label='Solution data')
    plt.title(N)
    plt.show()
    return fxj
SpectralSol(col1,col2,2,0.5,71)

for i in range(0,70):
    SpectralSol(col3,col4,-6,0.5,i)


