# -*- coding: utf-8 -*-
"""
Created on Sat Sep  1 22:32:56 2018

@author: s4395102
"""

#from math import *
import math as m
import scipy as sci
import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
from sympy.physics.mechanics import dynamicsymbols
#-------------------- step 5-------------------------
x0=0.
y0=1.
x3=3.
y3=1.
g=9.8
#case 1
#m1=1
#m2=1
#L01=0.9
#L12=0.9
#L23=0.9
#k01=100
#k12=0.1
#k23=100
#case 2
#m1=0.01
#m2=0.01
#L01=0.2
#L12=0.2
#L23=0.2
#k01=100
#k12=100
#k23=100
#-------------------- step 6-------------------------
x0=0.
y0=1.
x3=3.
y3=1.
g=9.8
m1=1
m2=3
L01=0.9
L12=0.8
L23=0.7
k01=90
k12=100
k23=80

x1,y1,x2,y2 =sym.symbols('x1 y1 x2 y2') #declares t u v w as symbols for differentiation
#t,u,v,w    
sym.init_printing(use_unicode=True)

var=['x1','y1','x2','y2']
"""initial guess"""
#x0 = sym.Matrix([1,1,1,1]) 
"""uses functions as itself of symbols rather than using given numerical values"""

#--------------method 1: simple substitute for forward difference-------------------------------------
#print("method 1 is using forward difference on simpler version of the code.")
f_t=m1*g*y1+m2*g*y2+(k01/2)*(((x1-x0)**2+(y1-y0)**2)**0.5-L01)**2+(k12/2)*(((x2-x1)**2+(y2-y1)**2)**0.5-L12)**2+(k23/2)*(((x3-x2)**2+(y3-y2)**2)**0.5-L23)**2
b={'x1':1.1,'y1':0.9,'x2':1.7,'y2':0.8}
fn=f_t.subs(b)
#print(fn)
dx=0.0001
b1={'x1':1.1+dx,'y1':0.9,'x2':1.7,'y2':0.8}
b2={'x1':1.1,'y1':0.9+dx,'x2':1.7,'y2':0.8}
b3={'x1':1.1,'y1':0.9,'x2':1.7+dx,'y2':0.8}
b4={'x1':1.1,'y1':0.9,'x2':1.7,'y2':0.8+dx}
fn1=(f_t.subs(b1)-fn)/dx
fn2=(f_t.subs(b2)-fn)/dx
fn3=(f_t.subs(b3)-fn)/dx
fn4=(f_t.subs(b4)-fn)/dx


#-----------------------method 2: using sympy to directly derive partial difference-----------
#length of elongation of each spring is named d1, d2, d3
PE = m1*g*(y1-y0) + m2 * g * (y2-y0)
d01 = (((x1-x0)**2 + (y1-y0)**2)**0.5- L01)**2
d12 = (((x2-x1)**2 + (y2-y1)**2)**0.5- L12)**2
d23 = (((x3-x2)**2 + (y3-y2)**2)**0.5- L23)**2
Ee = 0.5*(k01*d01 + k12*d12 + k23*d23)
f=PE + Ee
#print(type(f))
def part_der(f,l):
    
    """given an input function f and a list of variables l, returns the the list
       of partial derivative functions
       Input:
           f : function
           l: list of variables
       Output:
           list of partially derived functions
           """           
    part_der = []
    for i in l:
#        print(i)
        a=sym.diff(f,i)
        part_der.append(a)
    return part_der
a = part_der(f,var) #derives partial derivatives of total energy    

print(a[0])
print(a[1])
print(a[2])
print(a[3])
#x1=1.1
#x2=0,9
#y1=1.7
#y2=0.8
b={'x1':1.1,'y1':0.9,'x2':1.7,'y2':0.8}
#print(b)
#print(f.subs(b))
print("method 1: simple forward difference",fn1,fn2,fn3,fn4)
print("method 2: symbollic calculation using Sympy")
print("df/dx1 = ",a[0].subs(b))
print("df/dy1 = ",a[1].subs(b))
print("df/dx2 = ",a[2].subs(b))
print("df/dy2 = ",a[3].subs(b))
#--------------method 3: more elaborate forward difference-------------------------------------
def forward_diff(f,l,dx,iterations):
    f_orig  = f.subs(l)    #function with all values as 1
#    print("original f is ",f_orig)
    der1=[1, 1, 1, 1,] #matrix of forward difference values
    n=0
    while n < iterations+0.1:
        for i in l:
            i_val = float(l[i]) # i is the value. i_val is the value of i.        
#            print("for ",i," value is",i_val)
            test_value = b #matrix with added h ()
#            print("for iteration ",n," original test matrix is ",test_value)
            test_value[i] = i_val + dx
#            print("modified test matrix is ",test_value)
            df      = f.subs(test_value)    #function with modified unknown values
#            print(df)
            tot = (df-f_orig)/(dx)
            test_value[i] =test_value[i]-dx
#            print("at the end of iteration ",n,"variable values are: ",test_value)
            der1[list(l.keys()).index(i)] = tot
#            if n == 1.0:
#                print(der1)
            n += 0.25

    print("method 3: elaborate forward difference",der1)
    return der1
#b00={'x1':1,'y1':1,'x2':1.,'y2':1}
forward_diff(f,b,0.0001,8)
#----------------------step 7--------------------------------------
def newt(f1, f2, f3, f4, var,val):
    """takes first vector, 4 functions and a list specifying variables, runs newton's method and returns vectors of next iteration"""
    #define f(x)
#    for i in val:
#        print(type(i))
    x0=sym.Matrix([val[0],val[1],val[2],val[3]])
#    print(x0.shape)
    v={'x1':val[0],'y1':val[1],'x2':val[2],'y2':val[3]}
    fx = sym.Matrix([f1.subs(v),f2.subs(v),f3.subs(v),f4.subs(v)]) #f(x)
#    print(fx.shape)
    #define Jacobian matrix
    J_1=part_der(f1,var)
    J_1n=[]
    for i in J_1:
#        print(i)
        i_n=i.subs(v)
        J_1n.append(i_n)
#        print(i_n)
    J_2=part_der(f2,var)
    J_2n=[]
    for i in J_2:
        J_2n.append(i.subs(v))
    J_3=part_der(f3,var)
    J_3n=[]
    for i in J_3:
        J_3n.append(i.subs(v))
    J_4=part_der(f4,var)
    J_4n=[]
    for i in J_4:
        J_4n.append(i.subs(v))
    J =sym.Matrix([J_1,J_2,J_3,J_4]) #Jacobian matrix in symbollic form
    Jn=sym.Matrix([J_1n,J_2n,J_3n,J_4n])
#    print(Jn.shape)
#    print(Jn)
#    print(J)
    Jn_inv=Jn.inv() #inverse of Jacobian
    print(Jn_inv)
    x1=x0-Jn_inv*fx
    print(x1)
    print(x1.shape)
    return x1
#a=newt(a[0],a[1],a[2],a[3],var,[1.7,0.9,1.1,0.8])
    
#--------either uncomment line 168 or the lines below from 171 to 220.
def newt2(bn,val):
    x0=sym.Matrix([val[0],val[1],val[2],val[3]])
    v={'x1':val[0],'y1':val[1],'x2':val[2],'y2':val[3]}
    fx = sym.Matrix([a[0].subs(v),a[1].subs(v),a[2].subs(v),a[3].subs(v)]) #f(x)
    a0_dx1=part_der(a[0],var)[0]
    a0_dy1=part_der(a[0],var)[1]
    a0_dx2=part_der(a[0],var)[2]
    a0_dy2=part_der(a[0],var)[3]
    #print(part_der(a[0],var)[0].subs(bn))
    #print(a0_dx1)
    a0_dx1_n=a0_dx1.subs(bn)
    #print(a0_dx1_n)
    a0_dy1_n=a0_dy1.subs(bn)
    a0_dx2_n=a0_dx2.subs(bn)
    a0_dy2_n=a0_dy2.subs(bn)
    a0l=sym.Matrix([a0_dx1_n,a0_dy1_n,a0_dx2_n,a0_dy2_n])
    
    a1_dx1=part_der(a[1],var)[0]
    a1_dy1=part_der(a[1],var)[1]
    a1_dx2=part_der(a[1],var)[2]
    a1_dy2=part_der(a[1],var)[3]
    
    a1_dx1_n=a1_dx1.subs(bn)
    a1_dy1_n=a1_dy1.subs(bn)
    a1_dx2_n=a1_dx2.subs(bn)
    a1_dy2_n=a1_dy2.subs(bn)
    a1l=sym.Matrix([a1_dx1_n,a1_dy1_n,a1_dx2_n,a1_dy2_n])
    
    a2_dx1=part_der(a[2],var)[0]
    a2_dy1=part_der(a[2],var)[1]
    a2_dx2=part_der(a[2],var)[2]
    a2_dy2=part_der(a[2],var)[3]
    
    a2_dx1_n=a2_dx1.subs(bn)
    a2_dy1_n=a2_dy1.subs(bn)
    a2_dx2_n=a2_dx2.subs(bn)
    a2_dy2_n=a2_dy2.subs(bn)
    a2l=sym.Matrix([a2_dx1_n,a2_dy1_n,a2_dx2_n,a2_dy2_n])
    
    a3_dx1=part_der(a[3],var)[0]
    a3_dy1=part_der(a[3],var)[1]
    a3_dx2=part_der(a[3],var)[2]
    a3_dy2=part_der(a[3],var)[3]
    
    a3_dx1_n=a3_dx1.subs(bn)
    a3_dy1_n=a3_dy1.subs(bn)
    a3_dx2_n=a3_dx2.subs(bn)
    a3_dy2_n=a3_dy2.subs(bn)
    a3l=sym.Matrix([a3_dx1_n,a3_dy1_n,a3_dx2_n,a3_dy2_n])
    Jn=sym.Matrix([a0l,a1l,a2l,a3l])
    
    Jn2=sym.Matrix(4,4,Jn)
    print(Jn2)
    Jn_inv=Jn2.inv()
    
    
    x1=x0-Jn_inv*fx
#    print(x1)
#    print(x1[1,0])
    b1={'x1':x1[0,0],'y1':x1[1,0],'x2':x1[2,0],'y2':x1[3,0]}
#    print(b1)
    return [x1[0,0],x1[1,0],x1[2,0],x1[3,0],b1]
#    print(Jn2.shape)
b1=b
b1v=b1.values()
#print(b1v)
x1=newt2(b1,list(b1v))[4] #prints out x values after an iteration
#print(type(x1))
x2=newt2(x1,list(x1.values()))[4]
x3=newt2(x2,list(x2.values()))[4]
x4=newt2(x3,list(x3.values()))[4]
x5=newt2(x4,list(x4.values()))[4]
x6=newt2(x5,list(x5.values()))[4]
#print(x6)
x6v=list(x6.values()) #list with format: [x1,y1,x2,y2]
#print(x6v)
#---------plotting
x1=x6v[0]
#print(type(x1))
#print(x1)
y1=x6v[1]
#print(type(y1))
#print(y1)
x2=x6v[2]
#print(type(x2))
#print(x2)
y2=x6v[3]
#print(type(y2))
#print(y2)
x3=3.
xpos=np.array( [x0,x1,x2,x3] )
#print(type(xpos))
print("x positions are: ",xpos)
ypos=np.array( [y0,y1,y2,y3] )
print("y positions are: ",ypos)
plt.figure(1)
plt.ylim(-0.05, 1.05)
plt.xlabel("x position (m)")
plt.ylabel("y position (m)")
plt.xlim(-0.05, 3.05)
plt.plot(xpos,ypos,'o')
plt.show()