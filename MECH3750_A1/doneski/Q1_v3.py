# -*- coding: utf-8 -*-
"""
Spyder Editor

Author: Chul Lee
"""
from math import *
import scipy as sci
import numpy as np
#import plot.matplotlib as plt
import sympy as sym
from sympy.physics.mechanics import dynamicsymbols
"""Question 1"""
#-------------------------Step 1------------------------------------
t,u,v,w =sym.symbols('t u v w') #declares t u v w as symbols for differentiation
#t,u,v,w    
sym.init_printing(use_unicode=True)

var=['t','u','v','w']
"""initial guess"""
x0 = sym.Matrix([1,1,1,1]) #initial guess of x0
"""uses functions as itself of symbols rather than using given numerical values"""
#f1=t**4+u**4-1
#f2=t**2-u**2+1
#f3=v**4+w**4-1
#f4=v**2-w**2+1
#---------------------------------------step 3--------------------------------
f1=2*(u-t)**2-4*(u-t)+v**2+3*w**2+6*w+2
f2=(u-t)**2+v**2-2*v+2*w**2-5
f3=3*(u-t)**2-12*u**2+v**2+3*w**2+8
f4=u**2-v**2+t
f=[]
f.append(f1)
f.append(f2)
f.append(f3)
f.append(f4)
def part_der(f,l):
    
    """given an input function f and a list of variables l, returns the the list
       of partial derivative functions
       Input:
           f : function
           l: list of variables
       Output:
           list of partially derived functions
           """           
#    print(f)
    der = [] #empty list which derived functions will be appended onto
    for i in l:
        a=sym.diff(f,i)
        der.append(a)
    return der




def sub_num_df(func,t,u,v,w):
    """subs numerical values of t, u, v and w onto each partially derived functions in a row in outputted jacobian matrix."""
    mat = [] #empty list, numerical values of each function will be added
    tn = float(t) #prevents confusion caused by multiple variables stemming from same name
    un = float(u)
    vn = float(v)
    wn = float(w)
    for i in func:
        if i == 0:
            mat.append(0)
        else:
            mat.append(i.subs({'t':tn,'u':un,'v':vn,'w':wn}))
    
#    print(mat)
    return mat

         


def sub_num_f(func,t,u,v,w):
    """subs numerical values of t, u, v and w onto each vectors in a row in outputted jacobian matrix."""
    mat = [] #empty list, numerical values of each function will be added
    tn = float(t) #prevents confusion caused by multiple variables stemming from same name
    un = float(u)
    vn = float(v)
    wn = float(w)

    f = func.subs({'t':tn,'u':un,'v':vn,'w':wn})
    mat.append(f)
    return mat

def newt(x_n, f1, f2, f3, f4, var):
    """takes first vector, 4 functions and a list specifying variables, runs newton's method and returns vectors of next iteration"""
    print("for vector ",x_n)
    df1 = part_der (f1,var) #
    df2 = part_der (f2,var)
    df3 = part_der (f3,var)
    df4 = part_der (f4,var)
    J = sym.Matrix([df1,df2,df3,df4]) #jacobian matrix
#    print("jacobian is:  ",J)
    df1n = sub_num_df(df1,x0[0],x_n[1],x_n[2],x_n[3])
    df2n = sub_num_df(df2,x_n[0],x_n[1],x_n[2],x_n[3])
    df3n = sub_num_df(df3,x_n[0],x_n[1],x_n[2],x_n[3])
    df4n = sub_num_df(df4,x_n[0],x_n[1],x_n[2],x_n[3])   
    Jn = sym.Matrix([df1n,df2n,df3n,df4n]) #Jacobian matrix with numerical value
    print("Jn =   ",Jn)    
    
    f1n = sub_num_f(f1,x_n[0],x_n[1],x_n[2],x_n[3])
    f2n = sub_num_f(f2,x_n[0],x_n[1],x_n[2],x_n[3])
    f3n = sub_num_f(f3,x_n[0],x_n[1],x_n[2],x_n[3])
    f4n = sub_num_f(f4,x_n[0],x_n[1],x_n[2],x_n[3])
       
    f_x_n = sym.Matrix([f1n,f2n,f3n,f4n])
#    print("f_x0 = ",f_x0)
    x1 = x_n - Jn.inv() *f_x_n
#    print("after newton's method, approximation is:    ",x1)    
#    print(type(x1))
    return x1

counter = 0
#iterations = int(input("number of iterations?"))
iterations=2
x_n1 = []
x_n1.append(newt(x0,f1,f2,f3,f4,var))
print(x_n1)
while counter < iterations-1:
    x_n1.append(newt(x_n1[-1],f1,f2,f3,f4,var))
#    print(x_n1[-1])100
    counter += 1
print(x_n1[-1])
iterations=0
#-------------------------Step 2------------------------------------


b={'t':1,'u':1,'v':1,'w':1}
def forward_diff(f,l,dx,iterations):
    f_orig  = f.subs(l)    #function with all values as 1
#    print("original f is ",f_orig)
    der1=list(l.values()) #matrix of forward difference values
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

    print("for ",f," : ",der1)
    return der1
def newtf(tp):
    print("with given ",tp)
    fd1=forward_diff(f1,tp,0.02,8)
    fd2=forward_diff(f2,tp,0.02,8)
    fd3=forward_diff(f3,tp,0.02,8)
    fd4=forward_diff(f4,tp,0.02,8)
    Jac=sym.Matrix([fd1,fd2,fd3,fd4]) #jacobian matrix based on forward diff
    Jac_inv=Jac.inv()
    x_1=sym.Matrix(list(b.values()))-Jac_inv*sym.Matrix([f1.subs(b),f2.subs(b),f3.subs(b),f4.subs(b)])
#    print(x_1)
#    return list(x_1[0,0],x_1[1,0],x_1[2,0],x_1[3,0])
    result={'t':x_1[0,0],'u':x_1[1,0],'v':x_1[2,0],'w':x_1[3,0]}
    print(result)
    return result
it1=newtf(b) #iteration 1
#print(it1)
#it2=newtf(it1)
#print(it2)
itera=100
n=0
bn=b
while n < itera+0.1:
    itn=newtf(bn)
    bn = itn
    n += 1
    
#---------------------------------------step 3--------------------------------
#f5=2*(u-t)**2-4*(u-t)+v**2+3*w**2+6*w+2
#f6=(u-t)**2+v**2-2*v+2*w**2-5
#f7=3*(u-t)**2-12*u**2+v**2+3*w**2+8
#f8=u**2-v**2+t
#x_n2 = []
#x_n2.append(newt(x0,f5,f6,f7,f8,var))
#while counter < iterations-1:
#    print("step 3")
#    x_n2.append(newt(x_n2[-1],f5,f6,f7,f8,var))
##    print(x_n1[-1])100
#    counter += 1
#print(x_n2[-1])
#print(type(x_n2[-1]))


