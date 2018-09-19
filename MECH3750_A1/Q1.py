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
t,u,v,w =sym.symbols('t u v w') #declares t u v w as symbols for differentiation
#t,u,v,w    
sym.init_printing(use_unicode=True)

var=['t','u','v','w']
#z = [1.1,2,0.9,1] #initial guess
#t=z[0]
#u=z[1]
#v=z[2]
#w=z[3]
vec = [0,0,0,0] #sets up vector
vec[0]=t**4+u**4-1
vec[1]=t**2-u**2+1
vec[2]=v**4+w**4-1
vec[3]=v**2-w**2+1
#print(vec[0],vec[1],vec[2],vec[3])

"""uses functions as itself of symbols rather than using given numerical values"""
f1=t**4+u**4-1
f2=t**2-u**2+1
f3=v**4+w**4-1
f4=v**2-w**2+1

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

#df1 = part_der (f1,var) #
#df2 = part_der (f2,var)
#df3 = part_der (f3,var)
#df4 = part_der (f4,var)
df1 = part_der (vec[0],var) #
df2 = part_der (vec[1],var)
df3 = part_der (vec[2],var)
df4 = part_der (vec[3],var)
#for i in df1:
#    print(i.subs(t,5.))
J = sym.Matrix([df1,df2,df3,df4]) #jacobian matrix
print("jacobian is:  ",J)

x0 = [1,1,1,1] #initial guess of x0
#def value_apply(func,sym, num):
#    """applies numerical values to a function"""
#    for sym,num in enumerate(num):
#        sym =num
def sub_num(func,t,u,v,w):
    """subs numerical values of t, u, v and w onto each function in a row in outputted jacobian matrix."""
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
    
    print(mat)
    return mat
sub_num(df1,x0[0],x0[1],x0[2],x0[3])
            
        
    
    
    
    
    