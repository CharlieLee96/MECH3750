# -*- coding: utf-8 -*-
"""
Spyder Editor

Author: Chul Lee
"""
import scipy as sci
import numpy as np
#import plot.matplotlib as plt
import sympy as sym
from sympy.physics.mechanics import dynamicsymbols
"""Question 1"""
t,u,v,w =dynamicsymbols('t u v w') #declares t u v w as symbols for differentiation
t,u,v,w    
var=['t','u','v','w']
#z = [1,1,1,1] #initial guess
#t=z[0]
#u=z[1]
#v=z[2]
#w=z[3]

#vec = np.array(z) #sets up vector
#vec[0]=t**4+u**4-1
#vec[1]=t**2-u**2+1
#vec[2]=v**4+w**4-1
#vec[3]=v**2-w**2+1
f1=t**4+u**4-1
f2=t**2-u**2+1
f3=v**4+w**4-1
f4=v**2-w**2+1

print(sym.diff(f1,t))
def part_der(f,l):
    
    """given an input function f and a list of variables l, returns the the list
       of partial derivative functions
       Input:
           f : function
           l: list of variables
       Output:
           list of partially derived functions
           """
           
    der = [] #empty list which derived functions will be appended onto
    for i in l:
        
        der.append(sym.diff(f, lambda i: f))
        print(der)           
    return der

part_der (f1,var)
