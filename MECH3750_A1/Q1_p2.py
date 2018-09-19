# -*- coding: utf-8 -*-
"""
Created on Sat Sep  1 22:32:56 2018

@author: s4395102
"""

#from math import *
import math as m
import scipy as sci
import numpy as np
#import plot.matplotlib as plt
import sympy as sym
from sympy.physics.mechanics import dynamicsymbols
#--------------------uncomment below for step 6-------------------------
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

x1,x2,y1,y2 =sym.symbols('x1 y1 x2 y2') #declares t u v w as symbols for differentiation
#t,u,v,w    
sym.init_printing(use_unicode=True)

var=['x1','x2','y1','y2']
"""initial guess"""
#x0 = sym.Matrix([1,1,1,1]) 
"""uses functions as itself of symbols rather than using given numerical values"""
PE = m1*g*(y1-y0) + m2 * g * (y2-y0)

#length of elongation of each spring is named d1, d2, d3
d01 = ((x1-x0)**2 + (y1-y0)**2)**0.5- L01
d12 = ((x2-x1)**2 + (y2-y1)**2)**0.5- L12
d23 = ((x3-x2)**2 + (y3-y2)**2)**0.5- L23
Ee = 0.5*(k01*d01 + k12*d12 + k23*d23)
f=PE + Ee
#print(f)

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
#    x1=1.1
#    x2=0,9
#    y1=1.7
#    y2=0.8
    for i in l:
#        print(i)
        a=sym.diff(f,i)
        part_der.append(a)
    return part_der
a = part_der(f,var) #derives partial derivatives of total energy    
#print(a[0])
#print(a[1])
#print(a[2])
#print(a[3])
#x1=1.1
#x2=0,9
#y1=1.7
#y2=0.8
b={'x1':1.1,'y1':0.9,'x2':1.7,'y2':0.8}
#print(b)
#print(a[0].subs(b)) #works but gives wrong value
#print(a[1].subs(b))
#print(a[2].subs(b))
#print(a[3].subs(b))

def forward_diff(f,l):
    
    der=[] #matrix of forward difference values
    for i in l:
        i_val = float(l[i]) # i is the value. i_val is the value of i.        
        test_value = {'x1':1,'y1':1,'x2':1,'y2':1} #matrix with added h ()
        test_value[i] = i_val
#        print(test_value)
#        print(type(test_value))
        init_value = {'x1':1,'y1':1,'x2':1,'y2':1} #sets original f with all unknowns as 1
        f_orig  = f.subs(init_value)    #function with all values as 1
#        print(type(f_orig))
        df      = f.subs(test_value)    #function with modified unknown values
#        print(type(df))
        tot = (df-f_orig)/(i_val-1)
        der.append(tot)
    print(der)
    return der

forward_diff(f,b)
    