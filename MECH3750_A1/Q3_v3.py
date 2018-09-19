# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 10:07:45 2018
@author: s4395102
"""
from math import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
#import sympy as sym
#   f(x) is an even function -> bn = 0
nn=100
xrange1=np.linspace(-pi+0.01,0,nn/2,endpoint=False)
xrange2=np.linspace(0,pi,nn/2)
print(len(xrange1)+len(xrange2))

xrange=np.linspace(-pi+0.01,pi,nn,endpoint=False)
print(len(xrange))
count = np.linspace(1,nn,nn)
print(len(count))
ylist=[]
for i in xrange:
    if -pi <= i < 0:
        ylist.append(pi+i)
    if 0 <= i <= pi:
        ylist.append(pi-i)

fxlist=[]
for x in xrange1:
    fx=0
    fx += 2*pi**2 #a0 term
    ansum=0
    for i in count:
#        print("an= a",i)
        a_n_1 = ((2*sin(i*pi))/i)+(2*(pi*i*sin(pi*i) + cos(pi*i) - 1))/(pi*(i**2))
        a_n_2=cos(x*i)
#        print("an at ",x,"is",an)
        ansum += a_n_1*a_n_2
#        print(a_n)
#    print("ansum is ",ansum)
    fx += ansum
    fxlist.append(fx)
for x in xrange2:
    fx=0
    fx += 2*pi**2 #a0 term
    ansum=0
    for i in count:
#        print("an= a",i)
        a_n_1 = ((2*sin(i*pi))/i)-(2*(pi*i*sin(pi*i) + cos(pi*i) - 1))/(pi*(i**2))
        a_n_2=cos(x*i)
        an = a_n_1*a_n_2
#        print("an at ",x,"is",an)
        ansum += a_n_1*a_n_2
#        print(a_n)
#    print("ansum is ",ansum)
    fx += ansum
    fxlist.append(fx)
#print(fxlist)

#plt.plot(xrange, ylist)
#plt.show()


plt.plot(xrange,ylist,'r',xrange,fxlist,'g')
plt.show()