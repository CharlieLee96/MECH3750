# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 10:07:45 2018
@author: s4395102
"""
from math import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

xrange1=np.linspace(-pi+0.01,pi,(nn/10),endpoint=False)
xrange=np.linspace(-pi+0.01,pi,nn,endpoint=False)
count = np.linspace(1,nn,nn)
ylist=[]
for i in xrange1:
    if -pi <= i < 0:
        ylist.append(pi+i)
    if 0 <= i <= pi:
        ylist.append(pi-i)

fxlist=[]
for x in xrange:
    fx=0
    integrand1 = lambda x: (1/pi)*(pi-x)
    a_0=integrate.quad(integrand1,0,pi)[0]
    fx += a_0 #a0 term
    ansum=0
    for i in count:
        integrand=lambda x: -(2/pi)*(pi-x)*cos(i*x)
        a_n_1=integrate.quad(integrand,0,pi)[0]
        a_n_2=-cos(x*i)
        ansum += a_n_1*a_n_2
    fx += ansum
    fxlist.append(fx)


plt.plot(xrange,fxlist,'-')
plt.xlabel("x")
plt.xlabel("y")
plt.title("linear equation")
plt.show()
plt.plot(xrange1,ylist,'r')
plt.xlabel("x")
plt.xlabel("y")
plt.title("fourier representation")
plt.show()
plt.plot(xrange,fxlist,'-',xrange1,ylist,'r')
plt.xlabel("x")
plt.xlabel("y")
plt.title("linear eqn vs fourier series")
plt.show()