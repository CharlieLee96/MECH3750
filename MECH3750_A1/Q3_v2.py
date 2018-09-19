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
#   f(x) is an even function
xrange=np.linspace(-pi,pi,10)
#count = np.linspace(0,100,100)
#pi = pi
##x=sym.symbols('x')
#def f(x):
#    if -pi <= x < 0:
#        return pi+x
#    elif 0 <= x <= pi:
#        return pi-x 
#
##def cos(x):
##    return cos(x)
##def sin(x):
##    return sin(x)
#t_list=[]
#for x in xrange:
#    print(x)
#    fx = f(x)
#    fx_sum=0
#    for n in count:
#        if n == 0:
#            print(f(x))
#            a0=integrate.quad(f,a=-pi,b=pi,args=x)[0]*(1/pi)
#            
#            fx_sum += a0*cos(n*x)
#        else:
#            a=integrate.quad(f*cos(n*x),-pi,pi,args=x)*(1/pi)
#            asum=a*cos(n)
#            b=integrate.quad(f*sin(n*x),-pi,pi,args=x)*(1/pi)
#            fx_sum +=b*sin(x)
#    t_list.append(fx_sum)
#
#print(t_list)
from sympy import fourier_series, pi
from sympy.abc import x
from sympy.plotting import plot

s = fourier_series(pi+x, (x, -pi, pi))
print(type(s))
s.scale(2).truncate()
sp=plot(s)
#plt.show()
print(s)