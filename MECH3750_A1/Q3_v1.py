# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 10:07:45 2018
@author: s4395102
"""
from math import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
#   f(x) is an even function
xrange=np.linspace(-pi,pi,100)
count = np.linspace(0,100,100)
t_list=[]
for x in xrange:
    if -pi <x < 0:
        fx=pi+x
    elif 0 <= x < pi:
        fx=pi-x
        for n in count:
            fx_sum=0
            if n == 0:
                fx=fx*(2/pi)
                a=integrate.quad(fx*cos(n*x),-pi,pi,args=x)
                fx_sum += a*cos(x)
            else:
                fx=fx*(1/pi)
                a=integrate.quad(fx,-pi,pi,args=x)
                asum=a*cos(i)
                b=integrate.quad(fx*sin(n*x),-pi,pi,args=x)
                fx_sum +=b*sin(x)
            t_list.append(fx_sum)

print(t_list)