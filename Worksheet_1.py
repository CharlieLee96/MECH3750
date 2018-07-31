# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import matplotlib.pyplot as plt



def nfactorial(n):
    a = 1
    b = n
#    print("b= ",b)
    while b >= 1:
        a = a * b
        b += -1
#        print("b= ", b)
#    print(n,"! = ", a)
    return float(a)
    
def taylor(a,n,x):
    """a = evaluation point
    n = order of evaluation"""
    xlist = []
    ylist = []
    tot = 0
    b = 0
    xlist.append(b)
    b += 1
    tot=np.exp(a)
    ylist.append(tot)
    while b <= n:
        tot += (np.exp(a)*(x-a)**b)/nfactorial(b)
        ylist.append(tot)
        xlist.append(b)
        b += 1
    
    print(xlist)
    print(ylist)
    print("f(x) = ",tot)
    plt.plot(xlist,ylist,'r')
    plt.show()
    return [tot,ylist] #returns a list containing the final value, and the list of values up to the error order
def taylor_ln(a,n,x):
    xlist = []
    ylist = []
    tot = 0
    b = 0
    xlist.append(b)
    b += 1
    tot=np.log(a)
    ylist.append(tot)
    while b <= n:
        tot += ((x-a)**b)/((a**b)*b)
        ylist.append(tot)
        xlist.append(b)
        b += 1
    
    print(xlist)
    print(ylist)
    print("f(x) = ",tot)
    plt.plot(xlist,ylist,'r')
    plt.show()
    return [tot, ylist]


#names = ['bar', 'chocolate', 'chips']
#weights = [0.05, 0.1, 0.25]
#costs = [2.0, 5.0, 3.0]
#unit_costs = [40.0, 50.0, 12.0]
#
#
#titles = ['names', 'weights', 'costs', 'unit_costs']
#data = [titles] + list(zip(names, weights, costs, unit_costs))
#
#for i, d in enumerate(data):
#    line = '|'.join(str(x).ljust(12) for x in d)
#    print(line)
#    if i == 0:
#        print('-' * len(line))
