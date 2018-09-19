#MECH3750: Assignment 1
# Read in the data file and plot

from numpy import loadtxt
import matplotlib.pyplot as plt

data = loadtxt('data2.txt')
plt.figure(1)
plt.plot(data,'-')
plt.show()
