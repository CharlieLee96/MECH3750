#MECH3750: Tutorial Sheet 2
# Week 2 problem
# Updated: 06.08.18
# Maintainer: @TravisMitchell

from math import sqrt
from numpy import zeros,linalg,dot,transpose,linspace
from matplotlib.pyplot import plot,xlabel,ylabel,show

# data
U = [91.08, 85.04, 78.54, 68.51, 58.35, 49.63, 41.45, 34.32]
E = [0.8560, 0.8414, 0.8375, 0.8160, 0.7977, 0.7762, 0.7535, 0.7347]

m = len(U)
A = zeros((m,3), float) #Design Matrix
b = zeros((m), float)   
for k in range(m):
    A[k,0] = 1.0
    A[k,1] = U[k]
    A[k,2] = sqrt(U[k])
    b[k]   = E[k]

alpha = linalg.solve(dot(transpose(A),A), dot(transpose(A),b))
print("A = ",A)
print("b^T = ",b)
print("alpha = ",alpha)

plot(U,E, 'xr')
x_plot = linspace(34,92, 100)
y_plot = [alpha[0] +alpha[1]*x+ alpha[2]*sqrt(x) for x in x_plot]
plot(x_plot, y_plot, 'b:')
xlabel("U (m/s)")
ylabel("E (m/s)")
show()

