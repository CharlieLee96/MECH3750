# Author: Pr. Tony Roberts
#	Git Maintain: @TravisMitchell
#	    24-07-2018: + Created on git

#	Date: 24-07-2018

# TASK:
#	Use the numpy package to find the current 
# 	through an array of resistors.
#
#	    2 > 5
#	    ^	v
#	1 > 3 - 6 > 8	
#	    v   ^
#           4 > 7


import numpy
A=numpy.matrix( ' -2  1  0  1  0  0 ; 	\
                   1 -4  1  0  1  0 ; 	\
                   0  1 -2  0  0  1 ; 	\
                   1  0  0 -2  1  0 ; 	\
                   0  1  0  1 -4  1 ; 	\
                   0  0  1  0  1 -2 ' )

b=numpy.matrix(' 0 ; -5 ; 0 ; 0 ; -2 ; 0')

Ainv=numpy.linalg.inv(A);
print(Ainv)

x = numpy.dot(numpy.linalg.inv(A), b)
print(x)


phi3=x[1]
phi6=x[4]
Currentin=5-phi3
Currentout=phi6-2
print('Check Current in from node 1 is identical to current out a node 8')
print('I  in  5-phi3=%.4f') % (Currentin)
print('I out  phi6-2=%.4f') % (Currentout)

