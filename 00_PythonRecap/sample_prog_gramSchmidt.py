# Author: Pr. Tony Roberts
#	Git Maintain: @TravisMitchell
#	    24-07-2018: + Created on git

#	Date: 24-07-2018

# TASK:
#	Apply the Gram-Schmidt procedure to create an
#	orthogonal matrix.


import numpy as np 

#Based on Gram-Schmidt Wiki (also see MATH2001)
def proj(v,u):
    answer=np.dot(u,v)/np.dot(u,u)*u
    return answer


v1=np.array([ 1,1,0,0,0,1 ])
v2=np.array([ 1,0,1,0,0,1 ])
v3=np.array([ 1,1,1,0,0,1 ])
v4=np.array([ 1,0,1,1,0,1 ])
v5=np.array([ 1,1,1,1,1,1 ])
v6=np.array([ 0,1,1,0,1,1 ])

u1=v1
u2=v2-proj(v2,u1)
u3=v3-proj(v3,u1)-proj(v3,u2)
u4=v4-proj(v4,u3)-proj(v4,u2)-proj(v4,u1)
u5=v5-proj(v5,u4)-proj(v5,u3)-proj(v5,u2)-proj(v5,u1)
u6=v6-proj(v6,u5)-proj(v6,u4)-proj(v6,u3)-proj(v6,u2)-proj(v6,u1)

Q=np.vstack([u1, u2, u3, u4, u5])
QT=np.transpose(Q)

Check1=np.matmul(Q,QT)
np.set_printoptions(precision=3)
print(Check1)

